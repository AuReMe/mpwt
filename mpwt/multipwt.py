#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From genbank/gff files this script will create Pathway-Tools input data, then run Pathway-Tools's PathoLogic on them. It can also generate dat files.
The script takes a folder name as argument.

usage:
    mpwt -f=DIR [-o=DIR] [--patho] [--dat] [--md] [--cpu=INT] [-r] [-v] [--clean] [--log=FOLDER]
    mpwt --dat [-f=DIR] [-o=DIR] [--md] [--cpu=INT] [-v]
    mpwt -o=DIR [--md] [--cpu=INT] [-v]
    mpwt --clean [--cpu=INT] [-v]
    mpwt --delete=STR [--cpu=INT]
    mpwt --list

options:
    -h --help     Show help.
    -f=DIR     Working folder containing sub-folders with Genbank file.
    -o=DIR    Output folder path. Will create a output folder in this folder.
    --patho    Will run an inference of Pathologic on the input files.
    --dat    Will create BioPAX/attribute-value dat files from PGDB.
    --md    Move only the dat files into the output folder.
    --clean    Clean ptools-local folder, before any other operations.
    --delete=STR    Give a PGDB name and it will delete it (if multiple separe them with a ',', example: ecolicyc,athalianacyc).
    -r    Will delete files in ptools-local to reduce size of results when moving files to output folder (use it with -o).
    --cpu=INT     Number of cpu to use for the multiprocessing.
    --log=FOLDER     Create PathoLogic log files inside the given folder (use it with --patho).
    --list     List all PGDBs inside the ptools-local folder.
    -v     Verbose.

"""

import csv
import datetime
import docopt
import getpass
import os
import shutil
import subprocess
import sys

from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from gffutils.iterators import DataIterator

from mpwt import utils


def check_input_and_existing_pgdb(run_ids, input_folder, output_folder, verbose=None):
    """
    Check output folder for already existing PGDB, don't create them.
    Check if PGDBs are already in ptools-local folder.
    """
    # Check if there are files/folders inside the input folder.
    species_folders = [species_folder for species_folder in os.listdir(input_folder)]
    if len(species_folders) == 0:
        print("No folder containing genbank/gff file. In " + input_folder + " you must have sub-folders containing Genbank/GFF file.")
        return None

    # Check the structure of the input folder.
    invalid_characters = ['.', '/']
    for species_folder in species_folders:
        if os.path.isfile(input_folder+'/'+species_folder):
            print('Error: file inside the input_folder ({0}) instead of a subfolder. Check that you have a structure file of input_folder/species_1/species1.gbk and not input_folder/species_1.gbk.'.format(input_folder+'/'+species_folder))
            return None
        elif os.path.isdir(input_folder+'/'+species_folder):
            if any(char in invalid_characters for char in species_folder):
                print('Error: . or / in genbank/gff name {0} \nGenbank name is used as an ID in Pathway-Tools and Pathway-Tools does not create PGDB with . in ID.'.format(species_folder))
                return None

    if output_folder:
        already_present_outputs = [output_pgdb for output_pgdb in os.listdir(output_folder)]
        new_run_ids = set(run_ids) - set(already_present_outputs)
        new_run_ids = list(new_run_ids)

        if len(new_run_ids) == 0:
            print("All PGDBs are already present in the output folder. Remove them if you want a new inference.")
            return None

    else:
        new_run_ids = []
        for species_folder in species_folders:
            new_run_ids.append(species_folder)

    ptools_local_path = utils.find_ptools_path() + '/pgdbs/user/'
    already_present_pgdbs = [pgdb_species_folder[:-3] for pgdb_species_folder in os.listdir(ptools_local_path) if 'cyc' in pgdb_species_folder]
    if already_present_pgdbs != []:
        lower_case_new_run_ids = list(map(lambda x:x.lower(), new_run_ids))
        for pgdb in already_present_pgdbs:
            if pgdb in lower_case_new_run_ids:
                print("! PGDB {0} already in ptools-local, no inference will be launch on this species.".format(pgdb))
        lower_run_ids = dict(zip(lower_case_new_run_ids, new_run_ids))
        wo_ptools_run_ids = set(lower_case_new_run_ids) - set(already_present_pgdbs)
        new_run_ids = [lower_run_ids[run_id] for run_id in wo_ptools_run_ids]
        if len(new_run_ids) == 0:
            print("All PGDBs are already present in ptools-local. Remove them if you want a new inference, or use mpwt --dat only to create dat files from them.")
            return None

    return new_run_ids


def create_dat_creation_script(pgdb_id, lisp_pathname):
    """
    Input: a PGDB ID and an output path.
    Create a lisp script allowing dat extraction.
    """
    with open(lisp_pathname, 'w') as lisp_file:
        lisp_file.write("(in-package :ecocyc)")
        lisp_file.write('\n')
        lisp_file.write("(select-organism :org-id '" + pgdb_id + ")")
        lisp_file.write('\n')
        lisp_file.write('(let ((*progress-noter-enabled?* NIL))')
        lisp_file.write('\n')
        lisp_file.write("        (create-flat-files-for-current-kb))")


def create_dats_and_lisp(run_folder):
    """
    Read Genbank file and create Pathway Tools needed file.
    Create also a lisp file to create dat files from Pathway tools results.
    The name of the PGDB created by Pathway Tools will be the name of the species with '_' instead of space.

    Create organism-params.dat:
    ID  pgdb_id
    STORAGE FILE
    NCBI-TAXON-ID   taxon_id
    NAME    species_name

    Create genetic-elements.dats:
    NAME    
    ANNOT-FILE  gbk_name
    //

    Create dat_creation.lisp:
    (in-package :ecocyc)
    (select-organism :org-id 'pgdb_id)
    (create-flat-files-for-current-kb)
    """
    # Look for a Genbank file in the run folder.
    # PGDB ID corresponds to the name of the species folder.
    pgdb_id = run_folder.split('/')[-2]
    gbk_name = pgdb_id + ".gbk"
    gbk_pathname = run_folder + gbk_name
    gff_name = pgdb_id + ".gff"
    gff_pathname = run_folder + gff_name

    organism_dat = run_folder + 'organism-params.dat'
    genetic_dat = run_folder + 'genetic-elements.dat'

    taxon_id = ""
    species_name = ""

    if os.path.isfile(gbk_pathname):
        input_name = gbk_name
        # Take the species name and the taxon id from the genbank file.
        with open(gbk_pathname, "r") as gbk:
            # Take the first record of the genbank (first contig/chromosome) to retrieve the species name.
            first_seq_record = next(SeqIO.parse(gbk, "genbank"))
            try:
                species_name = first_seq_record.annotations['organism']
            except KeyError:
                raise KeyError('No organism in the Genbank. In the SOURCE you must have: ORGANISM  Species name')      

            # Take the source feature of the first record.
            # This feature contains the taxon ID in the db_xref qualifier.
            src_features = [feature for feature in first_seq_record.features if feature.type == "source"]
            try:
                for src_feature in src_features:
                    src_dbxref_qualifiers = src_feature.qualifiers['db_xref']
                    for src_dbxref_qualifier in src_dbxref_qualifiers:
                        if 'taxon:' in src_dbxref_qualifier:
                            taxon_id = src_dbxref_qualifier.replace('taxon:', '')
            except KeyError:
                raise KeyError('No taxon ID in the Genbank {0}. In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.'.format(gbk_pathname))

    elif os.path.isfile(gff_pathname):
        input_name = gff_name
        # Instead of parsing and creating a database from the GFF, parse the file and extract the first region feature.
        region_feature = [feature for feature in DataIterator(gff_pathname) if feature.featuretype == 'region'][0]
        if 'Dbxref' in region_feature.attributes:
            for dbxref in region_feature.attributes['Dbxref']:
                if 'taxon' in dbxref:
                    taxon_id = dbxref.replace('taxon:', '')
                else:
                    sys.exit('No taxon id in GFF file of {0}. GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))
        else:
            sys.exit('No Dbxref in GFF file of {0}. GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))
    else:
        sys.exit('Missing Genbank/GFF file. Check if you have a Genbank file and if it ends with .gbk or .gff')

    lisp_pathname = run_folder + "dat_creation.lisp"

    # Create the organism-params dat file.
    with open(organism_dat, 'w') as organism_file:
        organism_writer = csv.writer(organism_file, delimiter='\t', lineterminator='\n')
        organism_writer.writerow(['ID', pgdb_id])
        organism_writer.writerow(['STORAGE', "FILE"])
        organism_writer.writerow(['NCBI-TAXON-ID', taxon_id])
        organism_writer.writerow(['NAME', species_name])

    # Create the genetic-elements dat file.
    with open(genetic_dat, 'w') as genetic_file:
        genetic_writer = csv.writer(genetic_file, delimiter='\t', lineterminator='\n')
        genetic_writer.writerow(['NAME', ''])
        genetic_writer.writerow(['ANNOT-FILE', input_name])
        genetic_writer.writerow(['//'])

    # Create the lisp script.
    create_dat_creation_script(pgdb_id, lisp_pathname)


def pwt_input_files(run_folder):
    """
    Check if files needed by Pathway-Tools are available, if not create them.
    Check if there is a pathologic.log from a previous run. If yes, delete it.
    """
    required_files = set(['organism-params.dat','genetic-elements.dat','dat_creation.lisp'])
    files_in = set(next(os.walk(run_folder))[2])
    if global_verbose:
        species_folder = run_folder.split('/')[-2]

    if "pathologic.log" in files_in:
        os.remove(run_folder + "pathologic.log")

    missing_string = ""
    if required_files.issubset(files_in):
        if global_verbose:
            missing_string = "no missing files"
    else:
        if global_verbose:
            missing_string = "missing {0}".format("; ".join(required_files.difference(files_in))) + '. Inputs file created for {0}'.format(run_folder.split('/')[-2])
        create_dats_and_lisp(run_folder)
    if global_verbose:
        print("Checking inputs for {0}: {1}. ".format(species_folder, missing_string))


def create_lisp_script_PGDB():
    """
    Create a lisp script file for each PGDB in the ptools-local folder.
    Return a list containing all the path to the dat_creation.lisp.
    """
    ptools_local_path = utils.find_ptools_path()
    pgdb_folder = ptools_local_path + '/pgdbs/user/'
    tmp_folder = ptools_local_path + '/tmp/'

    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)

    lisp_folders = []
    for species_pgdb in os.listdir(pgdb_folder):
        if os.path.isdir(pgdb_folder + species_pgdb):
            pgdb_id = species_pgdb[:-3]
            pgdb_pathname = tmp_folder + pgdb_id + '/'
            os.mkdir(tmp_folder + pgdb_id)
            lisp_pathname = pgdb_pathname + "dat_creation.lisp"
            create_dat_creation_script(pgdb_id, lisp_pathname)
            lisp_folders.append(pgdb_pathname)

    return lisp_folders


def permission_change(folder_pathname):
    """
    Give permission to output files inside a folder.
    Used for log files and PGDB/dat files.
    """
    os.chmod(folder_pathname, 0o777)
    for root, subfolders, subfiles in os.walk(folder_pathname):
        for subfolder in subfolders:
            os.chmod(os.path.join(root, subfolder), 0o777)
        for subfile in subfiles:
            os.chmod(os.path.join(root, subfile), 0o777)


def check_pwt(genbank_paths, patho_log_folder):
    """
    Check PathoLogic's log.
    Create two log files (log_error.txt which contains Pathway-Tools log and resume_inference which contains summary of network).
    """
    if patho_log_folder:
        if os.path.exists(patho_log_folder) == False:
            print('No log directory, it will be created.')
            os.mkdir(patho_log_folder)

        patho_error_pathname = patho_log_folder + '/log_error.txt'
        patho_resume_pathname = patho_log_folder + '/resume_inference.tsv'

        patho_error_file = open(patho_error_pathname, 'w')
        patho_resume_file = open(patho_resume_pathname, 'w')
        patho_resume_writer = csv.writer(patho_resume_file, delimiter='\t', lineterminator='\n')
        patho_resume_writer.writerow(['species', 'gene_number', 'protein_number', 'pathway_number', 'reaction_number', 'compound_number'])

    failed_inferences = []
    passed_inferences = []

    for genbank_path in genbank_paths:
        species = genbank_path.split('/')[-2]
        patho_log = genbank_path + '/pathologic.log'

        if patho_log_folder:
            patho_error_file.write('------------ Species: ')
            patho_error_file.write(species)
            patho_error_file.write('\n')

        fatal_error_index = None

        if os.path.exists(patho_log):
            with open(patho_log, 'r') as input_file:
                for index, line in enumerate(input_file):
                    if 'fatal error' in line:
                        fatal_error_index = index
                        failed_inferences.append(species)
                        if patho_log_folder:
                            patho_error_file.write(line)
                            patho_resume_writer.writerow([species, 'ERROR', '', '', '', ''])

                    if fatal_error_index:
                        if index > fatal_error_index:
                            if patho_log_folder:
                                patho_error_file.write(line)

                    if 'Build done.' in  line:
                        if patho_log_folder:
                            patho_error_file.write(line)
                            resume_inference_line = next(input_file)
                            patho_error_file.write(resume_inference_line)
                            gene_number = int(resume_inference_line.split('PGDB contains ')[1].split(' genes')[0])
                            protein_number = int(resume_inference_line.split('genes, ')[1].split(' proteins')[0])
                            pathway_number = int(resume_inference_line.split('proteins, ')[1].split(' base pathways')[0])
                            reaction_number = int(resume_inference_line.split('base pathways, ')[1].split(' reactions')[0])
                            compound_number = int(resume_inference_line.split('reactions, ')[1].split(' compounds')[0])
                            patho_resume_writer.writerow([species, gene_number, protein_number, pathway_number, reaction_number, compound_number])

                        passed_inferences.append(species)
        else:
            if patho_log_folder:
                patho_error_file.write('No pathologic log, an error occured before PathoLogic run.\n')
                patho_resume_writer.writerow([species, 'ERROR', '', '', '', ''])
            print('No pathologic log for {0}, an error occured before PathoLogic run.'.format(species))

        if patho_log_folder:
            patho_error_file.write('------------\n\n')

    number_passed_inference = len(passed_inferences)
    number_failed_inference = len(failed_inferences)

    string_passed_build = 'build has' if number_passed_inference == 1 else 'builds have'
    string_failed_build = 'build has' if number_failed_inference == 1 else 'builds have'

    if number_passed_inference > 0:
        if global_verbose:
            print('\n{0} {1} passed!\n'.format(str(number_passed_inference), string_passed_build))
    if number_failed_inference > 0:
        if global_verbose:
            print('WARNING: {0} {1} failed! See the log for more information.\n'.format(str(number_failed_inference), string_failed_build))

    if patho_log_folder:
        patho_error_file.close()
        patho_resume_file.close()
        with open(patho_error_pathname,'r') as contents:
            save = contents.read()
        with open(patho_error_pathname, 'w') as output_file:
                output_file.write('Inference statistics:\n')
                if number_passed_inference > 0:
                    output_file.write('Build done: ' + str(number_passed_inference) + '\n')
                    output_file.write('Species: ' + ', '.join(passed_inferences) +  '\n\n')
                if number_failed_inference > 0:
                    output_file.write('Build failed: ' + str(number_failed_inference) + '\n')
                    output_file.write('Species: ' + ', '.join(failed_inferences) + '\n\n')
                output_file.write(save)

    if number_failed_inference > 0:
        sys.exit("Stop the inference.")

    if patho_log_folder:
        permission_change(patho_log_folder)


def extract_pgdb_pathname(run_folder):
    """
    Extract PGDB ID folder and path.
    Return a tuple with the name of the genbank and the pathname to the PGDB.
    """
    gbk_name = run_folder.split('/')[-2]

    ptools_local_path = utils.find_ptools_path()
    pgdb_path = ptools_local_path.replace('\n', '') + '/pgdbs/user/'

    # Replace all / by _ to ensure that there is no error with the path with gbk_name.
    pgdb_folder = pgdb_path + gbk_name.lower() + 'cyc/'

    return (gbk_name, pgdb_folder)


def pwt_error(genbank_path, subprocess_returncode, subprocess_stdout, subprocess_stderr):
    """
    Print error messages when there is a subprocess error during PathoLogic run.
    """
    print('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')
    species_name = genbank_path.split('/')[-2]
    print('Error for {0} with PathoLogic subprocess, return code: {1}'.format(species_name, str(subprocess_returncode)))
    if subprocess_stderr:
        print('An error occurred :' + subprocess_stderr.decode('utf-8'))

    print('\t', '=== Pathway-Tools log ===')
    for line in subprocess_stdout:
        print('\t', line, end='')
    print('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')


def run_pwt(genbank_path):
    """
    Create PGDB using files created during 'create_dats_and_lisp' ('organism-params.dat' and 'genetic-elements.dat').
    With verbose run check_output to retrieve the output of subprocess (and show when Pathway-Tools has been killed).
    Otherwise send the output to the null device.
    Command used:
    pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho
    """
    cmd_options = ['-no-web-cel-overview', '-no-cel-overview', '-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_pwt = ['pathway-tools', *cmd_options, '-patho', genbank_path]

    if global_verbose:
        print(' '.join(cmd_pwt))

    error_status = None
    errors = ['Restart actions (select using :continue):', 'Error']
    patho_lines = []

    try:
        patho_subprocess = subprocess.Popen(cmd_pwt, stdout=subprocess.PIPE, universal_newlines="")
        # Check internal error of Pathway-Tools (Error with Genbank).
        for patho_line in iter(patho_subprocess.stdout.readline, ""):
            patho_line = patho_line.decode('utf-8')
            if any(error in patho_line for error in errors):
                print('Error possibly with the genbank file.')
                error_status = True
                patho_subprocess.kill()

            patho_lines.append(patho_line)
            patho_subprocess.poll()
            return_code = patho_subprocess.returncode
            # Check if Pathway-Tools has been killed with returncode.
            # Also check if Pathway-Tools has finished PathoLogic inference (returncode 0).
            if return_code or return_code == 0:
                if return_code == 0:
                    patho_subprocess.stdout.close()
                    return error_status
                elif return_code != 0:
                    raise subprocess.CalledProcessError(return_code, cmd_pwt)

    except subprocess.CalledProcessError as subprocess_error:
        # Check error with subprocess (when process is killed).
        pwt_error(genbank_path, subprocess_error.returncode, patho_lines, patho_subprocess.stderr)
        error_status = True
    patho_subprocess.stdout.close()

    return error_status


def run_pwt_dat(genbank_path):
    """
    Create dat file using a lisp script created during 'create_dats_and_lisp'.
    Kill the subprocess when the command reach the Navigator Window opening proposition.
    If this proposition is not closed, the script can't continue.
    Command used:
    pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load
    """
    lisp_path = genbank_path + '/dat_creation.lisp'
    cmd_options = ['-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_dat = ['pathway-tools', *cmd_options, '-load', lisp_path]

    if global_verbose:
        print(' '.join(cmd_dat))

    error_status = None
    dat_creation_ends = ['Opening Navigator window.', 'No protein-coding genes with sequence data found.  Cannot continue.']
    load_errors = ['Error: Organism']
    load_lines = []

    try:
        load_subprocess = subprocess.Popen(cmd_dat, stdout=subprocess.PIPE, universal_newlines="")
        for load_line in iter(load_subprocess.stdout.readline, ""):
            load_line = load_line.decode('utf-8')
            if any(dat_end in load_line for dat_end in dat_creation_ends):
                load_subprocess.stdout.close()
                load_subprocess.kill()
                return
            if any(error in load_line for error in load_errors):
                error_status = True
                load_subprocess.kill()

            load_lines.append(load_line)
            load_subprocess.poll()
            return_code = load_subprocess.returncode
            if return_code:
                raise subprocess.CalledProcessError(return_code, cmd_dat)

    except subprocess.CalledProcessError as subprocess_error:
        # Check error with subprocess (when process is killed).
        pwt_error(genbank_path, subprocess_error.returncode, load_lines, load_subprocess.stderr)
        error_status = True
    load_subprocess.stdout.close()

    return error_status


def check_dat(pgdb_folder):
    """
    Check dats creation.
    """
    pgdb_folder_dbname = pgdb_folder[0].lower() + 'cyc'

    dats_path = pgdb_folder[1] +'/1.0/data/'

    dat_files = ["classes.dat", "compound-links.dat", "compounds.dat", "dnabindsites.dat", "enzrxns.dat", "gene-links.dat", "genes.dat", "pathway-links.dat",
                "pathways.dat", "promoters.dat", "protein-features.dat", "protein-links.dat", "proteins.dat", "protligandcplxes.dat", "pubs.dat",
                "reaction-links.dat", "reactions.dat", "regulation.dat", "regulons.dat", "rnas.dat", "species.dat", "terminators.dat", "transunits.dat"]

    dat_checks = []
    for dat_file in dat_files:
        dat_file_path = dats_path + '/' + dat_file
        if os.path.exists(dat_file_path):
            dat_checks.append(dat_file_path)
    if global_verbose:
        expected_dat_number = str(len(dat_files))
        found_dat_number = str(len(dat_checks))
        print('{0}: {1} on {2} dat files create.'.format(pgdb_folder_dbname, found_dat_number, expected_dat_number))


def run_move_pgdb(pgdb_folders):
    """
    Move the result files inside the shared folder containing the input data.
    pgdb_folder_dbname: ID of the species.
    pgdb_folder_path: path to the PGDB of the species (in ptools-local).
    """
    pgdb_folder_dbname = pgdb_folders[0]
    pgdb_folder_path = pgdb_folders[1]

    output_species = global_output_folder + '/' + pgdb_folder_dbname +'/'

    if global_dat_extraction:
        pgdb_folder_path = pgdb_folder_path + '/1.0/data'

    if global_size_reduction:
        for pgdb_file in os.listdir(pgdb_folder_path):
            file_to_move_pathname = pgdb_folder_path + '/' + pgdb_file
            output_file_pathname = output_species + pgdb_file
            if os.path.exists(output_species) == False:
                os.mkdir(output_species)
            if global_dat_extraction:
                if '.dat' in pgdb_file:
                    shutil.move(file_to_move_pathname, output_file_pathname)
            elif not global_dat_extraction:
                shutil.move(file_to_move_pathname, output_file_pathname)
        shutil.rmtree(pgdb_folder_path)
    else:
        shutil.copytree(pgdb_folder_path, output_species)
        if global_dat_extraction:
            for pgdb_file in os.listdir(output_species):
                if '.dat' not in pgdb_file:
                    os.remove(output_species+'/'+pgdb_file)


def multiprocess_pwt(input_folder=None, output_folder=None, patho_inference=None, dat_creation=None, dat_extraction=None, size_reduction=None, number_cpu=None, patho_log=None, verbose=None):
    """
    Function managing all the workflow (from the creatin of the input files to the results).
    Use it when you import mpwt in a script.
    """
    # Use a second verbose variable because a formal parameter can't be a global variable.
    # So if we want to use mpwt as a python import with this function we need to set a new global variable.
    # With this variable it is possible to set vervose in multiprocess function.
    global global_output_folder, global_dat_extraction, global_size_reduction, global_verbose
    global_output_folder = output_folder
    global_dat_extraction = dat_extraction
    global_size_reduction = size_reduction
    global_verbose = verbose

    # Use the number of cpu given by the user or all the cpu available.
    if number_cpu:
        number_cpu_to_use = int(number_cpu)
    else:
        number_cpu_to_use = cpu_count()
    mpwt_pool = Pool(processes=number_cpu_to_use)

    # Run folder contains sub-folders containing GBK/GFF file.
    if input_folder:
        run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
        if output_folder:
            if os.path.exists(output_folder) == False:
                if verbose:
                    print('No output directory, it will be created.')
                os.mkdir(output_folder)
        run_ids = check_input_and_existing_pgdb(run_ids, input_folder, output_folder, verbose)
        if not run_ids:
            return
        genbank_paths = [input_folder + "/" + run_id + "/" for run_id in run_ids]

        if verbose:
            print('~~~~~~~~~~Creation of input data from Genbank/GFF~~~~~~~~~~')
        mpwt_pool.map(pwt_input_files, genbank_paths)

        if patho_inference:
            if verbose:
                print('~~~~~~~~~~Inference on the data~~~~~~~~~~')
            error_status = mpwt_pool.map(run_pwt, genbank_paths)
            if verbose:
                print('~~~~~~~~~~Check inference~~~~~~~~~~')
            check_pwt(genbank_paths, patho_log)
            if any(error_status):
                sys.exit('Error during inference. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')

    # Create path for lisp if there is no folder given.
    if (dat_creation and not input_folder) or (output_folder and not input_folder):
        genbank_paths = create_lisp_script_PGDB()

    if verbose:
        print('~~~~~~~~~~Extraction of PGDB Pathname~~~~~~~~~~')
    pgdb_folders = {}
    for genbank_path in genbank_paths:
        pgdb_id_folder = extract_pgdb_pathname(genbank_path)
        pgdb_folders[genbank_path] = pgdb_id_folder

    if (input_folder and dat_creation) or dat_creation:
        if verbose:
            print('~~~~~~~~~~Creation of the .dat files~~~~~~~~~~')
        mpwt_pool.map(run_pwt_dat, genbank_paths)
        if verbose:
            print('~~~~~~~~~~Check .dat ~~~~~~~~~~')
        for genbank_path in pgdb_folders:
            check_dat(pgdb_folders[genbank_path])

    if (dat_creation and not input_folder) or (output_folder and not input_folder):
        ptools_local_path = utils.find_ptools_path()
        shutil.rmtree(ptools_local_path + '/tmp')

    if verbose:
        print('~~~~~~~~~~End of the Pathway-Tools Inference~~~~~~~~~~')
    if output_folder:
        if verbose:
            print('~~~~~~~~~~Moving result files~~~~~~~~~~')
        move_datas = []
        for genbank_path in pgdb_folders:
            move_datas.append(pgdb_folders[genbank_path])
        mpwt_pool.map(run_move_pgdb, move_datas)
        # Give access to the file for user outside the container.
        permission_change(output_folder)

    mpwt_pool.close()
    mpwt_pool.join()

    if verbose:
        print('~~~~~~~~~~The script have finished! Thank you for using it.')


def run_mpwt():
    """
    Function used with a mpwt call in the terminal.
    """
    args = docopt.docopt(__doc__)

    argument_number = len(sys.argv[1:])

    input_folder = args['-f']
    output_folder = args['-o']
    patho_inference = args['--patho']
    dat_creation = args['--dat']
    move_dat = args['--md']
    size_reduction = args['-r']
    number_cpu = args['--cpu']
    patho_log = args['--log']
    pgdb_to_deletes = args['--delete']
    pgdb_list = args['--list']
    verbose = args['-v']

    if pgdb_list:
        pgdbs = utils.list_pgdb()
        if pgdbs == []:
            print('No PGDB inside ptools-local.')
        else:
            print(str(len(pgdbs)) + ' PGDB inside ptools-local:\n' + '\t'.join(pgdbs))
        return

    # Delete PGDB if use of --delete argument.
    # Use a set to remove redudant PGDB.
    if pgdb_to_deletes:
        utils.remove_pgbds(list(set(pgdb_to_deletes.split(','))), number_cpu)
        return

    if args['--clean']:
        if verbose:
            print('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')
        utils.cleaning(number_cpu, verbose)
        if input_folder:
            utils.cleaning_input(input_folder, output_folder, verbose)
        if argument_number == 1 or (argument_number == 2 and verbose):
            sys.exit()

    multiprocess_pwt(input_folder, output_folder, patho_inference, dat_creation, move_dat, size_reduction, number_cpu, patho_log, verbose)


if __name__ == '__main__':
    run_mpwt()
