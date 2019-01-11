#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From genbank file this script will create Pathway-Tools input data, then run Pathway-Tools's PathoLogic on them and at last it will generate dat files for AuReMe.
The script takes a folder name as argument.

usage:
    mpwt -f=DIR [-o=DIR] [--patho] [--dat] [--cpu=INT] [-r] [-v] [--clean]
    mpwt --dat [-f=DIR] [-o=DIR] [-v]
    mpwt --clean [-v]
    mpwt --delete=STR

options:
    -h --help     Show help.
    -f=DIR     Working folder containing sub-folders with Genbank file.
    -o=DIR    Output folder path. Will create a output folder in this folder.
    --patho    Will run an inference of Pathologic on the input files.
    --dat    Will extract only dat files from Pathway-Tools results.
    --clean    Clean ptools-local folder, before any other operations.
    --delete=STR    Give a PGDB name and it will delete it (if multiple separe them with a ',', example: ecolicyc,athalianacyc).
    -r    Will delete files in ptools-local to reduce size of results.
    --cpu=INT     Number of cpu to use for the multiprocessing.
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


def run_mpwt():
    """
    Function used with a mpwt call in the terminal.
    """
    from mpwt.cleaning_pwt import cleaning, cleaning_input, delete_pgdb

    args = docopt.docopt(__doc__)

    argument_number = len(sys.argv[1:])

    # Delete PGDB if use of --delete argument.
    pgdb_to_deletes = args['--delete']
    if pgdb_to_deletes:
        for pgdb_to_delete in pgdb_to_deletes.split(','):
            delete_pgdb(pgdb_to_delete)
        return

    input_folder = args['-f']
    output_folder = args['-o']
    patho_inference = args['--patho']
    dat_extraction = args['--dat']
    size_reduction = args['-r']
    number_cpu = args['--cpu']
    verbose = args['-v']

    if args['--clean']:
        if verbose:
            print('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')
        cleaning(verbose)
        if input_folder:
            cleaning_input(input_folder, output_folder, verbose)
        if argument_number == 1 or (argument_number == 2 and verbose):
            sys.exit()

    multiprocess_pwt(input_folder, output_folder, patho_inference, dat_extraction, size_reduction, number_cpu, verbose)


def multiprocess_pwt(input_folder=None, output_folder=None, patho_inference=None, dat_extraction=None, size_reduction=None, number_cpu=None, verbose=None):
    """
    Function managing all the workflow (from the creatin of the input files to the results).
    Use it when you import mpwt in a script.
    """
    # Use a second verbose variable because a formal parameter can't be a global variable.
    # So if we want to use mpwt as a python import with this function we need to set a new global variable.
    # With this variable it is possible to set vervose in multiprocess function.
    global global_verbose
    global_verbose = verbose

    # Use the number of cpu entered by the user or all the cpu available.
    if number_cpu:
        number_cpu_to_use = int(number_cpu)
    else:
        number_cpu_to_use = cpu_count()
    p = Pool(processes=number_cpu_to_use)

    # Run folder contains sub-folders containing GBK file
    if input_folder:
        run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
        if output_folder:
            if os.path.exists(output_folder) == False:
                if verbose:
                    print('No output directory, it will be created.')
                os.mkdir(output_folder)
        run_ids = check_existing_pgdb(run_ids, input_folder, output_folder)
        if not run_ids:
            return
        genbank_paths = [input_folder + "/" + run_id + "/" for run_id in run_ids]
        if len(genbank_paths) == 0:
            sys.exit("No folder containing genbank file. In " + input_folder + " you must have sub-folders containing Genbank file.")

        if verbose:
            print('~~~~~~~~~~Creation of input data from Genbank~~~~~~~~~~')
        for genbank_path in genbank_paths:
            pwt_run(genbank_path)

        if patho_inference:
            if verbose:
                print('~~~~~~~~~~Inference on the data~~~~~~~~~~')
            error_status = p.map(run_pwt, genbank_paths)
            if any(error_status):
                sys.exit('Error during inference. Process stopped.')
            if verbose:
                print('~~~~~~~~~~Check inference~~~~~~~~~~')
            check_pwt(genbank_paths)

    # Create path for lisp if there is no folder given.
    if dat_extraction and not input_folder:
        genbank_paths = create_lisp_script_PGDB()

    if verbose:
        print('~~~~~~~~~~Extraction of PGDB Pathname~~~~~~~~~~')
    pgdb_folders = {}
    for genbank_path in genbank_paths:
        pgdb_id_folder = extract_pgdb_pathname(genbank_path)
        pgdb_folders[genbank_path] = pgdb_id_folder

    if (input_folder and dat_extraction) or dat_extraction:
        if verbose:
            print('~~~~~~~~~~Creation of the .dat files~~~~~~~~~~')
        p.map(run_pwt_dat, genbank_paths)
        if verbose:
            print('~~~~~~~~~~Check .dat ~~~~~~~~~~')
        for genbank_path in pgdb_folders:
            check_dat(pgdb_folders[genbank_path])
        if dat_extraction and not input_folder:
            ptools_local_path = ptools_path()
            shutil.rmtree(ptools_local_path + '/tmp')

    if verbose:
        print('~~~~~~~~~~End of the Pathway-Tools Inference~~~~~~~~~~')
    if output_folder:
        if verbose:
            print('~~~~~~~~~~Moving result files~~~~~~~~~~')
        move_datas = []
        for genbank_path in pgdb_folders:
            move_datas.append([genbank_path, pgdb_folders[genbank_path], output_folder, dat_extraction, size_reduction])
        p.map(run_move, move_datas)

    if verbose:
        print('~~~~~~~~~~The script have finished! Thank you for using it.')


def check_existing_pgdb(run_ids, input_folder, output_folder):
    """
    Check output folder for already existing PGDB, don't create them.
    Check if PGDBs are already in ptools-local folder.
    """
    if output_folder:
        already_present_outputs = [output_pgdb for output_pgdb in os.listdir(output_folder)]
        new_run_ids = set(run_ids) - set(already_present_outputs)
        new_run_ids = list(new_run_ids)

    else:
        new_run_ids = []
        invalid_characters = ['.', '/']
        for species_folder in os.listdir(input_folder):
            if any(char in invalid_characters for char in species_folder):
                print('Error: . or / in genbank name {0} \nGenbank name is used as an ID in Pathway-Tools and Pathway-Tools does not create PGDB with . in ID.'.format(species_folder))
                return None
            new_run_ids.append(species_folder)

    ptools_local_path = ptools_path() + '/pgdbs/user/'
    already_present_pgdbs = [pgdb_species_folder[:-3] for pgdb_species_folder in os.listdir(ptools_local_path) if 'cyc' in pgdb_species_folder]
    if already_present_pgdbs != []:
        lower_case_new_run_ids = map(lambda x:x.lower(),new_run_ids)
        for pgdb in already_present_pgdbs:
            if pgdb in lower_case_new_run_ids:
                print("! PGDB {0} already in ptools-local, no inference will be launch on this species.".format(pgdb))
        lower_run_ids = dict(zip(lower_case_new_run_ids, new_run_ids))
        wo_ptools_run_ids = set(lower_case_new_run_ids) - set(already_present_pgdbs)
        new_run_ids = [lower_run_ids[run_id] for run_id in wo_ptools_run_ids]

    if len(new_run_ids) == 0:
        print("All PGDBs are already present in the output folder. Remove them if you want a new inference.")
        return None

    return new_run_ids


def pwt_run(run_folder):
    """
    Check if files needed by Pathway-Tools are available, if not create them.
    Check if there is a pathologic.log from a previous run. If yes, delete it.
    """
    required_files = set(['organism-params.dat','genetic-elements.dat','dat_extraction.lisp'])
    files_in = set(next(os.walk(run_folder))[2])
    if global_verbose:
        species_folder = run_folder.split('/')[-2]
        print("Checking Pathway-Tools species_folder inputs for {0}:".format(species_folder))

    if "pathologic.log" in files_in:
        os.remove(run_folder + "pathologic.log")

    if required_files.issubset(files_in):
        if global_verbose:
            print("OK")
    else:
        if global_verbose:
            print("%s missing" %"; ".join(required_files.difference(files_in)))
        create_dats_and_lisp(run_folder)


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

    Create dat_extraction.lisp:
    (in-package :ecocyc)
    (select-organism :org-id 'pgdb_id)
    (create-flat-files-for-current-kb)
    """
    # Look for a Genbank file in the run folder.
    # PGDB ID corresponds to the name of the species folder.
    pgdb_id = run_folder.split('/')[-2]
    gbk_pathname = run_folder + pgdb_id + ".gbk"
    gbk_name = pgdb_id + ".gbk"

    # Check if a Genbank file have been found.
    try:
        gbk_pathname
    except NameError:
        raise NameError('Missing Genbank file. Check if you have a Genbank file and if it ends with .gbk or .gb or .gbff')

    organism_dat = run_folder + 'organism-params.dat'
    genetic_dat = run_folder + 'genetic-elements.dat'

    taxon_id = ""
    species_name = ""

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
        src_feature = [feature for feature in first_seq_record.features if feature.type == "source"][0]
        try:
            taxon_id = src_feature.qualifiers['db_xref'][0].replace('taxon:', '')
        except KeyError:
            raise KeyError('No taxon ID in the Genbank. In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.')               

    lisp_pathname = run_folder + "dat_extraction.lisp"

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
        genetic_writer.writerow(['ANNOT-FILE', gbk_name])
        genetic_writer.writerow(['//'])

    # Create the lisp script.
    create_dat_extraction_script(pgdb_id, lisp_pathname)

    print('Inputs file created for {0}.'.format(pgdb_id))


def create_dat_extraction_script(pgdb_id, lisp_pathname):
    """
    Input: a PGDB ID and an output path.
    Create a lisp script allowing dat extraction.
    """
    with open(lisp_pathname, 'w') as lisp_file:
        lisp_file.write("(in-package :ecocyc)")
        lisp_file.write('\n')
        lisp_file.write("(select-organism :org-id '" + pgdb_id + ")")
        lisp_file.write('\n')
        lisp_file.write("(create-flat-files-for-current-kb)")


def create_lisp_script_PGDB():
    """
    Create a lisp script file for each PGDB in the ptools-local folder.
    Return a list containing all the path to the dat_extraction.lisp.
    """
    ptools_local_path = ptools_path()
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
            lisp_pathname = pgdb_pathname + "dat_extraction.lisp"
            create_dat_extraction_script(pgdb_id, lisp_pathname)
            lisp_folders.append(pgdb_pathname)

    return lisp_folders


def check_pwt(genbank_paths):
    """
    Check PathoLogic's log.
    Create two log files (log_error.txt which contains Pathway-Tools log and resume_inference which contains summary of network).
    """
    failed_inferences = []
    passed_inferences = []

    with open('log_error.txt', 'w') as output_file:
        with open('resume_inference.tsv', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
            writer.writerow(['species', 'gene_number', 'protein_number', 'pathway_number', 'reaction_number', 'compound_number'])
            for genbank_path in genbank_paths:
                species = genbank_path.split('/')[-2]
                patho_log = genbank_path + '/pathologic.log'

                output_file.write('------------ Species: ')
                output_file.write(species)
                output_file.write('\n')

                fatal_error_index = None

                with open(patho_log, 'r') as input_file:
                    for index, line in enumerate(input_file):
                        if 'fatal error' in line:
                            fatal_error_index = index
                            output_file.write(line)
                            writer.writerow([species, 'ERROR', '', '', '', ''])
                            failed_inferences.append(species)
                        if fatal_error_index:
                            if index > fatal_error_index:
                                output_file.write(line)
                        if 'Build done.' in  line:
                            output_file.write(line)
                            resume_inference_line = next(input_file)
                            output_file.write(resume_inference_line)
                            gene_number = int(resume_inference_line.split('PGDB contains ')[1].split(' genes')[0])
                            protein_number = int(resume_inference_line.split('genes, ')[1].split(' proteins')[0])
                            pathway_number = int(resume_inference_line.split('proteins, ')[1].split(' base pathways')[0])
                            reaction_number = int(resume_inference_line.split('base pathways, ')[1].split(' reactions')[0])
                            compound_number = int(resume_inference_line.split('reactions, ')[1].split(' compounds')[0])
                            writer.writerow([species, gene_number, protein_number, pathway_number, reaction_number, compound_number])

                            passed_inferences.append(species)

                output_file.write('------------\n\n')

    with open('log_error.txt','r') as contents:
        save = contents.read()
    with open('log_error.txt', 'w') as output_file:
            output_file.write('Inference statistics:\n')
            if len(passed_inferences) > 0:
                if global_verbose:
                    print('\n' + str(len(passed_inferences)) + ' builds have passed!\n')   
                output_file.write('Build done: ' + str(len(passed_inferences)) + '\n')
                output_file.write('Species: ' + ', '.join(passed_inferences) +  '\n\n')
            if len(failed_inferences) > 0:
                if global_verbose:
                    print('WARNING: ' + str(len(failed_inferences)) + ' builds have failed! See the log for more information.\n')
                output_file.write('Build failed: ' + str(len(failed_inferences)) + '\n')
                output_file.write('Species: ' + ', '.join(failed_inferences) + '\n\n')
            output_file.write(save)
    
    if len(failed_inferences) > 0:
        sys.exit("Stop the inference.")

    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', 'log_error.txt'])
    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', 'resume_inference.tsv'])


def ptools_path():
    """
    Find the path of ptools using Pathway-Tools file.
    """
    pathway_tools_str = subprocess.check_output('type pathway-tools', shell=True)
    pathway_tools_path = pathway_tools_str.decode('UTF-8').split('is ')[1].strip('\n')

    pathway_tools_file = open(pathway_tools_path, 'r')
    ptools_local_str = [line for line in pathway_tools_file if 'PTOOLS_LOCAL_PATH' in line][0]
    ptools_local_path = ptools_local_str.split(';')[0].split('=')[1].replace('"', '').strip(' ') + '/ptools-local'
    pathway_tools_file.close()

    return ptools_local_path


def extract_pgdb_pathname(run_folder):
    """
    Extract PGDB ID folder and path.
    """
    gbk_name = run_folder.split('/')[-2]

    ptools_local_path = ptools_path()
    pgdb_path = ptools_local_path.replace('\n', '') + '/pgdbs/user/'

    # Replace all / by _ to ensure that there is no error with the path with gbk_name.
    pgdb_folder = pgdb_path + gbk_name.lower() + 'cyc/'
    pgdb_id_folder = (gbk_name, pgdb_folder)

    return pgdb_id_folder


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
    pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -load
    """
    lisp_path = genbank_path + '/dat_extraction.lisp'
    cmd_options = ['-no-web-cel-overview', '-no-cel-overview', '-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_dat = ['pathway-tools', *cmd_options, '-load', lisp_path]

    if global_verbose:
        print(' '.join(cmd_dat))

    error_status = None
    dat_creation_ends = ['Opening Navigator window.']
    load_lines = []

    try:
        load_subprocess = subprocess.Popen(cmd_dat, stdout=subprocess.PIPE, universal_newlines="")
        for load_line in iter(load_subprocess.stdout.readline, ""):
            load_line = load_line.decode('utf-8')
            if any(dat_creation_end in load_line for dat_creation_end in dat_creation_ends):
                print('End of creation of dat file.')
                load_subprocess.stdout.close()
                load_subprocess.kill()
                return

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


def run_move(move_data):
    genbank_path = move_data[0]
    pgdb_folder = move_data[1]
    output_folder = move_data[2]
    dat_extraction = move_data[3]
    size_reduction = move_data[4]
    move_pgdb(genbank_path, pgdb_folder, output_folder, dat_extraction, size_reduction)


def move_pgdb(genbank_path, pgdb_folder, output_folder, dat_extraction, size_reduction):
    """
    Move the result files inside the shared folder containing the input data.
    """
    pgdb_folder_dbname = pgdb_folder[0]
    pgdb_folder_path = pgdb_folder[1]

    output_species = output_folder + '/' + pgdb_folder_dbname +'/'

    if dat_extraction:
        pgdb_folder_path = pgdb_folder_path + '/1.0/data'

    if size_reduction:
        for pgdb_file in os.listdir(pgdb_folder_path):
            if os.path.exists(output_species) == False:
                os.mkdir(output_species)
            if dat_extraction:
                if '.dat' in pgdb_file:
                    shutil.move(pgdb_folder_path+'/'+pgdb_file, output_species+pgdb_file)
            elif not dat_extraction:
                shutil.move(pgdb_folder_path+'/'+pgdb_file, output_species+pgdb_file)
        shutil.rmtree(pgdb_folder_path)
    else:
        shutil.copytree(pgdb_folder_path, output_species)
        if dat_extraction:
            for pgdb_file in os.listdir(output_species):
                if '.dat' not in pgdb_file:
                    os.remove(output_species+'/'+pgdb_file)

    # Give access to the file for user outside the container.
    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_species])

    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_folder])


if __name__ == '__main__':
    run_mpwt()
