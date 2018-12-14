#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From genbank file this script will create Pathway-Tools input data, then run Pathway-Tools's PathoLogic on them and at last it will generate dat files for AuReMe.
The script takes a folder name as argument.
"""

import argparse
import csv
import datetime
import getpass
import os
import shutil
import subprocess
import sys

from Bio import SeqIO
from multiprocessing import Pool, cpu_count


def run():
    from mpwt.cleaning_pwt import cleaning, cleaning_input, delete_pgdb

    parser = argparse.ArgumentParser(usage="mpwt -f FOLDER [-o FOLDER] [-c INT] [-d] [-r ] [-v] [--clean] [--delete STRING]")
    parser.add_argument("-f", "--folder", dest = "folder", metavar = "FOLDER", help = "Folder containing sub-folders with Genbank file.")
    parser.add_argument("-o", "--output", dest = "output", metavar = "FOLDER", help = "Output folder path. Will create a output folder in this folder.", default=None)
    parser.add_argument("-c", "--cpu", dest = "cpu_number", metavar = "INT", help = "Number of cpu used by the multiprocess.", default=None)
    parser.add_argument("-d", "--dat", dest = "extract_dat", help = "Will extract only dat files from Pathway-Tools results.", action='store_true', default=None)
    parser.add_argument("-r", "--reduce", dest = "reduce_size", help = "Will delete files in ptools-local to reduce size of results.", action='store_true', default=None)
    parser.add_argument("-v", "--verbose", dest = "verbose", help = "mpwt will be more verbose.", action='store_true', default=None)
    parser.add_argument("--delete", dest = "delete", metavar = "STRING", help = "Give a PGDB name and it will delete it (if multiple separe them with a ',', example: ecolicyc,athalianacyc .", default=None)
    parser.add_argument("--clean", dest = "clean", help = "Arguments to clean ptools-local folder, before any other operations.", action='store_true', default=None)

    parser_args = parser.parse_args()

    # Print help and exit if no arguments.
    argument_number = len(sys.argv[1:])
    if argument_number == 0:
        parser.print_help()
        parser.exit()

    # Delete PGDB if use of --delete argument.
    pgdb_to_deletes = parser_args.delete
    if pgdb_to_deletes:
        for pgdb_to_delete in pgdb_to_deletes.split(','):
            delete_pgdb(pgdb_to_delete)
        return

    input_folder = parser_args.folder
    output_folder = parser_args.output
    dat_extraction = parser_args.extract_dat
    size_reduction = parser_args.reduce_size
    number_cpu = parser_args.cpu_number
    verbose = parser_args.verbose

    if parser_args.clean:
        if verbose:
            print('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')
        cleaning(verbose)
        if input_folder:
            cleaning_input(input_folder, output_folder, verbose)
        if argument_number == 1 or (argument_number == 2 and verbose):
            sys.exit()

    multiprocess_pwt(input_folder, output_folder, dat_extraction, size_reduction, number_cpu, verbose)


def multiprocess_pwt(input_folder, output_folder=None, dat_extraction=None, size_reduction=None, number_cpu=None, verbose=None):
    # Use a second verbose variable because a formal parameter can't be a global variable.
    # So if we want to use mpwt as a python import with this function we need to set a new global variable.
    # With this variable it is possible to set vervose in multiprocess function.
    global global_verbose
    global_verbose = verbose

    # Run folder contains sub-folders containing GBK file
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

    # Use the number of cpu entered by the user or all the cpu available.
    if number_cpu:
        number_cpu_to_use = int(number_cpu)
    else:
        number_cpu_to_use = cpu_count()
    p = Pool(processes=number_cpu_to_use)

    if verbose:
        print('~~~~~~~~~~Creation of input data from Genbank~~~~~~~~~~')
    for genbank_path in genbank_paths:
        pwt_run(genbank_path)

    if verbose:
        print('~~~~~~~~~~Inference on the data~~~~~~~~~~')
    p.map(run_pwt, genbank_paths)
    if verbose:
        print('~~~~~~~~~~Check inference~~~~~~~~~~')
    check_pwt(genbank_paths)

    if verbose:
        print('~~~~~~~~~~Extraction of PGDB Pathname~~~~~~~~~~')
    pgdb_folders = {}
    for genbank_path in genbank_paths:
        pgdb_id_folder = extract_pgdb_pathname(genbank_path)
        pgdb_folders[genbank_path] = pgdb_id_folder

    if dat_extraction:
        if verbose:
            print('~~~~~~~~~~Creation of the .dat files~~~~~~~~~~')
        p.map(run_pwt_dat, genbank_paths)
        if verbose:
            print('~~~~~~~~~~Check .dat ~~~~~~~~~~')
        for genbank_path in pgdb_folders:
            check_dat(pgdb_folders[genbank_path])

    if verbose:
        print('~~~~~~~~~~End of the Pathway-Tools Inference~~~~~~~~~~')
    if output_folder:
        if verbose:
            print('~~~~~~~~~~Moving result files~~~~~~~~~~')
        for genbank_path in pgdb_folders:
            move_pgdb(genbank_path, pgdb_folders[genbank_path], output_folder, dat_extraction, size_reduction)

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
        for pgdb in already_present_pgdbs:
            print("! PGDB {0} already in ptools-local, no inference will be launch on this species.".format(pgdb))
        lower_run_ids = dict(zip(map(lambda x:x.lower(),new_run_ids), new_run_ids))
        wo_ptools_run_ids = set(map(lambda x:x.lower(),new_run_ids)) - set(already_present_pgdbs)
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
    required_files = set(['organism-params.dat','genetic-elements.dat','script.lisp'])
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

    Create script.lisp:
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

    lisp_file = run_folder + "script.lisp"

    # Create the organism-params dat file.
    with open(organism_dat, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['ID', pgdb_id])
        writer.writerow(['STORAGE', "FILE"])
        writer.writerow(['NCBI-TAXON-ID', taxon_id])
        writer.writerow(['NAME', species_name])

    # Create the genetic-elements dat file.
    with open(genetic_dat, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['NAME', ''])
        writer.writerow(['ANNOT-FILE', gbk_name])
        writer.writerow(['//'])

    # Create the lisp script.
    with open(lisp_file, 'w') as file:
        file.write("(in-package :ecocyc)")
        file.write('\n')
        file.write("(select-organism :org-id '" + pgdb_id + ")")
        file.write('\n')
        file.write("(create-flat-files-for-current-kb)")

    print('Inputs file created for {0}.'.format(pgdb_id))


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


def pwt_error(genbank_path, subprocess_stdout, subprocess_stderr):
    if global_verbose:
        if subprocess_stderr != '':
            print('Error for {0}'.format(genbank_path))
            print('An error occurred :' + subprocess_stderr)

        errors = ['Error:', 'Killed', 'Exiting']
        index_error = None
        for index, line in enumerate(subprocess_stdout.split('\n')):
            if any(error in line for error in errors):
                index_error = index
        if index_error:
            print('----------------------------------------')
            print('Error for {0}'.format(genbank_path))
            print('\n'.join(subprocess_stdout.split('\n')[index_error:]))
            print('----------------------------------------')


def run_pwt(genbank_path):
    """
    Create PGDB using files created during 'create_dats_and_lisp'.
    With verbose run check_output to retrieve the output of subprocess (and show when Pathway-Tools has been killed).
    Otherwise send the output to the null device.
    """
    cmd_options = ['-no-web-cel-overview', '-no-cel-overview', '-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_pwt = ['pathway-tools', *cmd_options, '-patho', genbank_path]

    if global_verbose:
        print(' '.join(cmd_pwt))
    """
    patho_subprocess = subprocess.check_output(cmd_pwt)
    print(patho_subprocess.decode('utf-8'))
    patho_out, patho_err = patho_subprocess.communicate(input=b'(exit)')
    pwt_error(genbank_path, patho_out.decode("utf-8") , patho_err.decode("utf-8") )
    """
    try:
        subprocess.check_output(cmd_pwt)
    except subprocess.CalledProcessError as subprocess_error:
        print(subprocess_error.output)

def run_pwt_dat(genbank_path):
    """
    Create dat file using a lisp script created during 'create_dats_and_lisp'.
    Add an input to the subprocess call to close the Navigator Window opening proposition ('Enter name of X-window server to connect to (of the form HOST:N.M):').
    If this proposition is not closed the script can't continue.
    """
    lisp_path = genbank_path + '/script.lisp'
    cmd_options = ['-no-web-cel-overview', '-no-cel-overview', '-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_dat = ['pathway-tools', *cmd_options, '-load', lisp_path]

    if global_verbose:
        print(' '.join(cmd_dat))

    FNULL = open(os.devnull, 'w')

    load_subprocess = subprocess.Popen(cmd_dat, stdin=subprocess.PIPE, stdout=FNULL, stderr=subprocess.STDOUT)
    load_subprocess.communicate(input=b'(exit)')


def check_dat(pgdb_folder):
    """
    Check dats creation.
    Is it really useful?
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


def move_pgdb(genbank_path, pgdb_folder, output_folder, dat_extraction, size_reduction):
    """
    Move the result files inside the shared folder containing the input data.
    """
    pgdb_folder_dbname = pgdb_folder[0]
    pgdb_folder_path = pgdb_folder[1]

    output_species = output_folder + '/' + pgdb_folder_dbname +'/'

    if dat_extraction == True:
        pgdb_folder_path = pgdb_folder_path + '/1.0/data'

    if size_reduction is True:
        for pgdb_file in os.listdir(pgdb_folder_path):
            if os.path.exists(output_species) == False:
                os.mkdir(output_species)
            if dat_extraction == True:
                if '.dat' in pgdb_file:
                    shutil.move(pgdb_folder_path+'/'+pgdb_file, output_species+pgdb_file)
            elif not dat_extraction:
                shutil.move(pgdb_folder_path+'/'+pgdb_file, output_species+pgdb_file)
        shutil.rmtree(pgdb_folder_path)
    else:
        shutil.copytree(pgdb_folder_path, output_species)
        if dat_extraction == True:
            for pgdb_file in os.listdir(output_species):
                if '.dat' not in pgdb_file:
                    os.remove(output_species+'/'+pgdb_file)

    # Give access to the file for user outside the container.
    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_species])

    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_folder])


if __name__ == '__main__':
    run()
