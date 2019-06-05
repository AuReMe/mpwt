#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From genbank/gff files this script will create Pathway Tools input data, then run Pathway Tools's PathoLogic on them. It can also generate dat files.
The script takes a folder name as argument.

usage:
    mpwt -f=DIR [-o=DIR] [--patho] [--hf] [--dat] [--md] [--cpu=INT] [-r] [-v] [--clean] [--log=FOLDER]
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
    --hf    Use with --patho. Run the Hole Filler using Blast.
    --dat    Will create BioPAX/attribute-value dat files from PGDB.
    --md    Move only the dat files into the output folder.
    --clean    Clean ptools-local folder, before any other operations.
    --delete=STR    Give a PGDB name and it will delete it (if multiple separe them with a ',', example: ecolicyc,athalianacyc).
    -r    Will delete files in ptools-local to reduce size of results when moving files to output folder (use it with -o).
    --cpu=INT     Number of cpu to use for the multiprocessing (default=1).
    --log=FOLDER     Create PathoLogic log files inside the given folder (use it with --patho).
    --list     List all PGDBs inside the ptools-local folder.
    -v     Verbose.

"""

import csv
import datetime
import docopt
import getpass
import logging
import os
import shutil
import subprocess
import sys

from Bio import SeqIO
from multiprocessing import Pool
from gffutils.iterators import DataIterator

from mpwt import utils

logging.basicConfig(format='%(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def compare_input_ids_to_ptools_ids(compare_ids, ptools_run_ids, set_operation):
    """ Compare species IDs in input folder with IDs present in PGDB folder.

    Args:
        compare_ids (list): species IDs (folder and GBK/GFF file name)
        ptools_run_ids (list): PGDB IDs inside ptools-local
        set_operation (str): difference or intersection between the two inputs list
    Returns:
        list: ID intersecting or different between the two lists
    """
    # To compare them, lower case the species IDs, then run a difference or an intersection between set.
    lower_case_compare_ids = list(map(lambda x:x.lower(), compare_ids))

    lower_compare_ids = dict(zip(lower_case_compare_ids, compare_ids))

    # Difference to obtain all IDs which are not in PGDB folder (which need to be run on PathoLogic).
    if set_operation == 'difference':
        compare_ids_ptools = set(lower_case_compare_ids) - set(ptools_run_ids)

    # Intersection to obtain all IDs which are already in PGDB folder.
    # Which need only to create BioPAX/dat files and move them in output folder.
    elif set_operation == 'intersection':
        compare_ids_ptools = set(ptools_run_ids).intersection(set(lower_case_compare_ids))

    new_compare_ids = [lower_compare_ids[compare_id] for compare_id in compare_ids_ptools]

    return new_compare_ids


def check_input_and_existing_pgdb(run_ids, input_folder, output_folder, verbose=None):
    """ Check input structure and data in output folder and ptools-local.

    Args:
        run_ids (list): species IDs (folder and GBK/GFF file name)
        input_folder (str): pathname to the input folder
        output_folder (str): pathname to the output folder
        verbose (bool): boolean verose or not
    Returns:
        list: input IDs for PathoLogic and BioPAX/dat creation
        list: input IDs for BioPAX/dat creation
    """
    # Check if there are files/folders inside the input folder.
    # And remove hidden folder/file (beginning with '.').
    species_folders = [species_folder for species_folder in os.listdir(input_folder) if not species_folder.startswith('.')]
    if len(species_folders) == 0:
        logger.critical("No folder containing genbank/gff file. In {0} you must have sub-folders containing Genbank/GFF file.".format(input_folder))
        return None, None

    # Check if there is a Genbank or a GFF file inside each subfolder.
    input_extensions = ['.gbk', '.gff']
    species_folders = list(set([species_folder for species_folder in species_folders
                                    for species_file in os.listdir(input_folder+'/'+species_folder)
                                        if any(input_extension in species_file for input_extension in input_extensions)]))
    missing_input_files = list(set(run_ids) - set(species_folders))
    if len(species_folders) == 0:
        logger.critical('Missing Genbank/GFF file for: {0} \nCheck if you have a Genbank file and if it ends with .gbk or .gff'.format(' '.join(missing_input_files)))
        return None, None

    # Check the structure of the input folder.
    invalid_characters = ['.', '/']
    for species_folder in species_folders:
        if os.path.isfile(input_folder+'/'+species_folder):
            logger.critical('Error: file inside the input_folder ({0}) instead of a subfolder. Check that you have a structure file of input_folder/species_1/species1.gbk and not input_folder/species_1.gbk.'.format(input_folder+'/'+species_folder))
            return None, None
        elif os.path.isdir(input_folder+'/'+species_folder):
            if any(char in invalid_characters for char in species_folder):
                logger.critical('Error: . or / in genbank/gff name {0} \nGenbank name is used as an ID in Pathway Tools and Pathway Tools does not create PGDB with . in ID.'.format(species_folder))
                return None, None

    # Take run_ids and remove folder with error (with the intersection with species_folders) and if there is already present output.
    clean_run_ids = set(run_ids).intersection(set(species_folders))
    if output_folder:
        already_present_outputs = [output_pgdb for output_pgdb in os.listdir(output_folder)]
        new_run_ids = clean_run_ids - set(already_present_outputs)
        new_run_ids = list(new_run_ids)
        for pgdb in already_present_outputs:
            if pgdb in clean_run_ids:
                logger.warning("! PGDB {0} already in output folder {1}, no inference will be launched on this species.".format(pgdb, output_folder))

        if len(new_run_ids) == 0:
            logger.info("All PGDBs are already present in the output folder. Remove them if you want a new inference.")
            return None, None

    else:
        new_run_ids = []
        for species_folder in species_folders:
            new_run_ids.append(species_folder)

    # Check for PGDB in ptools-local to see if PGDB are already present but they haven't been exported.
    already_present_pgdbs = [pgdb_species_folder[:-3] for pgdb_species_folder in utils.list_pgdb()]
    if already_present_pgdbs != []:
        run_patho_dat_ids = compare_input_ids_to_ptools_ids(new_run_ids, already_present_pgdbs, 'difference')
        run_dat_ids = compare_input_ids_to_ptools_ids(new_run_ids, already_present_pgdbs, 'intersection')
        for run_dat_id in run_dat_ids:
            logger.info("! PGDB {0} already in ptools-local, no PathoLogic inference will be launched on this species.".format(run_dat_id))
        return run_patho_dat_ids, run_dat_ids

    return new_run_ids, None


def create_dat_creation_script(pgdb_id, lisp_pathname):
    """ Create a lisp script allowing dat extraction.

    Args:
        pgdb_id (str): ID of a PGDB
        lisp_pathname (str): pathname to the output list script
    Returns:
        bool: True if lisp_pathname has been created
    """
    with open(lisp_pathname, 'w') as lisp_file:
        lisp_file.write("(in-package :ecocyc)")
        lisp_file.write('\n')
        lisp_file.write("(select-organism :org-id '" + pgdb_id + ")")
        lisp_file.write('\n')
        lisp_file.write('(let ((*progress-noter-enabled?* NIL))')
        lisp_file.write('\n')
        lisp_file.write("        (create-flat-files-for-current-kb))")

    return os.path.isfile(lisp_pathname)


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

    Args:
        run_folder (str): ID of a species of the input folder
    Returns:
        list: boolean list, True if all files have been created
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
                raise KeyError('No organism in the Genbank {0} In the SOURCE you must have: ORGANISM  Species name'.format(pgdb_id))

            # Take the source feature of the first record.
            # This feature contains the taxon ID in the db_xref qualifier.
            src_features = [feature for feature in first_seq_record.features if feature.type == "source"]
            for src_feature in src_features:
                try:
                    src_dbxref_qualifiers = src_feature.qualifiers['db_xref']
                    for src_dbxref_qualifier in src_dbxref_qualifiers:
                        if 'taxon:' in src_dbxref_qualifier:
                            taxon_id = src_dbxref_qualifier.replace('taxon:', '')
                except KeyError:
                    raise KeyError('No taxon ID in the Genbank {0} In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.'.format(gbk_pathname))

    elif os.path.isfile(gff_pathname):
        input_name = gff_name
        # Check if there is a fasta file.
        try:
            with open(run_folder + input_name.replace('.gff', '.fasta'), 'r'):
                gff_fasta = input_name.replace('.gff', '.fasta')
        except FileNotFoundError:
            raise FileNotFoundError('No fasta file with the GFF of {0}'.format(pgdb_id))

        # Instead of parsing and creating a database from the GFF, parse the file and extract the first region feature.
        region_feature = [feature for feature in DataIterator(gff_pathname) if feature.featuretype == 'region'][0]
        try:
            region_feature.attributes['Dbxref']
        except KeyError:
            raise KeyError('No Dbxref in GFF file of {0} GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))

        for dbxref in region_feature.attributes['Dbxref']:
            if 'taxon' in dbxref:
                taxon_id = dbxref.split('taxon:')[1]
        if not taxon_id:
            raise Exception('Missing "taxon:" in GFF file of {0} GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))


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
        if os.path.isfile(gff_pathname):
            genetic_writer.writerow(['SEQ-FILE', gff_fasta])
        genetic_writer.writerow(['//'])

    # Create the lisp script.
    check_lisp_file = create_dat_creation_script(pgdb_id, lisp_pathname)

    return all([os.path.isfile(organism_dat), os.path.isfile(genetic_dat), check_lisp_file])


def pwt_input_files(multiprocess_input):
    """
    Check if files needed by Pathway Tools are available, if not create them.
    Check if there is a pathologic.log from a previous run. If yes, delete it.

    Args:
        multiprocess_input (dict): multiprocess dictionary input
    """
    run_folder = multiprocess_input['species_input_folder_path']
    verbose = multiprocess_input['verbose']

    required_files = set(['organism-params.dat', 'genetic-elements.dat', 'dat_creation.lisp'])
    files_in = set(next(os.walk(run_folder))[2])
    if verbose:
        species_folder = run_folder.split('/')[-2]

    if "pathologic.log" in files_in:
        os.remove(run_folder + "pathologic.log")

    missing_string = ""
    if required_files.issubset(files_in):
        if verbose:
            missing_string = "no missing files"
    else:
        if verbose:
            missing_string = "missing {0}".format("; ".join(required_files.difference(files_in))) + '. Inputs file created for {0}'.format(run_folder.split('/')[-2])
        check_datas_lisp = create_dats_and_lisp(run_folder)
        if not check_datas_lisp:
            raise Exception('Error with the creation of input files of {0}'.format(run_folder))

    if verbose:
        logger.info("Checking inputs for {0}: {1}. ".format(species_folder, missing_string))


def create_mpwt_input(run_ids, input_folder, pgdbs_folder_path, verbose=None, patho_hole_filler=None, dat_extraction=None, output_folder=None, size_reduction=None, only_dat_creation=None):
    """
    Create input list for all multiprocess function, containing one lsit for each input subfolder.
    All arguments are also stored.

    Args:
        run_ids (list): input species IDs
        input_folder (str): pathname to input folder
        pgdbs_folder_path (str): pathname to species PGDB in ptools-local
        verbose (bool): verbose argument
        patho_hole_filler (bool): PathoLogic Hole Filler argument
        dat_extraction (bool): BioPAX/attribute-values file extraction argument
        output_folder (str): pathname to output folder
        size_reduction (bool): ptools-local PGDB deletion after processing argument
        only_dat_creation (bool): only create BioPAX/attribute values argument
    Returns:
        dictionary: contain all these data for multiprocessing
    """
    multiprocess_inputs = []
    for run_id in run_ids:
        multiprocess_input = {}
        input_folder_path = input_folder + "/" + run_id + "/"
        species_pgdb_folder = pgdbs_folder_path + run_id.lower() + 'cyc/'
        pgdb_id_folders = (run_id, species_pgdb_folder)
        if only_dat_creation:
            multiprocess_input['pgdb_folders'] = retrieve_complete_id(pgdb_id_folders)
        else:
            multiprocess_input['pgdb_folders'] = pgdb_id_folders
        multiprocess_input['species_input_folder_path'] = input_folder_path
        multiprocess_input['verbose'] = verbose
        multiprocess_input['patho_hole_filler'] = patho_hole_filler
        multiprocess_input['dat_extraction'] = dat_extraction
        multiprocess_input['output_folder'] = output_folder
        multiprocess_input['size_reduction'] = size_reduction
        multiprocess_inputs.append(multiprocess_input)

    return multiprocess_inputs


def create_only_dat_lisp(pgdbs_folder_path, tmp_folder):
    """
    Create a lisp script file for each PGDB in the ptools-local folder.
    Return a generator with the PGDB IDs.

    Args:
        pgdbs_folder_path (str): pathname to species PGDB in ptools-local
        tmp_folder (str): termporary folder where lisp script will be stored
    Returns:
        generator: generator with the PGDB IDs.
    """
    for species_pgdb in os.listdir(pgdbs_folder_path):
        if os.path.isdir(pgdbs_folder_path + species_pgdb):
            pgdb_id = species_pgdb[:-3]
            pgdb_pathname = tmp_folder + pgdb_id + '/'
            os.mkdir(tmp_folder + pgdb_id)
            lisp_pathname = pgdb_pathname + "dat_creation.lisp"
            check_lisp_file = create_dat_creation_script(pgdb_id, lisp_pathname)
            if not check_lisp_file:
                raise Exception('Error with the creation of the lisp script for {0}'.format(species_pgdb))

            yield pgdb_id


def permission_change(folder_pathname):
    """
    Give permission to output files inside a folder.
    Used for log files and PGDB/dat files.

    Args:
        folder_pathname (str): pathname to the folder which permissions will be changed
    """
    os.chmod(folder_pathname, 0o777)
    for root, subfolders, subfiles in os.walk(folder_pathname):
        for subfolder in subfolders:
            os.chmod(os.path.join(root, subfolder), 0o777)
        for subfile in subfiles:
            os.chmod(os.path.join(root, subfile), 0o777)


def check_pwt(multiprocess_inputs, patho_log_folder):
    """
    Check PathoLogic's log.
    Create two log files (log_error.txt which contains Pathway Tools log and resume_inference which contains summary of network).

    Args:
        multiprocess_inputs (list): list of dictionary contaning multiprocess input data
        patho_log_folder (str): pathname to the PathoLogic log folder.
    """
    verbose = multiprocess_inputs[0]['verbose']

    if patho_log_folder:
        if not os.path.exists(patho_log_folder):
            logger.info('No log directory, it will be created.')
            os.mkdir(patho_log_folder)

        patho_error_pathname = patho_log_folder + '/log_error.txt'
        patho_resume_pathname = patho_log_folder + '/resume_inference.tsv'

        patho_error_file = open(patho_error_pathname, 'w')
        patho_resume_file = open(patho_resume_pathname, 'w')
        patho_resume_writer = csv.writer(patho_resume_file, delimiter='\t', lineterminator='\n')
        patho_resume_writer.writerow(['species', 'gene_number', 'protein_number', 'pathway_number', 'reaction_number', 'compound_number'])

    failed_inferences = []
    passed_inferences = []

    for multiprocess_input in multiprocess_inputs:
        species_input_folder_path = multiprocess_input['species_input_folder_path']
        species = species_input_folder_path.split('/')[-2]
        patho_log = species_input_folder_path + '/pathologic.log'

        if patho_log_folder:
            patho_error_file.write('------------ Species: ')
            patho_error_file.write(species)
            patho_error_file.write('\n')

        fatal_error_index = None

        if os.path.exists(patho_log):
            with open(patho_log, 'r') as input_file:
                for index, line in enumerate(input_file):
                    if 'fatal error' in line or 'Error' in line:
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
                if species not in passed_inferences and species not in failed_inferences:
                    failed_inferences.append(species)
                    if patho_log_folder:
                        patho_error_file.write('No build in PathoLogic inference.')
                        patho_resume_writer.writerow([species, 'ERROR', '', '', '', ''])
        else:
            if patho_log_folder:
                patho_error_file.write('No pathologic log, an error occured before PathoLogic run.\n')
                patho_resume_writer.writerow([species, 'ERROR', '', '', '', ''])
            logger.info('No pathologic log for {0}, an error occured before PathoLogic run.'.format(species))

        if patho_log_folder:
            patho_error_file.write('------------\n\n')

    number_passed_inference = len(passed_inferences)
    number_failed_inference = len(failed_inferences)

    string_passed_build = 'build has' if number_passed_inference == 1 else 'builds have'
    string_failed_build = 'build has' if number_failed_inference == 1 else 'builds have'

    if number_passed_inference > 0:
        if verbose:
            logger.info('\n{0} {1} passed!\n'.format(str(number_passed_inference), string_passed_build))
    if number_failed_inference > 0:
        if verbose:
            logger.critical('WARNING: {0} {1} failed! See the log for more information.\n'.format(str(number_failed_inference), string_failed_build))

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


def retrieve_complete_id(pgdb_id_folder):
    """
    Retrieve the ID of the PGDB from the genetic-elements.dat file.

    Args:
        pgdb_id_folder (tuple): second tuple argument is the pathname to the PGDB
    Returns:
        tuple: (new PGDB ID (according to input file), pathname to PGDB folder)
    """
    with open(pgdb_id_folder[1] + '/1.0/input/genetic-elements.dat') as organism_file:
        for line in organism_file:
            if 'ANNOT-FILE' in line and ';;' not in line:
                pgdb_id_complete = line.split('\t')[1].replace('.gff','').replace('.gbk','').strip()

    return (pgdb_id_complete, pgdb_id_folder[1])


def pwt_error(species_input_folder_path, subprocess_returncode, subprocess_stdout, subprocess_stderr, cmd):
    """
    Print error messages when there is a subprocess error during PathoLogic run.

    Args:
        species_input_folder_path (str): pathname to the species input folder
        subprocess_returncode (int): code returns by subprocess
        subprocess_stdout (str): stdout returns by subprocess
        cmd (str): command used which generates the error
    """
    logger.critical('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')
    species_name = species_input_folder_path.split('/')[-2]
    logger.critical('Error for {0} with PathoLogic subprocess, return code: {1}'.format(species_name, str(subprocess_returncode)))
    if subprocess_stderr:
        logger.critical('An error occurred :' + subprocess_stderr.decode('utf-8'))

    logger.critical('=== Pathway Tools log ===')
    for line in subprocess_stdout:
        if line != '':
            logger.critical('\t' + line)

    # Look for error in pathologic.log.
    if '-patho' in cmd:
        fatal_error_index = None
        with open(species_input_folder_path + '/pathologic.log', 'r') as pathologic_log:
            for index, line in enumerate(pathologic_log):
                if line != '':
                    if 'fatal error' in line or 'Error' in line and not fatal_error_index:
                        fatal_error_index = index
                        logger.critical('=== Error in Pathologic.log ===')
                        logger.critical('\t' + 'Error from the pathologic.log file: {0}'.format(species_input_folder_path + '/pathologic.log'))
                        logger.critical('\t' + line)
                    if fatal_error_index:
                        if index > fatal_error_index:
                            logger.critical('\t' + line)

    logger.critical('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')


def run_pwt(multiprocess_input):
    """
    Create PGDB using files created during 'create_dats_and_lisp' ('organism-params.dat' and 'genetic-elements.dat').
    With verbose run check_output to retrieve the output of subprocess (and show when Pathway Tools has been killed).
    Otherwise send the output to the null device.
    Command used:
    pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho

    Args:
        multiprocess_input (dictionary): contains multiprocess input (mpwt argument: input folder, output folder, ...)
    """
    species_input_folder_path = multiprocess_input['species_input_folder_path']
    verbose = multiprocess_input['verbose']
    patho_hole_filler = multiprocess_input['patho_hole_filler']

    cmd_options = ['-no-web-cel-overview', '-no-cel-overview', '-no-patch-download', '-disable-metadata-saving', '-nologfile']

    cmd_pwt = ['pathway-tools', *cmd_options, '-patho', species_input_folder_path]

    if patho_hole_filler:
        cmd_pwt.append('-hole-filler')

    if verbose:
        logger.info(' '.join(cmd_pwt))

    error_status = None
    errors = ['Restart actions (select using :continue):', 'Error']
    patho_lines = []

    try:
        patho_subprocess = subprocess.Popen(cmd_pwt, stdout=subprocess.PIPE, universal_newlines="")
        # Check internal error of Pathway Tools (Error with Genbank).
        for patho_line in iter(patho_subprocess.stdout.readline, ""):
            patho_line = patho_line.decode('utf-8')
            if any(error in patho_line for error in errors):
                logger.info('Error possibly with the genbank file.')
                error_status = True
                patho_subprocess.kill()

            patho_lines.append(patho_line)
            patho_subprocess.poll()
            return_code = patho_subprocess.returncode
            # Check if Pathway Tools has been killed with returncode.
            # Also check if Pathway Tools has finished PathoLogic inference (returncode 0).
            if return_code or return_code == 0:
                if return_code == 0:
                    patho_subprocess.stdout.close()
                    return error_status
                elif return_code != 0:
                    raise subprocess.CalledProcessError(return_code, cmd_pwt)

    except subprocess.CalledProcessError as subprocess_error:
        # Check error with subprocess (when process is killed).
        pwt_error(species_input_folder_path, subprocess_error.returncode, patho_lines, patho_subprocess.stderr, cmd_pwt)
        error_status = True
    patho_subprocess.stdout.close()

    return error_status


def run_pwt_dat(multiprocess_input):
    """
    Create dat file using a lisp script created during 'create_dats_and_lisp'.
    Kill the subprocess when the command reach the Navigator Window opening proposition.
    If this proposition is not closed, the script can't continue.
    Command used:
    pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load

    Args:
        multiprocess_input (dictionary): contains multiprocess input (mpwt argument: input folder, output folder, ...)
    """
    species_input_folder_path = multiprocess_input['species_input_folder_path']
    verbose = multiprocess_input['verbose']

    lisp_path = species_input_folder_path + '/dat_creation.lisp'
    cmd_options = ['-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_dat = ['pathway-tools', *cmd_options, '-load', lisp_path]

    if verbose:
        logger.info(' '.join(cmd_dat))

    error_status = None
    dat_creation_ends = ['Opening Navigator window.', 'No protein-coding genes with sequence data found.  Cannot continue.']
    load_errors = ['Error']
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
        pwt_error(species_input_folder_path, subprocess_error.returncode, load_lines, load_subprocess.stderr, cmd_dat)
        error_status = True
    load_subprocess.stdout.close()

    return error_status


def check_dat(multiprocess_input):
    """
    Check dats creation.

    Args:
        multiprocess_input (dictionary): contains multiprocess input (mpwt argument: input folder, output folder, ...)
    """
    pgdb_folder = multiprocess_input['pgdb_folders']
    verbose = multiprocess_input['verbose']
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
    if verbose:
        expected_dat_number = str(len(dat_files))
        found_dat_number = str(len(dat_checks))
        logger.info('{0}: {1} out of {2} dat files create.'.format(pgdb_folder_dbname, found_dat_number, expected_dat_number))


def run_move_pgdb(move_data):
    """
    Move the result files inside the shared folder containing the input data.
    pgdb_folder_dbname: ID of the species.
    pgdb_folder_path: path to the PGDB of the species (in ptools-local).

    Args:
        move_data (dictionary): contains multiprocess input (PGDB ID, pttols-local PGDB pathname, ...)
    """
    pgdb_folder_dbname = move_data['pgdb_folders'][0]
    pgdb_folder_path = move_data['pgdb_folders'][1]
    dat_extraction = move_data['dat_extraction']
    output_folder = move_data['output_folder']
    size_reduction = move_data['size_reduction']

    output_species = output_folder + '/' + pgdb_folder_dbname +'/'

    if dat_extraction:
        pgdb_folder_path = pgdb_folder_path + '/1.0/data'

    if size_reduction:
        if not os.path.exists(output_species):
            os.mkdir(output_species)
        for pgdb_file in os.listdir(pgdb_folder_path):
            file_to_move_pathname = pgdb_folder_path + '/' + pgdb_file
            output_file_pathname = output_species + pgdb_file
            if dat_extraction:
                if '.dat' in pgdb_file:
                    shutil.move(file_to_move_pathname, output_file_pathname)
            else:
                shutil.move(file_to_move_pathname, output_file_pathname)
        shutil.rmtree(pgdb_folder_path)
    else:
        shutil.copytree(pgdb_folder_path, output_species)
        if dat_extraction:
            for pgdb_file in os.listdir(output_species):
                if '.dat' not in pgdb_file:
                    os.remove(output_species+'/'+pgdb_file)


def multiprocess_pwt(input_folder=None, output_folder=None, patho_inference=None, patho_hole_filler=None, dat_creation=None, dat_extraction=None, size_reduction=None, number_cpu=None, patho_log=None, verbose=None):
    """
    Function managing all the workflow (from the creatin of the input files to the results).
    Use it when you import mpwt in a script.

    Args:
        input_folder (str): pathname to input folder
        output_folder (str): pathname to output folder
        patho_inference (bool): pathologic boolean (True/False)
        patho_hole_filler (bool): pathologic hole filler boolean (True/False)
        dat_creation (bool): BioPAX/attributes-values files creation boolean (True/False)
        dat_extraction (bool): BioPAX/attributes-values files extraction boolean (True/False)
        size_reduction (bool): Delete ptools-local data at the end boolean (True/False)
        number_cpu (int): number of CPU used (default=1)
        patho_log (str): pathname to mpwt log folder
        verbose (bool): verbose argument
    """
    # Check if Pathway Tools is in the path.
    # Find PGDB folder path.
    ptools_local_path = utils.find_ptools_path()
    pgdbs_folder_path = ptools_local_path + '/pgdbs/user/'

    # Check if patho_hole_filler or patho_log are launched with patho_inference.
    if (patho_hole_filler and not patho_inference) or (patho_log and not patho_inference):
        sys.exit('To use either --hf/patho_hole_filler or --log/patho_log, you need to add the --patho/patho_inference argument.')

    # Use the number of cpu given by the user or all the cpu available.
    if number_cpu:
        number_cpu_to_use = int(number_cpu)
    else:
        number_cpu_to_use = 1
    mpwt_pool = Pool(processes=number_cpu_to_use)

    # Check input folder and create input files for PathoLogic.
    if input_folder:
        run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
        if output_folder:
            if not os.path.exists(output_folder):
                if verbose:
                    logger.info('No output directory, it will be created.')
                os.mkdir(output_folder)
        run_patho_dat_ids, run_dat_ids = check_input_and_existing_pgdb(run_ids, input_folder, output_folder, verbose)

        # Check if some inputs need to be process by PathoLogic.
        if run_patho_dat_ids:
            # Create the list containing all the data used by the multiprocessing call.
            multiprocess_inputs = create_mpwt_input(run_patho_dat_ids, input_folder, pgdbs_folder_path, verbose, patho_hole_filler, dat_extraction, output_folder, size_reduction)

            if verbose:
                logger.info('~~~~~~~~~~Creation of input data from Genbank/GFF~~~~~~~~~~')
            mpwt_pool.map(pwt_input_files, multiprocess_inputs)

            # Launch PathoLogic.
            if patho_inference:
                if verbose:
                    logger.info('~~~~~~~~~~Inference on the data~~~~~~~~~~')
                error_status = mpwt_pool.map(run_pwt, multiprocess_inputs)
                if verbose:
                    logger.info('~~~~~~~~~~Check inference~~~~~~~~~~')
                check_pwt(multiprocess_inputs, patho_log)
                if any(error_status):
                    sys.exit('Error during inference. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')
        else:
            multiprocess_inputs = []

    # Create path for lisp if there is no folder given.
    # Create the input for the creaetion of BioPAX/attribute-values files.
    if (dat_creation and not input_folder) or (output_folder and not input_folder):
        only_dat_creation = True
        # Create a temporary folder in ptools-local where list script will be stored.
        tmp_folder = ptools_local_path + '/tmp/'
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)

        # Create a lisp script file for each PGDB in the ptools-local folder.
        dat_run_ids = create_only_dat_lisp(pgdbs_folder_path, tmp_folder)

        multiprocess_inputs = create_mpwt_input(dat_run_ids, tmp_folder, pgdbs_folder_path, verbose, patho_hole_filler, dat_extraction, output_folder, size_reduction, only_dat_creation)

    # Add species that have data in PGDB but are not present in output folder.
    if input_folder:
        if run_dat_ids:
            for run_dat_id in run_dat_ids:
                create_dat_creation_script(run_dat_id, input_folder + "/" + run_dat_id + "/" + "dat_creation.lisp")
            multiprocess_dat_inputs = create_mpwt_input(run_dat_ids, input_folder, pgdbs_folder_path, verbose, patho_hole_filler, dat_extraction, output_folder, size_reduction)
            multiprocess_inputs.extend(multiprocess_dat_inputs)

    # Create BioPAX/attributes-values dat files.
    if (input_folder and dat_creation) or dat_creation:
        if verbose:
            logger.info('~~~~~~~~~~Creation of the .dat files~~~~~~~~~~')
        mpwt_pool.map(run_pwt_dat, multiprocess_inputs)
        if verbose:
            logger.info('~~~~~~~~~~Check .dat~~~~~~~~~~')
        for multiprocess_input in multiprocess_inputs:
            check_dat(multiprocess_input)

    if (dat_creation and not input_folder) or (output_folder and not input_folder):
        ptools_local_path = utils.find_ptools_path()
        shutil.rmtree(ptools_local_path + '/tmp')

    if verbose:
        logger.info('~~~~~~~~~~End of the Pathway Tools Inference~~~~~~~~~~')

    # Move PGDBs files.
    if output_folder:
        if verbose:
            logger.info('~~~~~~~~~~Moving result files~~~~~~~~~~')
        mpwt_pool.map(run_move_pgdb, multiprocess_inputs)
        # Give access to the file for user outside the container.
        permission_change(output_folder)

    mpwt_pool.close()
    mpwt_pool.join()

    if verbose:
        logger.info('~~~~~~~~~~mpwt has finished! Thank you for using it.')


def run_mpwt():
    """
    Function used with a mpwt call in the terminal.
    """
    args = docopt.docopt(__doc__)

    argument_number = len(sys.argv[1:])

    input_folder = args['-f']
    output_folder = args['-o']
    patho_inference = args['--patho']
    patho_hole_filler = args['--hf']
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
            logger.info('No PGDB inside ptools-local.')
        else:
            logger.info(str(len(pgdbs)) + ' PGDB inside ptools-local:\n' + '\t'.join(pgdbs))
        return

    # Delete PGDB if use of --delete argument.
    # Use a set to remove redudant PGDB.
    if pgdb_to_deletes:
        utils.remove_pgbds(list(set(pgdb_to_deletes.split(','))), number_cpu)
        return

    if args['--clean']:
        if verbose:
            logger.info('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')
        utils.cleaning(number_cpu, verbose)
        if input_folder:
            utils.cleaning_input(input_folder, verbose)
        if argument_number == 1 or (argument_number == 2 and verbose) or (argument_number == 3 and verbose and number_cpu):
            sys.exit()

    multiprocess_pwt(input_folder,
                    output_folder,
                    patho_inference,
                    patho_hole_filler,
                    dat_creation,
                    move_dat,
                    size_reduction,
                    number_cpu,
                    patho_log,
                    verbose)


if __name__ == '__main__':
    run_mpwt()
