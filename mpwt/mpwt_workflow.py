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
    -f=DIR     Working folder containing sub-folders with Genbank/GFF/PF files.
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

import logging
import os
import shutil
import sys

from mpwt import utils
from mpwt.pwt_wrapper import run_pwt, run_pwt_dat, run_move_pgdb
from mpwt.results_check import check_dat, check_pwt, permission_change
from mpwt.pathologic_input import check_input_and_existing_pgdb, create_mpwt_input, pwt_input_files, create_only_dat_lisp, create_dat_creation_script
from multiprocessing import Pool

logging.basicConfig(format='%(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


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
                logger.info('~~~~~~~~~~Creation of input data from Genbank/GFF/PF~~~~~~~~~~')
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
