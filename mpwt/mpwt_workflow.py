"""
Workflow of mpwt:
    -create input files (pathologic_input)
    -launch Pathway Tools in multiprocess (pwt_wrapper)
    -check the results (results_check)
"""

import logging
import os
import shutil
import sys
import time

from mpwt import utils
from mpwt.pwt_wrapper import run_pwt, run_pwt_dat, run_move_pgdb
from mpwt.results_check import check_dat, check_pwt, permission_change
from mpwt.pathologic_input import check_input_and_existing_pgdb, create_mpwt_input, pwt_input_files, create_only_dat_lisp, create_dat_creation_script
from multiprocessing import Pool

logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
logger = logging.getLogger('mpwt')
logger.setLevel(logging.CRITICAL)


def multiprocess_pwt(input_folder=None, output_folder=None, patho_inference=None,
                     patho_hole_filler=None, dat_creation=None, dat_extraction=None,
                     size_reduction=None, number_cpu=None, patho_log=None,
                     ignore_error=None, taxon_file=None, verbose=None):
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
    if verbose:
        logger.setLevel(logging.DEBUG)

    start_time = time.time()
    times = []
    steps = []
    times.append(start_time)
    steps.append('start')

    # Check if Pathway Tools is in the path.
    # Find PGDB folder path.
    ptools_local_path = utils.find_ptools_path()
    pgdbs_folder_path = ptools_local_path + '/pgdbs/user/'

    # Check if patho_hole_filler or patho_log are launched with patho_inference.
    if (patho_hole_filler and not patho_inference) or (patho_log and not patho_inference):
        sys.exit('To use either --hf/patho_hole_filler or --log/patho_log, you need to add the --patho/patho_inference argument.')

    # Check if size_reduction is used with output_folder.
    if size_reduction and not output_folder:
        sys.exit('To use -r/size_reduction, you need to give an output folder (-o/output_folder).')

    # Check if ignore_error is used with patho_inference.
    if ignore_error and not patho_inference:
        sys.exit('To use --ignore-error/ignore_error, you need to use the --patho/patho_inference argument.')

    # Check if taxon_file is used with patho_inference.
    if taxon_file and not patho_inference:
        sys.exit('To use --taxon-file/taxon_file, you need to use the --patho/patho_inference argument.')

    # Use the number of cpu given by the user or 1 CPU.
    if number_cpu:
        try:
            number_cpu_to_use = int(number_cpu)
        except ValueError:
            raise ValueError('The number of CPU must be an integer.')
    else:
        number_cpu_to_use = 1
    mpwt_pool = Pool(processes=number_cpu_to_use)

    # Check input folder and create input files for PathoLogic.
    if input_folder:
        run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
        if output_folder:
            if not os.path.exists(output_folder):
                logger.info('No output directory, it will be created.')
                os.mkdir(output_folder)
        run_patho_dat_ids, run_dat_ids = check_input_and_existing_pgdb(run_ids, input_folder, output_folder)

        # Launch PathoLogic inference on species with no PGDBs.
        if run_patho_dat_ids:
            # Create the list containing all the data used by the multiprocessing call.
            multiprocess_inputs = create_mpwt_input(run_ids=run_patho_dat_ids, input_folder=input_folder, pgdbs_folder_path=pgdbs_folder_path,
                                                    patho_hole_filler=patho_hole_filler, dat_extraction=dat_extraction, output_folder=output_folder,
                                                    size_reduction=size_reduction, only_dat_creation=None, taxon_file=taxon_file)

            logger.info('~~~~~~~~~~Creation of input data from Genbank/GFF/PF~~~~~~~~~~')
            mpwt_pool.map(pwt_input_files, multiprocess_inputs)

            input_time = time.time()
            times.append(input_time)
            steps.append('pwt input creation')
            logger.info('----------End of creation of input data from Genbank/GFF/PF: {0:.2f}s----------'.format(times[-1] - times[-2]))

            # Launch PathoLogic.
            if patho_inference:
                logger.info('~~~~~~~~~~Inference on the data~~~~~~~~~~')
                error_status = mpwt_pool.map(run_pwt, multiprocess_inputs)

                # Check PathoLogic build.
                logger.info('~~~~~~~~~~Check inference~~~~~~~~~~')
                passed_inferences = check_pwt(multiprocess_inputs, patho_log)
                if any(error_status):
                    if ignore_error:
                        logger.critical('Error during inference. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')
                    else:
                        sys.exit('Error during inference. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')

            patho_time = time.time()
            times.append(patho_time)
            steps.append('PathoLogic inference')
            logger.info('----------End of PathoLogic inference: {0:.2f}s----------'.format(times[-1] - times[-2]))
        else:
            multiprocess_inputs = []

    # Create path for lisp if there is no folder given.
    # Create the input for the creation of BioPAX/attribute-values files.
    if (dat_creation and not input_folder) or (output_folder and not input_folder):
        only_dat_creation = True
        # Create a temporary folder in ptools-local where list script will be stored.
        tmp_folder = ptools_local_path + '/tmp/'
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)

        # Create a lisp script file for each PGDB in the ptools-local folder.
        dat_run_ids = create_only_dat_lisp(pgdbs_folder_path, tmp_folder)

        multiprocess_inputs = create_mpwt_input(run_ids=dat_run_ids, input_folder=tmp_folder, pgdbs_folder_path=pgdbs_folder_path,
                                                patho_hole_filler=patho_hole_filler, dat_extraction=dat_extraction, output_folder=output_folder,
                                                size_reduction=size_reduction, only_dat_creation=only_dat_creation, taxon_file=taxon_file)
    # Add species that have data in PGDB but are not present in output folder.
    # Or if ignore_error has been used, select only PathoLogic build that have succeed + species in input with PGDB and not in output.
    if input_folder:
        if ignore_error:
            multiprocess_inputs = []
            tmp_run_dat_ids = list(set(passed_inferences).intersection(set(run_patho_dat_ids)))
            tmp_run_dat_ids.extend(run_dat_ids)
            run_dat_ids = tmp_run_dat_ids
        if run_dat_ids:
            for run_dat_id in run_dat_ids:
                create_dat_creation_script(run_dat_id, input_folder + "/" + run_dat_id + "/" + "dat_creation.lisp")
            multiprocess_dat_inputs = create_mpwt_input(run_ids=run_dat_ids, input_folder=input_folder, pgdbs_folder_path=pgdbs_folder_path,
                                                        patho_hole_filler=patho_hole_filler, dat_extraction=dat_extraction, output_folder=output_folder,
                                                        size_reduction=size_reduction, only_dat_creation=None, taxon_file=taxon_file)
            multiprocess_inputs.extend(multiprocess_dat_inputs)

    # Create BioPAX/attributes-values dat files.
    if (input_folder and dat_creation) or dat_creation:
        logger.info('~~~~~~~~~~Creation of the .dat files~~~~~~~~~~')
        dat_error_status = mpwt_pool.map(run_pwt_dat, multiprocess_inputs)
        logger.info('~~~~~~~~~~Check .dat~~~~~~~~~~')
        for multiprocess_input in multiprocess_inputs:
            check_dat(multiprocess_input)
        if any(dat_error_status):
            if ignore_error:
                logger.critical('Error during dat creation. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')
            else:
                sys.exit('Error during dat creation. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')

        dat_time = time.time()
        times.append(dat_time)
        steps.append('BioPAX/attribute-value dat files creation')
        logger.info('----------End of dat files creation: {0:.2f}s----------'.format(times[-1] - times[-2]))

    if (dat_creation and not input_folder) or (output_folder and not input_folder):
        ptools_local_path = utils.find_ptools_path()
        shutil.rmtree(ptools_local_path + '/tmp')

    logger.info('~~~~~~~~~~End of Pathway Tools~~~~~~~~~~')

    # Move PGDBs or attribute-values/dat files.
    if output_folder:
        logger.info('~~~~~~~~~~Moving result files~~~~~~~~~~')
        mpwt_pool.map(run_move_pgdb, multiprocess_inputs)
        # Give access to the file for user outside the container.
        permission_change(output_folder)

        move_time = time.time()
        times.append(move_time)
        steps.append('Moving results files')
        logger.info('----------End of moving fimes: {0:.2f}s----------'.format(times[-1] - times[-2]))


    mpwt_pool.close()
    mpwt_pool.join()

    end_time = time.time()
    times.append(end_time)
    steps.append('mpwt')

    # Write each step time in log file.
    if patho_log:
        if patho_log:
            if not os.path.exists(patho_log):
                logger.info('No log directory, it will be created.')
                os.mkdir(patho_log)
        patho_error_pathname = patho_log + '/log_error.txt'
        with open(patho_error_pathname, 'a') as input_file:
            input_file.write('\n\n---------Time---------\n')
            for index, step_time in enumerate(times):
                if index != 0:
                    if index + 1 == len(times):
                        step_duration = step_time - times[0]
                    else:
                        step_duration = step_time - times[index-1]
                    input_file.write('Step {0} takes: {1:.2f}s.\n'.format(steps[index] , step_duration))

        permission_change(patho_log)

    logger.info('----------mpwt has finished in {0:.2f}s! Thank you for using it.'.format(end_time - start_time))
