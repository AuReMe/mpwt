# Copyright (C) 2018-2021 Arnaud Belcour - Inria Dyliss
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

"""
Workflow of mpwt:
    - create PathoLogic input files: organisms.dat and genetic-elements.dat (pathologic_input)
    - mutliprocessing of PathoLogic (run_pwt)
    - create flat files (run_pwt_flat)
    - check the results (results_check)
    - move the results to output folder (run_move_pgdb)
"""

import csv
import logging
import os
import shutil
import stat
import sys
import time

from mpwt import utils
from mpwt.pwt_wrapper import run_pwt, run_pwt_flat, run_move_pgdb
from mpwt.results_check import check_dat, check_pwt, check_mpwt_pathologic_runs
from mpwt.pathologic_input import check_input_and_existing_pgdb, pwt_input_files, create_only_flat_lisp, create_flat_creation_script, read_taxon_id
from multiprocessing import Pool

logger = logging.getLogger(__name__)


def multiprocess_pwt(input_folder=None, output_folder=None, patho_inference=None,
                     patho_hole_filler=None, patho_operon_predictor=None,
                     patho_transporter_inference=None, no_download_articles=None,
                     flat_creation=None, dat_extraction=None, xml_extraction=None,
                     owl_extraction=None, col_extraction=None, size_reduction=None,
                     number_cpu=None, patho_log=None, pathway_score=None,
                     taxon_file=None, verbose=None, permission=None):
    """
    Function managing all the workflow (from the creatin of the input files to the results).
    Use it when you import mpwt in a script.

    Args:
        input_folder (str): pathname to input folder
        output_folder (str): pathname to output folder
        patho_inference (bool): PathoLogic inference (True/False)
        patho_hole_filler (bool): PathoLogic hole filler (True/False)
        patho_operon_predictor (bool): PathoLogic operon predictor (True/False)
        patho_transporter_inference (bool): PathoLogic Transport Inference Parser (True/False)
        no_download_articles (bool): turning off loading of PubMed citations (True/False)
        flat_creation (bool): BioPAX/attributes-values flat files creation (True/False)
        dat_extraction (bool): move BioPAX/attributes-values files to output folder (True/False)
        xml_extraction (bool): move metabolic-reactions.xml to output folder (True/False)
        owl_extraction (bool): move owl files to output folder (True/False)
        col_extraction (bool): move tabular files to output folder (True/False)
        size_reduction (bool): delete ptools-local data at the end (True/False)
        number_cpu (int): number of CPU used (default=1)
        patho_log (str): pathname to mpwt log folder
        pathway_score (float): score between 0 and 1 to accept or reject pathway
        taxon_file (str): pathname to the mpwt taxon ID file
        verbose (bool): verbose argument
        permission (str): Choose permission access to PGDB in ptools-local and output files, either 'all' or 'group' (by default it is user).
    """
    if verbose:
        logger.setLevel(logging.DEBUG)

    # Check if ptools-local is accessible.
    error = utils.check_ptools_local_pwt()
    if error:
        sys.exit(1)

    # Check if patho_hole_filler or patho_log are launched with patho_inference.
    if (patho_hole_filler and not patho_inference) or (patho_log and not patho_inference):
        sys.exit('To use either --hf/patho_hole_filler or --log/patho_log, you need to add the --patho/patho_inference argument.')

    # Check if size_reduction is used with output_folder.
    if size_reduction and not output_folder:
        sys.exit('To use -r/size_reduction, you need to give an output folder (-o/output_folder).')

    # Check if taxon_file is used with patho_inference.
    if (taxon_file and not patho_inference) and (taxon_file and not input_folder):
        sys.exit('To use --taxon-file/taxon_file, you need to use the --patho/patho_inference argument. Or you can use it with the -f argument to create the taxon file from data.')

    # Check if patho_operon_predictor is used with patho_inference.
    if patho_operon_predictor and not patho_inference:
        sys.exit('To use --op/patho_operon_predictor, you need to use the --patho/patho_inference argument.')

    # Check if no_download_articles is used with patho_inference.
    if no_download_articles and not patho_inference:
        sys.exit('To use --nc/no_download_articles, you need to use the --patho/patho_inference argument.')

    # Check if patho_transporter_inference is used with patho_inference.
    if patho_transporter_inference and not patho_inference:
        sys.exit('To use --tp/patho_transporter_inference, you need to use the --patho/patho_inference argument.')

    # Check if pathway_score is used with patho_inference.
    if pathway_score and not patho_inference:
        sys.exit('To use -p/pathway_score, you need to use the --patho/patho_inference argument.')

    # Check if pathway_score is a float between 0 and 1.
    if pathway_score:
        try:
            pathway_score = float(pathway_score)
        except ValueError:
            sys.exit("{0} is not a float. Pathway score prediction must be a flaot between 0.0 and 1.0.".format(pathway_score))

        if pathway_score < 0.0 or pathway_score > 1.0:
            sys.exit("{0} is not a float between 0.0 and 1.0.".format(pathway_score))

    # Check if dat_extraction is used with output_folder.
    if dat_extraction and not output_folder:
        sys.exit('To use --md/dat_extraction, you need to use the -o/output_folder argument.')

    # Check if xml_extraction is used with output_folder.
    if xml_extraction and not output_folder:
        sys.exit('To use --mx/xml_extraction, you need to use the -o/output_folder argument.')

    # Check if owl_extraction is used with output_folder.
    if owl_extraction and not output_folder:
        sys.exit('To use --mo/owl_extraction, you need to use the -o/output_folder argument.')

    # Check if col_extraction is used with output_folder.
    if col_extraction and not output_folder:
        sys.exit('To use --mc/col_extraction, you need to use the -o/output_folder argument.')

    # Check permission.
    if permission and permission not in ['all', 'group']:
        sys.exit('--permission/permission must be either "group" or "all".')

    # Use the number of cpu given by the user or 1 CPU.
    if number_cpu:
        try:
            number_cpu_to_use = int(number_cpu)
        except ValueError:
            raise ValueError('The number of CPU must be an integer.')
    else:
        number_cpu_to_use = 1

    independent_mpwt(input_folder, output_folder, patho_inference,
                        patho_hole_filler, patho_operon_predictor,
                        patho_transporter_inference, no_download_articles,
                        flat_creation, dat_extraction, xml_extraction,
                        owl_extraction, col_extraction, size_reduction,
                        number_cpu_to_use, patho_log, pathway_score,
                        taxon_file, permission)


def close_mpwt(mpwt_pool, no_download_articles, pathway_score=None, old_pathway_score=None):
    """End multiprocessing Pool and restore ptools-init.dat

    mpwt_pool (multiprocessing Pool): mpwt multiprocessing Pool
    no_download_articles (bool): turning off loading of PubMed citations (True/False)
    pathway_score (float): score between 0 and 1 to accept or reject pathway
    """
    mpwt_pool.close()
    mpwt_pool.join()

    # Turn on loading of pubmed entries.
    if no_download_articles:
        utils.pubmed_citations(activate_citations=True)

    # Remodify the pathway score to its original value.
    if pathway_score:
        utils.modify_pathway_score(old_pathway_score, comment_line=True)


def classic_mpwt(input_folder, output_folder=None, patho_inference=None,
                     patho_hole_filler=None, patho_operon_predictor=None,
                     patho_transporter_inference=None, no_download_articles=None,
                     flat_creation=None, dat_extraction=None, xml_extraction=None,
                     owl_extraction=None, col_extraction=None, size_reduction=None,
                     number_cpu_to_use=None, patho_log=None, ignore_error=None,
                     pathway_score=None, taxon_file=None, permission=None):
    """
    Run mpwt in a "dependent" way (the classic run of mpwt).
    This means that all the PathoLogic must succeed to move to the flat files creation steps and then to moving the output files.

    Args:
        input_folder (str): pathname to input folder
        output_folder (str): pathname to output folder
        patho_inference (bool): PathoLogic inference (True/False)
        patho_hole_filler (bool): PathoLogic hole filler (True/False)
        patho_operon_predictor (bool): PathoLogic operon predictor (True/False)
        patho_transporter_inference (bool): PathoLogic Transport Inference Parser (True/False)
        no_download_articles (bool): turning off loading of PubMed citations (True/False)
        flat_creation (bool): BioPAX/attributes-values flat files creation (True/False)
        dat_extraction (bool): move BioPAX/attributes-values files to output folder (True/False)
        xml_extraction (bool): move metabolic-reactions.xml to output folder (True/False)
        owl_extraction (bool): move owl files to output folder (True/False)
        col_extraction (bool): move tabular files to output folder (True/False)
        size_reduction (bool): delete ptools-local data at the end (True/False)
        number_cpu (int): number of CPU used (default=1)
        patho_log (str): pathname to mpwt log folder
        ignore_error (bool): Ignore error during PathoLogic inference (True/False)
        pathway_score (float): score between 0 and 1 to accept or reject pathway
        taxon_file (str): pathname to the mpwt taxon ID file
        permission (str): Choose permission access to PGDB in ptools-local and output files, either 'all' or 'group' (by default it is user).
    """
    start_time = time.time()
    times = []
    steps = []
    times.append(start_time)
    steps.append('start')

    # Check if Pathway Tools is in the path.
    # Find PGDB folder path.
    ptools_local_path = utils.find_ptools_path()
    pgdbs_folder_path = os.path.join(*[ptools_local_path, 'pgdbs', 'user'])

    mpwt_pool = Pool(processes=number_cpu_to_use)

    if input_folder:
        if not os.path.exists(input_folder):
            logger.critical('mpwt can not run: ' + input_folder + ' does not exist.')
            return
        if not os.path.isdir(input_folder):
            logger.critical('mpwt can not run: ' + input_folder + ' is not a directory.')
            return

    # Create taxon file in the input folder.
    if taxon_file and input_folder and not patho_inference:
        taxon_file_pathname = os.path.join(input_folder, 'taxon_id.tsv')
        if os.path.exists(taxon_file_pathname):
            sys.exit('taxon ID file (' + taxon_file_pathname + ') already exists.')
        else:
            taxon_ids = read_taxon_id(input_folder)
            with open(taxon_file_pathname, 'w', encoding='utf-8') as taxon_id_file:
                taxon_id_writer = csv.writer(taxon_id_file, delimiter='\t')
                taxon_id_writer.writerow(['species', 'taxon_id'])
                for species, taxon_id in taxon_ids.items():
                    taxon_id_writer.writerow([species, taxon_id])

    # Turn off loading of pubmed entries.
    if no_download_articles:
        utils.pubmed_citations(activate_citations=False)

    # Modify pathway prediction score.
    if pathway_score:
        utils.modify_pathway_score(pathway_score)

    # Check input folder and create input files for PathoLogic.
    if input_folder:
        run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
        run_patho_flat_ids, run_flat_ids = check_input_and_existing_pgdb(run_ids, input_folder, output_folder, number_cpu_to_use)

        # Launch PathoLogic inference on species with no PGDBs.
        if run_patho_flat_ids:
            # Create the list containing all the data used by the multiprocessing call.
            multiprocess_pwt_input_files = []
            multiprocess_run_pwts = []
            multiprocess_run_pwt_flats = []
            multiprocess_run_move_pgdbs = []
            for run_patho_flat_id in run_patho_flat_ids:
                input_folder_path = os.path.join(input_folder, run_patho_flat_id)
                species_pgdb_folder = os.path.join(pgdbs_folder_path, run_patho_flat_id.lower() + 'cyc')
                input_run_move_pgdbs = [run_patho_flat_id, species_pgdb_folder]
                input_run_move_pgdbs.extend([output_folder, dat_extraction, size_reduction, xml_extraction, owl_extraction, col_extraction])
                multiprocess_pwt_input_files.append([input_folder_path, taxon_file])
                multiprocess_run_pwts.append([input_folder_path, patho_hole_filler, patho_operon_predictor, patho_transporter_inference])
                multiprocess_run_pwt_flats.append([input_folder_path])
                multiprocess_run_move_pgdbs.append(input_run_move_pgdbs)

            logger.info('~~~~~~~~~~Creation of input data from Genbank/GFF/PF~~~~~~~~~~')
            input_error_status = mpwt_pool.starmap(pwt_input_files, multiprocess_pwt_input_files)
            if any(input_error_status):
                close_mpwt(mpwt_pool, no_download_articles, pathway_score)
                sys.exit('Error during PathoLogic input files creation.')

            input_time = time.time()
            times.append(input_time)
            steps.append('pwt input creation')
            logger.info('----------End of creation of input data from Genbank/GFF/PF: {0:.2f}s----------'.format(times[-1] - times[-2]))

            # Launch PathoLogic.
            if patho_inference:
                logger.info('~~~~~~~~~~Inference on the data~~~~~~~~~~')
                error_status = mpwt_pool.starmap(run_pwt, multiprocess_run_pwts)

                # Check PathoLogic build.
                logger.info('~~~~~~~~~~Check inference~~~~~~~~~~')
                passed_inferences = check_pwt(multiprocess_run_pwts, patho_log)
                if any(error_status):
                    if ignore_error:
                        logger.critical('Error during inference. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')
                    else:
                        close_mpwt(mpwt_pool, no_download_articles, pathway_score)
                        sys.exit('Error during inference. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')

            patho_time = time.time()
            times.append(patho_time)
            steps.append('PathoLogic inference')
            logger.info('----------End of PathoLogic inference: {0:.2f}s----------'.format(times[-1] - times[-2]))
        else:
            multiprocess_run_pwt_flats = []
            multiprocess_run_move_pgdbs = []
            passed_inferences = []

    # Create path for lisp if there is no folder given.
    # Create the input for the creation of BioPAX/attribute-values files.
    if (flat_creation and not input_folder) or (output_folder and not input_folder):
        # Create a temporary folder in ptools-local where lisp script will be stored.
        tmp_folder = os.path.join(ptools_local_path, 'tmp')
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)

        # Create a lisp script file for each PGDB in the ptools-local folder.
        flat_run_ids = create_only_flat_lisp(pgdbs_folder_path, tmp_folder)
        multiprocess_run_pwt_flats = []
        multiprocess_run_move_pgdbs = []
        for flat_run_id in flat_run_ids:
            input_tmp_folder_path = os.path.join(tmp_folder, flat_run_id)
            species_pgdb_folder = os.path.join(pgdbs_folder_path, flat_run_id.lower() + 'cyc')
            input_run_move_pgdbs = [flat_run_id, species_pgdb_folder]

            input_run_move_pgdbs.extend([output_folder, dat_extraction, size_reduction, xml_extraction, owl_extraction, col_extraction])
            multiprocess_run_pwt_flats.append([input_tmp_folder_path])
            multiprocess_run_move_pgdbs.append(input_run_move_pgdbs)

    # Add species that have data in PGDB but are not present in output folder.
    # Or if ignore_error has been used, select only PathoLogic build that have succeed + species in input with PGDB and not in output.
    if input_folder:
        if ignore_error:
            multiprocess_run_pwt_flats = []
            multiprocess_run_move_pgdbs = []
            if run_patho_flat_ids:
                if passed_inferences:
                    tmp_run_flat_ids = list(set(passed_inferences).intersection(set(run_patho_flat_ids)))
                else:
                    tmp_run_flat_ids = []
            else:
                tmp_run_flat_ids = []
            if run_flat_ids:
                tmp_run_flat_ids.extend(run_flat_ids)
            run_flat_ids = tmp_run_flat_ids
        if run_flat_ids:
            for run_flat_id in run_flat_ids:
                flat_creation_path = os.path.join(*[input_folder, run_flat_id, 'flat_files_creation.lisp'])
                create_flat_creation_script(run_flat_id, flat_creation_path)
                input_folder_path = os.path.join(input_folder, run_flat_id)
                species_pgdb_folder = os.path.join(pgdbs_folder_path, run_flat_id.lower() + 'cyc')
                multiprocess_run_pwt_flats.append([input_folder_path])
                input_run_move_pgdbs = [run_flat_id, species_pgdb_folder]
                input_run_move_pgdbs.extend([output_folder, dat_extraction, size_reduction, xml_extraction, owl_extraction, col_extraction])
                multiprocess_run_move_pgdbs.append(input_run_move_pgdbs)

    if not multiprocess_run_pwt_flats:
        close_mpwt(mpwt_pool, no_download_articles, pathway_score)
        logger.critical('No PGDB to export to move to output folder.')
        return

    if not multiprocess_run_move_pgdbs:
        close_mpwt(mpwt_pool, no_download_articles, pathway_score)
        logger.critical('No PGDB to export in flat files or to move to output folder.')
        return

    # Create BioPAX/attributes-values flat files.
    if (input_folder and flat_creation) or flat_creation:
        logger.info('~~~~~~~~~~Creation of the flat files~~~~~~~~~~')
        flat_error_status = mpwt_pool.starmap(run_pwt_flat, multiprocess_run_pwt_flats)
        logger.info('~~~~~~~~~~Check .dat~~~~~~~~~~')
        for multiprocess_run_move_pgdb in multiprocess_run_move_pgdbs:
            check_dat(multiprocess_run_move_pgdb[0], multiprocess_run_move_pgdb[1])
        if any(flat_error_status):
            if ignore_error:
                logger.critical('Error during flat files creation. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')
            else:
                close_mpwt(mpwt_pool, no_download_articles, pathway_score)
                sys.exit('Error during flat files creation. Process stopped. Look at the command log. Also by using --log argument, you can have additional information.')

        flat_time = time.time()
        times.append(flat_time)
        steps.append('BioPAX/attribute-value flat files creation')
        logger.info('----------End of flat files creation: {0:.2f}s----------'.format(times[-1] - times[-2]))

    if (flat_creation and not input_folder) or (output_folder and not input_folder):
        ptools_local_path = utils.find_ptools_path()
        ptools_local_tmp_path = os.path.join(ptools_local_path, 'tmp')
        shutil.rmtree(ptools_local_tmp_path)

    logger.info('~~~~~~~~~~End of Pathway Tools~~~~~~~~~~')

    # Move PGDBs or attribute-values/dat files.
    if output_folder:
        logger.info('~~~~~~~~~~Moving result files~~~~~~~~~~')
        if not os.path.exists(output_folder):
            logger.info('No output directory, it will be created.')
            os.mkdir(output_folder)
        mpwt_pool.starmap(run_move_pgdb, multiprocess_run_move_pgdbs)

        move_time = time.time()
        times.append(move_time)
        steps.append('Moving results files')
        logger.info('----------End of moving files: {0:.2f}s----------'.format(times[-1] - times[-2]))


    close_mpwt(mpwt_pool, no_download_articles, pathway_score)

    end_time = time.time()
    times.append(end_time)
    steps.append('mpwt')

    # Write each step time in log file.
    if patho_log:
        if patho_log:
            if not os.path.exists(patho_log):
                logger.info('No log directory, it will be created.')
                os.mkdir(patho_log)
        patho_error_pathname = os.path.join(patho_log, 'log_error.txt')
        with open(patho_error_pathname, 'a') as input_file:
            input_file.write('\n\n---------Time---------\n')
            for index, step_time in enumerate(times):
                if index != 0:
                    if index + 1 == len(times):
                        step_duration = step_time - times[0]
                    else:
                        step_duration = step_time - times[index-1]
                    input_file.write('Step {0} takes: {1:.2f}s.\n'.format(steps[index] , step_duration))

    logger.info('----------mpwt has finished in {0:.2f}s! Thank you for using it.'.format(end_time - start_time))


def give_permission(folder, permission):
    """ Give permission to group or all for folder.

    Args:
        folder (str): pathname to the input folder
        permission (str): level of permission (either group or all)
    """
    if permission == 'group':
        for dirpath, _dirnames, filenames in os.walk(folder):
            os.chmod(dirpath, stat.S_IRGRP | ~stat.S_IWGRP | ~stat.S_IXGRP)
            for filename in filenames:
                os.chmod(os.path.join(dirpath, filename), ~stat.S_IRGRP | ~stat.S_IWGRP | ~stat.S_IXGRP)
    elif permission == 'all':
        for dirpath, _dirnames, filenames in os.walk(folder):
            os.chmod(dirpath, stat.S_IROTH | stat.S_IWOTH | stat.S_IXOTH)
            for filename in filenames:
                os.chmod(os.path.join(dirpath, filename), stat.S_IROTH | stat.S_IWOTH | stat.S_IXOTH)
    else:
        logger.critical(f'Invalid permission "{permission}", permission must be "group" or "all"')


def run_mpwt(run_folder=None, input_folder=None, run_input_files_creation=None,
                run_output_folder=None, output_folder=None,
                run_patho_inference=None, pathologic_options=None,
                run_flat_creation=None, move_options=None,
                taxon_file=None, permission=None):
    """ Single run of mpwt on one folder.
    Used in multiprocessing in independent_mpwt.

    Args:
        run_folder (str): name of the folder containg input files
        input_folder (str): pathname to input folder
        run_input_files_creation (bool): if True runs creation of PathoLogic input files
        run_output_folder (bool): if True moves the output file to the ouput folder
        output_folder (str): pathname to output folder
        run_patho_inference (bool): if True PathoLogic is run on the input folder
        pathologic_options (list): list of bool for: patho_hole_filler, patho_operon_predictor, patho_transporter_inference
        run_flat_creation (bool): if True flat files will be created
        move_options (list): list of bool for: dat_extraction, size_reduction, xml_extraction, owl_extraction, col_extraction
        taxon_file (str): pathname to the mpwt taxon ID file
        permission (str): Choose permission access to PGDB in ptools-local and output files, either 'all' or 'group' (by default it is user).
    """
    ptools_local_path = utils.find_ptools_path()
    pgdbs_folder_path = os.path.join(*[ptools_local_path, 'pgdbs', 'user'])
    species_pgdb_folder = os.path.join(pgdbs_folder_path, run_folder.lower() + 'cyc')

    patho_error_status = False
    flat_error_status = False
    move_error_status = False

    if input_folder:
        run_folder_path = os.path.join(input_folder, run_folder)

    if run_input_files_creation:
        pwt_input_files(run_folder_path, taxon_file)

    if run_patho_inference:
        patho_error_status = run_pwt(run_folder_path, *pathologic_options)
        if patho_error_status:
            return run_folder, patho_error_status, flat_error_status, move_error_status

    if run_flat_creation:
        flat_error_status = run_pwt_flat(run_folder_path)
        check_dat(run_folder_path, species_pgdb_folder)
        if flat_error_status:
            return run_folder, patho_error_status, flat_error_status, move_error_status

    if permission:
        give_permission(permission, species_pgdb_folder)

    if run_output_folder:
        move_error_status = run_move_pgdb(run_folder, species_pgdb_folder, output_folder, *move_options)
        if move_error_status:
            return run_folder, patho_error_status, flat_error_status, move_error_status

    if permission and output_folder:
        give_permission(permission, output_folder)

    return run_folder, patho_error_status, flat_error_status, move_error_status


def independent_mpwt(input_folder, output_folder=None, patho_inference=None,
                     patho_hole_filler=None, patho_operon_predictor=None,
                     patho_transporter_inference=None, no_download_articles=None,
                     flat_creation=None, dat_extraction=None, xml_extraction=None,
                     owl_extraction=None, col_extraction=None, size_reduction=None,
                     number_cpu_to_use=None, patho_log=None, pathway_score=None,
                     taxon_file=None, permission=None):
    """
    Function managing the workflow for independent run of mpwt.
    Each process of Pathway Tools on an organism are run separatly so if one failed the other that passed will succeed.

    Args:
        input_folder (str): pathname to input folder
        output_folder (str): pathname to output folder
        patho_inference (bool): PathoLogic inference (True/False)
        patho_hole_filler (bool): PathoLogic hole filler (True/False)
        patho_operon_predictor (bool): PathoLogic operon predictor (True/False)
        patho_transporter_inference (bool): PathoLogic Transport Inference Parser (True/False)
        no_download_articles (bool): turning off loading of PubMed citations (True/False)
        flat_creation (bool): BioPAX/attributes-values flat files creation (True/False)
        dat_extraction (bool): move BioPAX/attributes-values files to output folder (True/False)
        xml_extraction (bool): move metabolic-reactions.xml to output folder (True/False)
        owl_extraction (bool): move owl files to output folder (True/False)
        col_extraction (bool): move tabular files to output folder (True/False)
        size_reduction (bool): delete ptools-local data at the end (True/False)
        number_cpu (int): number of CPU used (default=1)
        patho_log (str): pathname to mpwt log folder
        pathway_score (float): score between 0 and 1 to accept or reject pathway
        taxon_file (str): pathname to the mpwt taxon ID file
        permission (str): Choose permission access to PGDB in ptools-local and output files, either 'all' or 'group' (by default it is user).
    """
    logger.info('---------- Launching mpwt ----------')
    ptools_local_path = utils.find_ptools_path()
    pgdbs_folder_path = os.path.join(*[ptools_local_path, 'pgdbs', 'user'])

    start_time = time.time()

    # Check if input folder exists and is a folder.
    if input_folder:
        if not os.path.exists(input_folder):
            logger.critical('mpwt can not run: ' + input_folder + ' does not exist.')
            return
        if not os.path.isdir(input_folder):
            logger.critical('mpwt can not run: ' + input_folder + ' is not a directory.')
            return

    # If output_folder does not exists, creates it.
    if output_folder:
        if not os.path.exists(output_folder):
            logger.info('No output directory, it will be created.')
            os.mkdir(output_folder)

    # Turn off loading of pubmed entries.
    if no_download_articles:
        utils.pubmed_citations(activate_citations=False)

    # Modify pathway prediction score.
    if pathway_score:
        old_pathway_score = utils.extract_pathway_score()
        utils.modify_pathway_score(pathway_score)

    mpwt_pool = Pool(processes=number_cpu_to_use)

    # Check input subfolder.
    if input_folder:
        run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
        run_patho_flat_ids, run_flat_ids = check_input_and_existing_pgdb(run_ids, input_folder, output_folder, number_cpu_to_use)

    # Create path for lisp if there is no folder given.
    # Create the input for the creation of BioPAX/attribute-values files.
    if (flat_creation and not input_folder) or (output_folder and not input_folder):
        # Create a temporary folder in ptools-local where lisp script will be stored.
        tmp_folder = os.path.join(ptools_local_path, 'tmp')
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)

        # Create a lisp script file for each PGDB in the ptools-local folder.
        run_ids = list(create_only_flat_lisp(pgdbs_folder_path, tmp_folder))
        if flat_creation:
            run_flat_ids = list(run_ids)
        else:
            run_flat_ids = None
        run_patho_flat_ids = None

    pathologic_options = [patho_hole_filler, patho_operon_predictor, patho_transporter_inference]
    move_options = [dat_extraction, size_reduction, xml_extraction, owl_extraction, col_extraction]

    # Create data for multiprocessing.
    # For each organism, find if a PathoLogic inference must be run, if a flat creation must be performed and if the output files must be moved.
    multiprocess_run_mpwts = []
    run_input_files_creation = False
    run_patho_inference = False
    run_flat_creation = False
    run_output_folder = False
    tmp_folder = False

    for run_id in run_ids:
        # For species without PGDB in ptools-local, launch input files creations, PathoLogic reconstruction, flat files creation and moving output files (according to user input)
        if run_patho_flat_ids and run_id in run_patho_flat_ids:
            run_input_files_creation = True
            if patho_inference:
                run_patho_inference = True
            if flat_creation:
                run_flat_creation = True
            if output_folder:
                run_output_folder = True
        # For speccies with PGDB in ptools-local, launch only flat files creation and moving output files (according to user input)
        if run_flat_ids and run_id in run_flat_ids:
            if flat_creation:
                run_flat_creation = True
            if output_folder:
                run_output_folder = True

        # If no input_folder, PGDBs from ptools-local will be used.
        if not input_folder:
            tmp_folder = True
            # If flat_creation, flat files of these PGDBs will be created and moved to the output folder.
            if flat_creation:
                run_flat_creation = True
                input_folder = os.path.join(ptools_local_path, 'tmp')
            if output_folder:
                if not os.path.exists(os.path.join(output_folder, run_id)):
                    run_output_folder = True
                else:
                    logger.info('{0} contains already {1}, output files will not be moved'.format(output_folder, run_id))

        multiprocess_run_mpwt = [run_id, input_folder, run_input_files_creation, run_output_folder, output_folder, run_patho_inference, pathologic_options,
                                run_flat_creation, move_options, taxon_file, permission]

        multiprocess_run_mpwts.append(multiprocess_run_mpwt)

    results = mpwt_pool.starmap(run_mpwt, multiprocess_run_mpwts)

    if patho_log:
        if not os.path.exists(patho_log):
            os.mkdir(patho_log)
        input_folders = [os.path.join(multiprocess_run_mpwt[1], multiprocess_run_mpwt[0]) for multiprocess_run_mpwt in multiprocess_run_mpwts]
        check_mpwt_pathologic_runs(input_folders, patho_log)

    logger.info('-------------- Checking mpwt runs --------------')
    nb_failed_runs = len([result for result in results if any(result[1:])])
    nb_total_runs = len(results)
    if nb_failed_runs > 0:
        if nb_failed_runs == 1:
            failed_str = '{0} failed run'.format(nb_failed_runs)
        else:
            failed_str = '{0} failed runs'.format(nb_failed_runs)
        if nb_total_runs == 1:
            total_str = '{0} run'.format(nb_total_runs)
        else:
            total_str = '{0} runs'.format(nb_total_runs)

        logger.info('/!\\ {0} on a total of {1}.'.format(failed_str, total_str))
    else:
        logger.info('All runs are successful.')

    for result in results:
        run_id = result[0]
        if any(result[1:]):
            if result[1]:
                logger.info('/!\\ Error in {0} during PathoLogic inference step.'.format(run_id))
            if result[2]:
                logger.info('/!\\ Error in {0} during Flat files creation step.'.format(run_id))
            if result[3]:
                logger.info('/!\\ Error in {0} during Moving output files step.'.format(run_id))

    # Remove tmp folder in ptools-local.
    if tmp_folder:
        ptools_local_tmp_path = os.path.join(ptools_local_path, 'tmp')
        shutil.rmtree(ptools_local_tmp_path)

    # Close multiprocessing Pool
    if pathway_score:
        close_mpwt(mpwt_pool, no_download_articles, pathway_score, old_pathway_score)
    else:
        close_mpwt(mpwt_pool, no_download_articles)

    end_time = time.time()

    logger.info('-------------- mpwt has finished in {0:.2f}s! Thank you for using it. --------------'.format(end_time - start_time))

