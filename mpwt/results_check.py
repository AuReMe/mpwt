"""
Check results from Pathway Tools command:
-check_pwt: results from PathoLogic (by looking at pathologic.log)
-check_dat: attribute-values dat files
"""

import csv
import os
import logging
import sys

logger = logging.getLogger(__name__)


def check_pwt(multiprocess_run_pwts, patho_log_folder):
    """
    Check PathoLogic's log.
    Create two log files (log_error.txt which contains Pathway Tools log and resume_inference.tsv which contains summary of metabolic networks).

    Args:
        multiprocess_run_pwts (list): list of dictionary contaning multiprocess input data
        patho_log_folder (str): pathname to the PathoLogic log folder.

    Returns:
        list: Species with successful build.
    """
    if patho_log_folder:
        if not os.path.exists(patho_log_folder):
            logger.info('No log directory, it will be created.')
            os.mkdir(patho_log_folder)

        patho_error_pathname = os.path.join(patho_log_folder, 'log_error.txt')
        patho_resume_pathname = os.path.join(patho_log_folder, 'resume_inference.tsv')

        patho_error_file = open(patho_error_pathname, 'w', encoding='utf-8')
        patho_resume_file = open(patho_resume_pathname, 'w', encoding='utf-8')
        patho_resume_writer = csv.writer(patho_resume_file, delimiter='\t', lineterminator='\n')
        patho_resume_writer.writerow(['species', 'gene_number', 'protein_number', 'pathway_number', 'reaction_number', 'compound_number', 'pwt_non_fatal_error', 'pwt_warning'])

    failed_inferences = []
    passed_inferences = []

    for multiprocess_run_pwt in multiprocess_run_pwts:
        species_input_folder_path = multiprocess_run_pwt[0]
        species = os.path.basename(species_input_folder_path)
        patho_log = os.path.join(species_input_folder_path, 'pathologic.log')

        if patho_log_folder:
            patho_error_file.write('------------ Species: ')
            patho_error_file.write(species)
            patho_error_file.write('\n')

        fatal_error_index = None
        non_fatal_error_count = 0
        warning_count = 0

        if os.path.exists(patho_log):
            pgdb_build_done = False
            with open(patho_log, 'r') as input_file:
                for index, line in enumerate(input_file):
                    if ';;; Error:' in line:
                        non_fatal_error_count += 1
                    if 'Warning:' in line:
                        warning_count += 1
                    if 'fatal error' in line:
                        fatal_error_index = index
                        if species not in failed_inferences:
                            failed_inferences.append(species)
                        if patho_log_folder:
                            patho_error_file.write(line)
                            patho_resume_writer.writerow([species, 'ERROR', '', '', '', ''])

                    if fatal_error_index:
                        if index > fatal_error_index:
                            if patho_log_folder:
                                patho_error_file.write(line)

                    if 'Build done.' in line or 'PGDB build done.' in line:
                        pgdb_build_done = True
                        if patho_log_folder:
                            patho_error_file.write(line)
                            resume_inference_line = next(input_file)
                            patho_error_file.write(resume_inference_line)
                            if non_fatal_error_count > 0:
                                non_fatal_error_line = 'Number of non fatal errors: ' + str(non_fatal_error_count) + '. More information in ' + patho_log + '.\n'
                                patho_error_file.write(non_fatal_error_line)
                            if warning_count > 0:
                                warning_line = 'Number of warning: ' + str(warning_count) + '. More information in ' + patho_log + '.\n'
                                patho_error_file.write(warning_line)
                            gene_number = int(resume_inference_line.split('PGDB contains ')[1].split(' genes')[0])
                            protein_number = int(resume_inference_line.split('genes, ')[1].split(' proteins')[0])
                            pathway_number = int(resume_inference_line.split('proteins, ')[1].split(' base pathways')[0])
                            reaction_number = int(resume_inference_line.split('base pathways, ')[1].split(' reactions')[0])
                            compound_number = int(resume_inference_line.split('reactions, ')[1].split(' compounds')[0])

                    if 'Done' in line:
                        passed_inferences.append(species)
                        if patho_log_folder and pgdb_build_done:
                            patho_resume_writer.writerow([species, gene_number, protein_number, pathway_number, reaction_number, compound_number, non_fatal_error_count, warning_count])

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
        logger.info('\n{0} {1} passed!\n'.format(str(number_passed_inference), string_passed_build))
    if number_failed_inference > 0:
        logger.critical('WARNING: {0} {1} failed! See the log for more information.\n'.format(str(number_failed_inference), string_failed_build))

    if patho_log_folder:
        patho_error_file.close()
        patho_resume_file.close()
        with open(patho_error_pathname,'r') as contents:
            save = contents.read()
        with open(patho_error_pathname, 'w', encoding='utf-8') as output_file:
                output_file.write('Inference statistics:\n')
                if number_passed_inference > 0:
                    output_file.write('Build done: ' + str(number_passed_inference) + '\n')
                    output_file.write('Species: ' + ', '.join(passed_inferences) +  '\n\n')
                if number_failed_inference > 0:
                    output_file.write('Build failed: ' + str(number_failed_inference) + '\n')
                    output_file.write('Species: ' + ', '.join(failed_inferences) + '\n\n')
                output_file.write(save)

    return passed_inferences

def check_dat(run_dat_id, species_pgdb_folder):
    """
    Check dats creation.

    Args:
        run_dat_id (str): species ID
        species_pgdb_folder (str): path to species PGDB folder
    """
    pgdb_folder_dbname = run_dat_id.lower() + 'cyc'

    dats_path = os.path.join(*[species_pgdb_folder, '1.0', 'data'])

    dat_files = ["classes.dat", "compound-links.dat", "compounds.dat", "dnabindsites.dat", "enzrxns.dat", "gene-links.dat", "genes.dat", "pathway-links.dat",
                "pathways.dat", "promoters.dat", "protein-features.dat", "protein-links.dat", "proteins.dat", "protligandcplxes.dat", "pubs.dat",
                "reaction-links.dat", "reactions.dat", "regulation.dat", "regulons.dat", "rnas.dat", "species.dat", "terminators.dat", "transunits.dat"]

    dat_checks = []
    for dat_file in dat_files:
        dat_file_path = os.path.join(dats_path, dat_file)
        if os.path.exists(dat_file_path):
            dat_checks.append(dat_file_path)

    expected_dat_number = str(len(dat_files))
    found_dat_number = str(len(dat_checks))
    logger.info('{0}: {1} out of {2} dat files created.'.format(pgdb_folder_dbname, found_dat_number, expected_dat_number))

