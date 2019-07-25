"""
Check results from Pathway Tools command:
-check_pwt: results from PathoLogic (by looking at pathologic.log)
-check_dat: attribute-values dat files
"""

import csv
import os
import logging
import sys

from mpwt.utils import permission_change

logger = logging.getLogger(__name__)


def check_pwt(multiprocess_inputs, patho_log_folder):
    """
    Check PathoLogic's log.
    Create two log files (log_error.txt which contains Pathway Tools log and resume_inference.tsv which contains summary of metabolic networks).

    Args:
        multiprocess_inputs (list): list of dictionary contaning multiprocess input data
        patho_log_folder (str): pathname to the PathoLogic log folder.

    Returns:
        list: Species with successful build.
    """
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
                        if species not in failed_inferences:
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
        logger.info('\n{0} {1} passed!\n'.format(str(number_passed_inference), string_passed_build))
    if number_failed_inference > 0:
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

    if patho_log_folder:
        permission_change(patho_log_folder)

    return passed_inferences

def check_dat(multiprocess_input):
    """
    Check dats creation.

    Args:
        multiprocess_input (dictionary): contains multiprocess input (mpwt argument: input folder, output folder, ...)
    """
    pgdb_folder = multiprocess_input['pgdb_folders']
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

    expected_dat_number = str(len(dat_files))
    found_dat_number = str(len(dat_checks))
    logger.info('{0}: {1} out of {2} dat files create.'.format(pgdb_folder_dbname, found_dat_number, expected_dat_number))

