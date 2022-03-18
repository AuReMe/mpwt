# Copyright (C) 2018-2022 Arnaud Belcour - Inria Dyliss
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
Check results from Pathway Tools command:
-check_mpwt_pathologic_runs: results from PathoLogic (by looking at pathologic.log)
-check_dat: attribute-values dat files
"""

import csv
import os
import logging
import re

logger = logging.getLogger(__name__)


def extract_pathologic(patho_log):
    """ Read PathoLogic's log and extract informations such as if the build has passed, if there is warning/error.
    And also if the build has passed the number of gene/protein/reaction/compound in the draft network.

    Args:
        patho_log (str): pathname to the pathologic.log file.

    Returns:
        non_fatal_error_count (int): number of non fatal errors.
        warning_count (int): number of warning message.
        fatal_error_index (None, int): either None (if no fatal error) or the index of the line in pathologic.log with fatal error
        passed_inferences (None, bool): either None (if inference has no passed) or True otherwise
        pgdb_build_done (None, bool): either None (if inference has no passed) or True otherwise
        gene_number (None, int): either None (if build has failed) or int with the number of genes in network
        protein_number (None, int): either None (if build has failed) or int with the number of proteins in network
        pathway_number (None, int): either None (if build has failed) or int with the number of pathways in network
        reaction_number (None, int): either None (if build has failed) or int with the number of reactions in network
        compound_number (None, int): either None (if build has failed) or int with the number of compounds in network
        log_str (str): string containing informaitons about the PathoLogic build behaves
        log_resume_list (list): list with number of gene/protein/reaction/compound if build has passed
    """
    non_fatal_error_count = 0
    warning_count = 0
    fatal_error_index = None
    passed_inferences = False
    pgdb_build_done = False
    gene_number = None
    protein_number = None
    pathway_number = None
    reaction_number = None
    compound_number = None

    log_str = ''
    log_resume_list = []

    input_folder = os.path.dirname(patho_log)
    organism_name = os.path.basename(input_folder)

    log_str += '------------ Species: '
    log_str += organism_name + '\n'

    if not os.path.exists(patho_log):
        log_str += 'No pathologic log, an error occured before PathoLogic run.\n'
        log_resume_list.append([organism_name, 'ERROR', '', '', '', ''])
        logger.info('No pathologic log for {0}, an error occured before PathoLogic run.'.format(organism_name))
        return

    with open(patho_log, 'r') as input_file:
        for index, line in enumerate(input_file):
            if ';;; Error:' in line:
                non_fatal_error_count += 1
            if 'Warning:' in line:
                warning_count += 1
            if 'fatal error' in line:
                fatal_error_index = index
                log_str += line
                log_resume_list.append([organism_name, 'ERROR', '', '', '', ''])

            if fatal_error_index:
                if index > fatal_error_index:
                    log_str += line

            # Search for Build done line and its following line which look like (for Pathway Tools 25.0):
            # PGDB contains XXXX genes, XXXX polypeptides, XXXX base pathways, XXXX reactions, XXXX compounds
            # In Pathway Tools 25.5 it looks like this:
            # PGDB contains XXXX classes and XXXX instances: XXXX genes, XXXX polypeptides, XXXX base pathways, XXXX reactions, XXXX compounds, XXXX publications
            if 'Build done.' in line or 'PGDB build done.' in line:
                log_str += line
                if non_fatal_error_count > 0:
                    log_str += 'Number of non fatal errors: ' + str(non_fatal_error_count) + '. More information in ' + patho_log + '.\n'
                if warning_count > 0:
                    log_str += 'Number of warning: ' + str(warning_count) + '. More information in ' + patho_log + '.\n'

                resume_inference_line = next(input_file)
                log_str += resume_inference_line
                pgdb_build_done = True
                # Search the PGDB stat line and use regex to extract informations.
                # # This is done by searching for association like (digit word) like(XXXX genes).
                pgdb_stats = {}
                pgdb_line_re = r'(?P<stat_nb>[\d]+)\ (?P<variable_name>[\w]+(\ pathways)?)'
                # Create a dictionary containing the stat  as: {'genes': XXXX, 'reactions': XXXX, ...}
                for match_found in re.finditer(pgdb_line_re, resume_inference_line):
                    pgdb_stats[match_found.group('variable_name')] = int(match_found.group('stat_nb'))

                gene_number = pgdb_stats['genes']
                # proteins is listed in pathologic.log for Pathway Tools inferior to 25.0
                # Since the 25.0 polypeptides replace proteins
                if 'proteins' in pgdb_stats:
                    protein_number = pgdb_stats['proteins']
                elif 'polypeptides' in pgdb_stats:
                    protein_number = pgdb_stats['polypeptides']
                pathway_number = pgdb_stats['base pathways']
                reaction_number = pgdb_stats['reactions']
                compound_number = pgdb_stats['compounds']

            if 'Done' in line:
                passed_inferences = True
                if pgdb_build_done:
                    log_resume_list.append([organism_name, gene_number, protein_number,
                                            pathway_number, reaction_number, compound_number,
                                            non_fatal_error_count, warning_count])

    if passed_inferences is None and fatal_error_index is None:
        log_str += 'No build in PathoLogic inference.'
        log_resume_list.append([organism_name, 'ERROR', '', '', '', ''])

    log_str += '------------\n\n'

    return organism_name, non_fatal_error_count, warning_count, fatal_error_index, \
            passed_inferences, pgdb_build_done, gene_number, \
            protein_number, pathway_number, reaction_number, \
            compound_number, log_str, log_resume_list


def check_mpwt_pathologic_runs(species_input_folder_paths, patho_log_folder):
    """
    Check PathoLogic's log.
    Create two log files (log_error.txt which contains Pathway Tools log and resume_inference.tsv which contains summary of metabolic networks).

    Args:
        species_input_folder_paths (list): list of input folder pathname
        patho_log_folder (str): pathname to the PathoLogic log folder.
    """
    mpwt_pathologic_informations = []

    failed_inferences = []
    passed_inferences = []
    no_pathologic_files = []
    for species_input_folder_path in species_input_folder_paths:
        patho_log = os.path.join(species_input_folder_path, 'pathologic.log')
        if os.path.exists(patho_log):
            species_pathologic_informations = extract_pathologic(patho_log)

            mpwt_pathologic_informations.append(species_pathologic_informations)
            if species_pathologic_informations[3] is not None:
                failed_inferences.append(species_pathologic_informations[0])
            elif species_pathologic_informations[4] is not None:
                passed_inferences.append(species_pathologic_informations[0])
            elif species_pathologic_informations[3] is None and species_pathologic_informations[4] is not None:
                failed_inferences.append(species_pathologic_informations[0])
        else:
            logger.info('|Output Check|WARNING: No pathologic.log file for {0}, could not write log.'.format(species_input_folder_path))
            base_name = os.path.basename(species_input_folder_path)
            no_pathologic_files.append(base_name)
            log_str = ''
            log_str += '------------ Species: '
            log_str += base_name + '\n'
            log_str += 'No pathologic.log\n'
            log_str += '------------\n\n'
            mpwt_pathologic_informations.append([base_name, *['']*10, log_str, [['No pathologic.log', '', '', '', '']]])


    number_passed_inference = len(passed_inferences)
    number_failed_inference = len(failed_inferences)
    number_total_build = number_passed_inference + number_failed_inference

    string_passed_build = 'build has' if number_passed_inference == 1 else 'builds have'
    string_failed_build = 'build has' if number_failed_inference == 1 else 'builds have'

    if number_passed_inference > 0:
        logger.info('|Output Check| {0} on {1} {2} passed!'.format(str(number_passed_inference), str(number_total_build), string_passed_build))
    if number_failed_inference > 0:
        logger.critical('|Output Check|WARNING: {0} on {1} {2} failed: {3}! See the log for more information.'.format(str(number_failed_inference), str(number_total_build), string_failed_build, ' '.join(failed_inferences)))

    if patho_log_folder:
        patho_error_pathname = os.path.join(patho_log_folder, 'log_error.txt')
        patho_resume_pathname = os.path.join(patho_log_folder, 'resume_inference.tsv')

        patho_log_file = open(patho_error_pathname, 'w', encoding='utf-8')
        patho_resume_file = open(patho_resume_pathname, 'w', encoding='utf-8')

        patho_resume_writer = csv.writer(patho_resume_file, delimiter='\t', lineterminator='\n')
        patho_resume_writer.writerow(['species', 'gene_number', 'protein_number', 'pathway_number', 'reaction_number', 'compound_number', 'pwt_non_fatal_error', 'pwt_warning'])

        if number_passed_inference > 0:
            patho_log_file.write('Build done: ' + str(number_passed_inference) + '\n')
            patho_log_file.write('Species: ' + ', '.join(passed_inferences) +  '\n\n')
        if number_failed_inference > 0:
            patho_log_file.write('Build failed: ' + str(number_failed_inference) + '\n')
            patho_log_file.write('Species: ' + ', '.join(failed_inferences) + '\n\n')

        for species_pathologic_informations in mpwt_pathologic_informations:
            patho_log_file.write(species_pathologic_informations[11])
            patho_resume_writer.writerow(*species_pathologic_informations[12])


def check_dat(run_dat_id, species_pgdb_folder):
    """
    Check dats creation.

    Args:
        run_dat_id (str): species ID
        species_pgdb_folder (str): path to species PGDB folder
    """
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
    logger.info('|Output Check|{0}| {1} out of {2} dat files created.'.format(run_dat_id, found_dat_number, expected_dat_number))
