# Copyright (C) 2018-2022 Arnaud Belcour- Inria Dyliss
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
Uselful functions for mpwt.
"""

import logging
import os
import shutil
import sys

from multiprocessing import Pool

logger = logging.getLogger(__name__)


def find_ptools_path():
    """
    Check if Pathway Tools is in PATH. If not stop the script.
    If it finds it, return the path of ptools-local by reading Pathway Tools file.

    Args:
    Returns:
        string: pathname to ptools-local folder
    """
    pathway_tools_path = shutil.which('pathway-tools')
    if not pathway_tools_path:
        logger.critical('Pathway Tools is not in the Path, mpwt can not work without it.')
        sys.exit(1)

    pathway_tools_file = open(pathway_tools_path, 'r')
    ptools_local_str = [line for line in pathway_tools_file if 'PTOOLS_LOCAL_PATH' in line][0]
    ptools_local_path = os.path.join(ptools_local_str.split(';')[0].split('=')[1].replace('"', '').strip(' '), 'ptools-local')
    pathway_tools_file.close()

    return ptools_local_path


def check_ptools_local_pwt():
    """
    Check if ptools-init.dat exists in ptools-local folder.

    Args:
    Returns:
        bool: True if ptools-init.dat is missing
    """
    ptools_path = find_ptools_path()

    error = None

    ptools_init_path = os.path.join(ptools_path, 'ptools-init.dat')

    if not os.path.exists(ptools_init_path):
        logger.critical('Missing ptools-init.dat file in ptools-local folder. Use "pathway-tools -config" to recreate it.')
        error = True

    return error


def list_pgdb():
    """
    List all the PGDB inside the ptools-local folder.
    Return a list of their IDs.

    Args:
    Returns:
        list: PGDB IDs inside ptools-local
    """
    ptools_local_path = find_ptools_path()
    pgdb_folder = os.path.join(*[ptools_local_path, 'pgdbs', 'user'])

    return [species_pgdb for species_pgdb in os.listdir(pgdb_folder) if 'cyc' in species_pgdb]


def delete_pgdb(pgdb_name):
    """
    Remove a specific PGDB.

    Args:
        pgdb_name (str): PGDB ID to delete
    """
    ptools_local_path = find_ptools_path()
    pgdb_path = os.path.join(*[ptools_local_path.replace('\n', ''), 'pgdbs', 'user', pgdb_name])
    if os.path.isdir(pgdb_path):
        shutil.rmtree(pgdb_path)
        logger.info('{0} (at {1}) has been removed.'.format(pgdb_name, pgdb_path))


def remove_pgdbs(to_delete_pgdbs, number_cpu=None):
    """
    Delete all PGDB inside to_delete_pgdbs using multiprocessing.
    Check if there is a Pool and if not spawn one.

    Args:
        to_delete_pgdbs (list): PGDB IDs to delete
        number_cpu (str): number of CPUs to use (default = 1)
    """
    # Use the number of cpu given by the user or all the cpu available.
    if number_cpu:
        number_cpu_to_use = int(number_cpu)
    else:
        number_cpu_to_use = 1
    mpwt_pool = Pool(processes=number_cpu_to_use)

    mpwt_pool.map(delete_pgdb, to_delete_pgdbs)

    mpwt_pool.close()
    mpwt_pool.join()


def cleaning(number_cpu=None, verbose=None):
    """
    Clean Pathway Tools PGDB's folder.
    The script will delete folders and files in ptools-local/pgdbs/user.

    Args:
        number_cpu (str): number of CPUs to use
        verbose (bool): verbose argument
    """
    if verbose:
        logger.setLevel(logging.DEBUG)

    ptools_local_path = find_ptools_path()
    file_path = os.path.join(*[ptools_local_path.replace('\n', ''), 'pgdbs', 'user'])

    pgdb_metadata_path = os.path.join(file_path, 'PGDB-METADATA.ocelot')
    if os.path.isfile(pgdb_metadata_path):
        os.remove(pgdb_metadata_path)
        logger.info('PGDB-METADATA.ocelot has been removed.')

    pgdb_counter_path = os.path.join(file_path, 'PGDB-counter.dat')
    if os.path.isfile(pgdb_counter_path):
        os.remove(pgdb_counter_path)
        logger.info('PGDB-counter.dat has been removed.')

    # Extract all pgdbs inside ptools-local. Then delete them.
    all_pgdbs = os.listdir(file_path)
    remove_pgdbs(all_pgdbs, number_cpu)


def cleaning_input(input_folder, verbose=None):
    """
    Remove flat_files_creation.lisp, pathologic.log, genetic-elements.dat and organism-params.dat in a genbank folder.

    Args:
        input_folder (str): pathname to input folder
        verbose (bool): verbose argument
    """
    if verbose:
        logger.setLevel(logging.DEBUG)

    if not os.path.exists(input_folder):
        sys.exit('mpwt can not run: ' + input_folder + ' does not exist.')
    if not os.path.isdir(input_folder):
        sys.exit('mpwt can not run: ' + input_folder + ' is not a directory.')

    run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]

    input_paths = [os.path.join(input_folder, run_id) for run_id in run_ids]

    for input_path in input_paths:
        if os.path.isdir(input_path):
            lisp_script = os.path.join(input_path, 'flat_files_creation.lisp')
            patho_log = os.path.join(input_path, 'pathologic.log')
            pwt_log = os.path.join(input_path, 'pwt_terminal.log')
            flat_log = os.path.join(input_path, 'flat_files_creation.log')
            genetic_dat = os.path.join(input_path, 'genetic-elements.dat')
            organism_dat = os.path.join(input_path, 'organism-params.dat')
            if os.path.exists(lisp_script):
                os.remove(lisp_script)
            if os.path.exists(patho_log):
                os.remove(patho_log)
            if os.path.exists(pwt_log):
                os.remove(pwt_log)
            if os.path.exists(flat_log):
                os.remove(flat_log)
            if os.path.exists(genetic_dat):
                os.remove(genetic_dat)
            if os.path.exists(organism_dat):
                os.remove(organism_dat)
            species = os.path.basename(input_path)
            logger.info('Remove ' + species + ' temporary datas.')


def pubmed_citations(activate_citations):
    """
    Activate or deactivate loading of PubMed citations.

    Args:
        activate_citations (bool): boolean to indicate if you want to activate or not the downlaod of Pubmed entries.
    """
    ptools_init_filepath = os.path.join(find_ptools_path(), 'ptools-init.dat')
    new_ptools_file = ""

    download_pubmed_entries_parameter = None
    with open(ptools_init_filepath, 'r') as ptools_init_file:
        for line in ptools_init_file.read().split('\n'):
            if 'Batch-PathoLogic-Download-Pubmed-Entries?' in line:
                if '#' in line:
                    line = line.replace('#', '')
                download_pubmed_entries_parameter = True
                if activate_citations:
                    line = line.replace('NIL', 'T')
                else:
                    line = line.replace('T', 'NIL')
            if line != '':
                new_ptools_file = new_ptools_file + line + '\n'
            else:
                new_ptools_file = new_ptools_file + line

    if not download_pubmed_entries_parameter:
        sys.exit('There is no Batch-PathoLogic-Download-Pubmed-Entries parameter in ' + ptools_init_filepath +'. To use --nc/no_download_articles, mpwt needs Pathway Tools 23.5 or higher.')

    with open(ptools_init_filepath, 'w', encoding='utf-8') as ptools_init_file:
        ptools_init_file.write(new_ptools_file)


def extract_pathway_score():
    """
    Get the Pathway-Prediction-Score-Cutoff of ptools-init.dat

    Returns:
        pathway_score (float): score between 0 and 1 to accept or reject pathways
    """
    ptools_init_filepath = os.path.join(find_ptools_path() ,'ptools-init.dat')

    pathway_prediction_score_cutoff = None
    with open(ptools_init_filepath, 'r') as ptools_init_file:
        for line in ptools_init_file.read().split('\n'):
            if 'Pathway-Prediction-Score-Cutoff' in line:
                pathway_prediction_score_cutoff = line.split(' ')[1]

    if not pathway_prediction_score_cutoff:
        sys.exit('There is no Pathway-Prediction-Score-Cutoff parameter in ' + ptools_init_filepath +'.')

    return pathway_prediction_score_cutoff


def modify_pathway_score(pathway_score, comment_line=None):
    """
    Modify the Pathway-Prediction-Score-Cutoff of ptools-init.dat

    Args:
        pathway_score (float): score between 0 and 1 to accept or reject pathways
        comment_line (bool): boolean if True comment Pathway-Prediction-Score-Cutoff line
    """
    ptools_init_filepath = os.path.join(find_ptools_path() ,'ptools-init.dat')
    new_ptools_file = ""

    pathway_prediction_score_cutoff = None
    with open(ptools_init_filepath, 'r') as ptools_init_file:
        for line in ptools_init_file.read().split('\n'):
            if 'Pathway-Prediction-Score-Cutoff' in line:
                pathway_prediction_score_cutoff = True
                if '#' in line:
                    line = line.replace('#', '')
                if pathway_score:
                    if comment_line:
                        line = '###' + line.split(' ')[0] + ' ' + str(pathway_score)
                    else:
                        line = line.split(' ')[0] + ' ' + str(pathway_score)
            if line != '':
                new_ptools_file = new_ptools_file + line + '\n'
            else:
                new_ptools_file = new_ptools_file + line

    if not pathway_prediction_score_cutoff:
        sys.exit('There is no Pathway-Prediction-Score-Cutoff parameter in ' + ptools_init_filepath +'.')

    with open(ptools_init_filepath, 'w', encoding='utf-8') as ptools_init_file:
        ptools_init_file.write(new_ptools_file)