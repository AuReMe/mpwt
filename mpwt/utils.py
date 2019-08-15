"""
Uselful functions for mpwt.
"""

import csv
import logging
import os
import shutil
import sys

from Bio import SeqIO
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
    ptools_local_path = ptools_local_str.split(';')[0].split('=')[1].replace('"', '').strip(' ') + '/ptools-local'
    pathway_tools_file.close()

    return ptools_local_path


def list_pgdb():
    """
    List all the PGDB inside the ptools-local folder.
    Return a list of their IDs.

    Args:
    Returns:
        list: PGDB IDs inside ptools-local
    """
    ptools_local_path = find_ptools_path()
    pgdb_folder = ptools_local_path + '/pgdbs/user/'

    return [species_pgdb for species_pgdb in os.listdir(pgdb_folder) if 'cyc' in species_pgdb]


def delete_pgdb(pgdb_name):
    """
    Remove a specific PGDB.

    Args:
        pgdb_name (str): PGDB ID to delete
    """
    ptools_local_path = find_ptools_path()
    pgdb_path = ptools_local_path.replace('\n', '') +'/pgdbs/user/' + pgdb_name
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
    file_path = ptools_local_path.replace('\n', '') +'/pgdbs/user/'

    pgdb_metadata_path = file_path + 'PGDB-METADATA.ocelot'
    if os.path.isfile(pgdb_metadata_path):
        os.remove(pgdb_metadata_path)
        logger.info('PGDB-METADATA.ocelot has been removed.')

    pgdb_counter_path = file_path + 'PGDB-counter.dat'
    if os.path.isfile(pgdb_counter_path):
        os.remove(pgdb_counter_path)
        logger.info('PGDB-counter.dat has been removed.')

    # Extract all pgdbs inside ptools-local. Then delete them.
    all_pgdbs = os.listdir(file_path)
    remove_pgdbs(all_pgdbs, number_cpu)


def cleaning_input(input_folder, verbose=None):
    """
    Remove dat_creation.lisp, pathologic.log, genetic-elements.dat and organism-params.dat in a genbank folder.

    Args:
        input_folder (str): pathname to input folder
        verbose (bool): verbose argument
    """
    if verbose:
        logger.setLevel(logging.DEBUG)

    run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]

    genbank_paths = [input_folder + "/" + run_id + "/" for run_id in run_ids]

    for genbank_path in genbank_paths:
        if os.path.isdir(genbank_path):
            lisp_script = genbank_path + 'dat_creation.lisp'
            patho_log = genbank_path + 'pathologic.log'
            genetic_dat = genbank_path + 'genetic-elements.dat'
            organism_dat = genbank_path + 'organism-params.dat'
            if os.path.exists(lisp_script):
                os.remove(lisp_script)
            if os.path.exists(patho_log):
                os.remove(patho_log)
            if os.path.exists(genetic_dat):
                os.remove(genetic_dat)
            if os.path.exists(organism_dat):
                os.remove(organism_dat)
            species = genbank_path.split('/')[-2]
            logger.info('Remove ' + species + ' temporary datas.')


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


def create_pathologic_file(input_folder, output_folder):

    for genbank_folder in os.listdir(input_folder):
        taxon_id = None
        input_path = input_folder + '/' + genbank_folder + '/' + genbank_folder + '.gbk'
        output_path = output_folder + '/' + genbank_folder
        if not os.path.exists(output_path):
            os.makedirs(output_path)

    # Add taxon ID in taxon_id.tsv if available.
    with open(input_path, "r") as gbk:
        first_seq_record = next(SeqIO.parse(gbk, "genbank"))
        src_features = [feature for feature in first_seq_record.features if feature.type == "source"]
        for src_feature in src_features:
            try:
                src_dbxref_qualifiers = src_feature.qualifiers['db_xref']
                for src_dbxref_qualifier in src_dbxref_qualifiers:
                    if 'taxon:' in src_dbxref_qualifier:
                        taxon_id = src_dbxref_qualifier.replace('taxon:', '')
            except KeyError:
                logger.info('No taxon ID in the Genbank {0} In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.'.format(genbank_folder))
        if taxon_id:
            if not os.path.exists(output_folder + '/taxon_id.tsv'):
                with open(output_folder + '/taxon_id.tsv', 'w') as taxon_id_file:
                    taxon_writer = csv.writer(taxon_id_file, delimiter='\t')
                    taxon_writer.writerow(['species', 'taxon_id'])
                    taxon_writer.writerow([genbank_folder, taxon_id])
            else:
                with open(output_folder + '/taxon_id.tsv', 'a') as taxon_id_file:
                    taxon_writer = csv.writer(taxon_id_file, delimiter='\t')
                    taxon_writer.writerow([genbank_folder, taxon_id])

        for record in SeqIO.parse(input_path, 'genbank'):
            element_id = record.id
            records = [record]
            SeqIO.write(records, output_path + '/' + element_id + '.fasta', 'fasta')
            with open(output_path + '/' + element_id + '.pf', 'w') as element_file:
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                element_file.write(';; ' + element_id + '\n')
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                for feature in record.features:
                    if feature.type == 'CDS':
                        gene_id = feature.qualifiers['locus_tag'][0]
                        element_file.write('ID\t' + gene_id + '\n')
                        element_file.write('NAME\t' + gene_id + '\n')
                        element_file.write('STARTBASE\t' + str(feature.location.start) + '\n')
                        element_file.write('ENDBASE\t' + str(feature.location.end) + '\n')
                        if 'function' in feature.qualifiers:
                            for function in feature.qualifiers['function']:
                                element_file.write('FUNCTION\t' + function + '\n')
                        if 'EC_number' in feature.qualifiers:
                            for ec in feature.qualifiers['EC_number']:
                                element_file.write('EC\t' + ec + '\n')
                        if 'go_component' in feature.qualifiers:
                            for go in feature.qualifiers['go_component']:
                                element_file.write('DBLINK\t' + go + '\n')
                        element_file.write('PRODUCT-TYPE\tP' + '\n')
                        element_file.write('PRODUCT-ID\tprot ' + gene_id + '\n')
                        element_file.write('//\n\n')
