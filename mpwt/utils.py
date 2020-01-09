"""
Uselful functions for mpwt.
"""

import csv
import gffutils
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


def create_pathologic_file(input_folder, output_folder, number_cpu=None):
    """
    Create PathoLogic file from Genbank or GFF files.

    Args:
        input_folder (str): pathname to the folder containing Genbanks or GFFs
        output_folder (str): pathname to the output folder containing the PathoLogic files
        number_cpu (str): number of CPU
    """
    if number_cpu:
        number_cpu_to_use = int(number_cpu)
    else:
        number_cpu_to_use = 1

    multiprocessing_input_data = []

    mpwt_pool = Pool(processes=number_cpu_to_use)

    input_names = os.listdir(input_folder)

    if 'taxon_id.tsv' in input_names:
        taxon_ids = {}
        input_names.remove('taxon_id.tsv')
        with open(input_folder + '/taxon_id.tsv') as taxon_file:
            for row in csv.reader(taxon_file, delimiter='\t'):
                taxon_ids[row[0]] = row[1]
    else:
        taxon_ids = None

    for input_name in input_names:
        input_path_gbk = input_folder + '/' + input_name + '/' + input_name + '.gbk'
        input_path_gff = input_folder + '/' + input_name + '/' + input_name + '.gff'
        if os.path.exists(input_path_gbk):
            input_path = input_path_gbk
        elif os.path.exists(input_path_gff):
            input_path = input_path_gff
        elif all([True for species_file in os.listdir(input_folder + '/' + input_name + '/') if '.pf' in species_file or '.fasta' in species_file]):
            input_path = input_folder + '/' + input_name + '/'
        else:
            sys.exit('No .gff or .gbk file in ' + input_folder + '/' + input_name)

        output_path = output_folder + '/' + input_name

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        multiprocessing_dict = {'input_path': input_path, 'output_path': output_path,
                                'output_folder': output_folder, 'input_name': input_name}
        if taxon_ids:
            if input_name in taxon_ids:
                multiprocessing_dict['taxon_id'] = taxon_ids[input_name]

        multiprocessing_input_data.append(multiprocessing_dict)

    mpwt_pool.map(run_create_pathologic_file, multiprocessing_input_data)

    mpwt_pool.close()
    mpwt_pool.join()


def write_taxon_id_file(input_name, taxon_id, output_folder):
    if not os.path.exists(output_folder + '/taxon_id.tsv'):
        with open(output_folder + '/taxon_id.tsv', 'w') as taxon_id_file:
            taxon_writer = csv.writer(taxon_id_file, delimiter='\t')
            taxon_writer.writerow(['species', 'taxon_id'])
            taxon_writer.writerow([input_name, taxon_id])
    else:
        with open(output_folder + '/taxon_id.tsv', 'a') as taxon_id_file:
            taxon_writer = csv.writer(taxon_id_file, delimiter='\t')
            taxon_writer.writerow([input_name, taxon_id])


def run_create_pathologic_file(multiprocessing_input_data):
    """
    Create PathoLogic files from a Genbank or a GFF file.

    Args:
        multiprocess_input (dictionary): contains multiprocess input (input folder, output_path, output folder and input_name)
    """
    input_path = multiprocessing_input_data['input_path']
    output_folder = multiprocessing_input_data['output_folder']
    output_path = multiprocessing_input_data['output_path']
    input_name = multiprocessing_input_data['input_name']
    taxon_id = None
    # Add taxon ID in taxon_id.tsv if available.
    if input_path.endswith('.gbk'):

        if not os.path.exists(output_path):
            os.makedirs(output_path)

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
                    logger.info('No taxon ID in the Genbank {0} In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.'.format(input_path))
            if taxon_id:
                write_taxon_id_file(input_name, taxon_id, output_folder)

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
                        gene_name = None
                        gene_id = None
                        if 'locus_tag' in feature.qualifiers:
                            gene_id = feature.qualifiers['locus_tag'][0]
                        if 'gene' in feature.qualifiers:
                            gene_name = feature.qualifiers['gene'][0]
                        if not gene_id and not gene_name:
                            logger.critical('No locus_tag and no gene qualifiers in feature of record: ' + record.id)
                            pass
                        if gene_id:
                            element_file.write('ID\t' + gene_id + '\n')
                        else:
                            if gene_name:
                                element_file.write('ID\t' + gene_name + '\n')
                        if gene_name:
                            element_file.write('NAME\t' + gene_name + '\n')
                        else:
                            if gene_id:
                                element_file.write('NAME\t' + gene_id + '\n')
                        element_file.write('STARTBASE\t' + str(feature.location.start+1) + '\n')
                        element_file.write('ENDBASE\t' + str(feature.location.end) + '\n')
                        if 'function' in feature.qualifiers:
                            for function in feature.qualifiers['function']:
                                element_file.write('FUNCTION\t' + function + '\n')
                        if 'EC_number' in feature.qualifiers:
                            for ec in feature.qualifiers['EC_number']:
                                element_file.write('EC\t' + ec + '\n')
                        if 'go_component' in feature.qualifiers:
                            for go in feature.qualifiers['go_component']:
                                element_file.write('GO\t' + go + '\n')
                        if 'go_function' in feature.qualifiers:
                            for go in feature.qualifiers['go_function']:
                                element_file.write('GO\t' + go + '\n')
                        if 'go_process' in feature.qualifiers:
                            for go in feature.qualifiers['go_process']:
                                element_file.write('GO\t' + go + '\n')
                        element_file.write('PRODUCT-TYPE\tP' + '\n')
                        if gene_id:
                            element_file.write('PRODUCT-ID\tprot ' + gene_id + '\n')
                        else:
                            if gene_name:
                                element_file.write('PRODUCT-ID\tprot ' + gene_name + '\n')
                        element_file.write('//\n\n')

    elif input_path.endswith('.gff'):

        if not os.path.exists(output_path):
            os.makedirs(output_path)

        gff_database = gffutils.create_db(input_path, ':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        regions = list(set([region.chrom for region in gff_database.features_of_type('region')]))
        try:
            first_region = list(set([region for region in gff_database.features_of_type('region')]))[0]
        except IndexError:
            raise IndexError("No 'region' feature in the GFF " + input_path)
        if 'Dbxref' in first_region.attributes:
            for dbxref in first_region.attributes['Dbxref']:
                if 'taxon' in dbxref:
                    taxon_id = dbxref.replace('taxon:', '')
        if taxon_id:
            write_taxon_id_file(input_name, taxon_id, output_folder)

        for record in SeqIO.parse(input_path.replace('.gff', '.fasta'), 'fasta'):
            output_fasta = output_path + '/' + record.id + '.fasta'
            SeqIO.write(record, output_fasta, 'fasta')
        for region in regions:
            with open(output_path + '/' + region + '.pf', 'w') as element_file:
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                element_file.write(';; ' + region + '\n')
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                for feature in gff_database.features_of_type(tuple(gff_database.featuretypes())):
                    if feature.featuretype == 'gene':
                        if feature.chrom == region:
                            element_file.write('ID\t' + feature.id + '\n')
                            element_file.write('NAME\t' + feature.id + '\n')
                            element_file.write('STARTBASE\t' + str(feature.start) + '\n')
                            element_file.write('ENDBASE\t' + str(feature.stop) + '\n')
                            element_file.write('PRODUCT-TYPE\tP' + '\n')
                            element_file.write('PRODUCT-ID\tprot ' + feature.id + '\n')
                            for child in gff_database.children(feature.id):
                                if 'product' in child.attributes:
                                    for product in child.attributes['product']:
                                        element_file.write('FUNCTION\t' + product + '\n')
                                if 'ec_number' in child.attributes:
                                    for ec in child.attributes['ec_number']:
                                        element_file.write('EC\t' + ec + '\n')
                            element_file.write('//\n\n')

    elif all([True for species_file in os.listdir(input_path) if '.pf' in species_file or '.fasta' in species_file]):
        taxon_id = multiprocessing_input_data['taxon_id']
        write_taxon_id_file(input_name, taxon_id, output_folder)
        shutil.copytree(input_path, output_path)


def pubmed_citations(activate_citations):
    """
    Activate or deactivate loading of PubMed citations.

    Args:
    activate_citations (bool): boolean to indicate if you want to activate or not the downlaod of Pubmed entries.
    """
    ptools_init_filepath = find_ptools_path() + '/ptools-init.dat'
    new_ptools_file = ""

    download_pubmed_entries_parameter = None
    with open(ptools_init_filepath, 'r') as ptools_init_file:
        for line in ptools_init_file.read().split('\n'):
            if 'Batch-PathoLogic-Download-Pubmed-Entries?' in line:
                if '#' in line:
                    line = line.replace('#', '')
                download_pubmed_entries_parameter = True
                if activate_citations:
                    line = line.replace('nil', 'T')
                else:
                    line = line.replace('T', 'nil')
            if line != '':
                new_ptools_file = new_ptools_file + line + '\n'
            else:
                new_ptools_file = new_ptools_file + line

    if not download_pubmed_entries_parameter:
        sys.exit('There is no Batch-PathoLogic-Download-Pubmed-Entries parameter in ' + ptools_init_filepath +'. To use --nc/no_download_articles, mpwt needs Pathway Tools 23.5 or higher.')

    with open(ptools_init_filepath, 'w') as ptools_init_file:
        ptools_init_file.write(new_ptools_file)
