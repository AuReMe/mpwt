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
Scripts to convert GenBank files and GFF files into PathoLogic Files.
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

    if not os.path.exists(input_folder):
        logger.critical('mpwt can not run: ' + input_folder + ' does not exist.')
        return
    if not os.path.isdir(input_folder):
        logger.critical('mpwt can not run: ' + input_folder + ' is not a directory.')
        return

    input_names = os.listdir(input_folder)

    # If there is a taxon_id file in input read it and extract informations.
    if 'taxon_id.tsv' in input_names:
        taxon_ids = {}
        input_names.remove('taxon_id.tsv')
        taxon_id_path = os.path.join(input_folder, 'taxon_id.tsv')
        with open(taxon_id_path, 'r') as taxon_file:
            for row in csv.reader(taxon_file, delimiter='\t'):
                taxon_ids[row[0]] = row[1]
    else:
        taxon_ids = None

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Initiate the output taxon_id.tsv file.
    output_taxon_id_path = os.path.join(output_folder, 'taxon_id.tsv')
    with open(output_taxon_id_path, 'w', encoding='utf-8') as taxon_id_file:
        taxon_writer = csv.writer(taxon_id_file, delimiter='\t')
        taxon_writer.writerow(['species', 'taxon_id'])

    # For each input, search for a GenBank, a GFF or PathoLogic files.
    for input_name in input_names:
        input_path_gbk = os.path.join(*[input_folder, input_name, input_name + '.gbk'])
        input_path_gbff = os.path.join(*[input_folder, input_name, input_name + '.gbff'])
        input_path_gff = os.path.join(*[input_folder, input_name, input_name + '.gff'])
        if os.path.exists(input_path_gbk):
            input_path = input_path_gbk
        elif os.path.exists(input_path_gbff):
            input_path = input_path_gbff
        elif os.path.exists(input_path_gff):
            input_path = input_path_gff
        elif any([species_file.endswith('.pf') for species_file in os.listdir(os.path.join(input_folder, input_name))]):
            input_path = os.path.join(input_folder, input_name)
        else:
            sys.exit('No .gff/.gbk/.gbff/.pf file in ' + os.path.join(input_folder, input_name))

        output_path = os.path.join(output_folder, input_name)

        if taxon_ids:
            if input_name in taxon_ids:
                taxon_id = taxon_ids[input_name]
            else:
                taxon_id = None
        else:
            taxon_id = None

        multiprocessing_input_data.append([input_path, output_path, output_folder, input_name, taxon_id])

    # Convert GenBank/GFF files to PathoLogic files.
    check_boolean = mpwt_pool.starmap(run_create_pathologic_file, multiprocessing_input_data)

    mpwt_pool.close()
    mpwt_pool.join()

    if not all(check_boolean):
        sys.exit('Issue during PathoLogic files creation, check errors.')


def write_taxon_id_file(input_name, taxon_id, output_folder):
    """
    Add species to taxon_id.tsv file.

    Args:
        input_name (str): name of the species
        taxon_id (str): taxon_id linked to the species
        output_folder (str): path to output folder
    """
    taxon_id_path = os.path.join(output_folder, 'taxon_id.tsv')

    with open(taxon_id_path, 'a') as taxon_id_file:
        taxon_writer = csv.writer(taxon_id_file, delimiter='\t')
        taxon_writer.writerow([input_name, taxon_id])


def run_create_pathologic_file(input_path, output_path, output_folder, input_name, taxon_id):
    """
    Create PathoLogic files from a Genbank or a GFF file.

    Args:
        input_path (str): path to species input folder
        output_path (str): path to output species folder
        output_folder (str): path to output folder
        input_name (str): species name
        taxon_id (dictionary): dictionary with the taxon_id for each species, if taxon_id.tsv does not exit None
    """
    # Add taxon ID in taxon_id.tsv if available.
    if input_path.endswith('.gbk') or input_path.endswith('.gbff'):
        logger.info('Creating PathoLogic file for ' + input_path)
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
            output_fasta_path = os.path.join(output_path, element_id + '.fasta')
            output_pf_path = os.path.join(output_path, element_id + '.pf')
            SeqIO.write(records, output_fasta_path, 'fasta')
            with open(output_pf_path, 'w', encoding='utf-8') as element_file:
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                element_file.write(';; ' + element_id + '\n')
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                for feature in record.features:
                    if feature.type in ['rRNA', 'tRNA', 'ncRNA', 'misc_RNA', 'CDS']:
                        gene_name = None
                        gene_id = None
                        start_location = str(feature.location.start+1).replace('<', '').replace('>', '')
                        end_location = str(feature.location.end).replace('<', '').replace('>', '')
                        if 'locus_tag' in feature.qualifiers:
                            gene_id = feature.qualifiers['locus_tag'][0]
                        if 'gene' in feature.qualifiers:
                            gene_name = feature.qualifiers['gene'][0]
                        if not gene_id and not gene_name:
                            logger.critical('No locus_tag and no gene qualifiers in feature of record: ' + record.id + ' at position ' + start_location + '-' +end_location)
                            pass
                        if gene_name:
                            if len(gene_name) > 40:
                                logger.critical('Critical warning: gene ID ' + gene_id + ' of ' + feature.type + ' of file ' + input_path + 'is too long (more than 40 characters), this will cause errors in Pathway Tools.')
                            element_file.write('NAME\t' + gene_name + '\n')
                        else:
                            if gene_id:
                                if len(gene_id) > 40:
                                    logger.critical('Critical warning: gene ID ' + gene_id + ' of ' + feature.type + ' of file ' + input_path + 'is too long (more than 40 characters), this will cause errors in Pathway Tools.')
                                element_file.write('NAME\t' + gene_id + '\n')
                        if gene_id and gene_id != gene_name:
                            element_file.write('ID\t' + gene_id + '\n')
                        element_file.write('STARTBASE\t' + start_location + '\n')
                        element_file.write('ENDBASE\t' + end_location + '\n')
                        if 'product' in feature.qualifiers:
                            for function in feature.qualifiers['product']:
                                element_file.write('FUNCTION\t' + function + '\n')
                        if 'db_xref' in feature.qualifiers:
                            for db_xref in feature.qualifiers['db_xref']:
                                if ':' in db_xref:
                                    element_file.write('DBLINK\t' + db_xref + '\n')
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
                        if 'gene_synonym' in feature.qualifiers:
                            for gene_synonym in feature.qualifiers['gene_synonym']:
                                element_file.write('SYNONYM\t' + gene_synonym + '\n')
                        if feature.type == 'rRNA':
                            if 'pseudo' in feature.qualifiers:
                                element_file.write('PRODUCT-TYPE\tPSEUDO' + '\n')
                            else:
                                element_file.write('PRODUCT-TYPE\tRRNA' + '\n')
                            if gene_id:
                                element_file.write('PRODUCT-ID\trnra ' + gene_id + '\n')
                            else:
                                if gene_name:
                                    element_file.write('PRODUCT-ID\trnra ' + gene_name + '\n')
                            element_file.write('//\n\n')
                        if feature.type == 'misc_RNA' or feature.type == 'ncRNA':
                            if 'pseudo' in feature.qualifiers:
                                element_file.write('PRODUCT-TYPE\tPSEUDO' + '\n')
                            else:
                                element_file.write('PRODUCT-TYPE\tMISC-RNA' + '\n')
                            if gene_id:
                                element_file.write('PRODUCT-ID\tmiscnra ' + gene_id + '\n')
                            else:
                                if gene_name:
                                    element_file.write('PRODUCT-ID\tmiscnra ' + gene_name + '\n')
                            element_file.write('//\n\n')
                        if feature.type == 'tRNA':
                            if 'pseudo' in feature.qualifiers:
                                element_file.write('PRODUCT-TYPE\tPSEUDO' + '\n')
                            else:
                                element_file.write('PRODUCT-TYPE\tTRNA' + '\n')
                            if gene_id:
                                element_file.write('PRODUCT-ID\ttnra ' + gene_id + '\n')
                            else:
                                if gene_name:
                                    element_file.write('PRODUCT-ID\ttnra ' + gene_name + '\n')
                            element_file.write('//\n\n')
                        if feature.type == 'CDS':
                            if 'pseudo' in feature.qualifiers:
                                element_file.write('PRODUCT-TYPE\tPSEUDO' + '\n')
                            else:
                                element_file.write('PRODUCT-TYPE\tP' + '\n')
                            if gene_id:
                                element_file.write('PRODUCT-ID\tprot ' + gene_id + '\n')
                            else:
                                if gene_name:
                                    element_file.write('PRODUCT-ID\tprot ' + gene_name + '\n')
                            element_file.write('//\n\n')

                    elif feature.type == 'mRNA':
                         if 'pseudo' in feature.qualifiers:
                            start_location = str(feature.location.start+1).replace('<', '').replace('>', '')
                            end_location = str(feature.location.end).replace('<', '').replace('>', '')
                            if 'locus_tag' in feature.qualifiers:
                                gene_id = feature.qualifiers['locus_tag'][0]
                            if 'gene' in feature.qualifiers:
                                gene_name = feature.qualifiers['gene'][0]
                            if not gene_id and not gene_name:
                                logger.critical('No locus_tag and no gene qualifiers in feature of record: ' + record.id + ' at position ' + start_location + '-' +end_location)
                                pass
                            if gene_id:
                                if len(gene_id) > 40:
                                    logger.critical('Critical warning: gene ID ' + gene_id + ' of ' + feature.type + ' of file ' + input_path + 'is too long (more than 40 characters), this will cause errors in Pathway Tools.')
                                element_file.write('ID\t' + gene_id + '\n')
                            else:
                                if gene_name:
                                    if len(gene_name) > 40:
                                        logger.critical('Critical warning: gene ID ' + gene_id + ' of ' + feature.type + ' of file ' + input_path + 'is too long (more than 40 characters), this will cause errors in Pathway Tools.')
                                    element_file.write('ID\t' + gene_name + '\n')
                            if gene_name:
                                element_file.write('NAME\t' + gene_name + '\n')
                            else:
                                if gene_id:
                                    element_file.write('NAME\t' + gene_id + '\n')
                            element_file.write('STARTBASE\t' + start_location + '\n')
                            element_file.write('ENDBASE\t' + end_location + '\n')
                            element_file.write('PRODUCT-TYPE\tPSEUDO' + '\n')
                            if 'product' in feature.qualifiers:
                                for function in feature.qualifiers['product']:
                                    element_file.write('FUNCTION\t' + function + '\n')
                            if 'db_xref' in feature.qualifiers:
                                for db_xref in feature.qualifiers['db_xref']:
                                    if ':' in db_xref:
                                        element_file.write('DBLINK\t' + db_xref + '\n')
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

        gff_fasta = None
        if os.path.exists(input_path.replace('.gff', '.fasta')):
            gff_fasta = input_path.replace('.gff', '.fasta')
        elif os.path.exists(input_path.replace('.gff', '.fsa')):
            gff_fasta = input_path.replace('.gff', '.fsa')
        if not gff_fasta:
            logger.critical('Error: No fasta file (.fasta or .fsa) for GFF file ' + input_path)
            return None

        for record in SeqIO.parse(gff_fasta, 'fasta'):
            gff_fasta_extension = os.path.splitext(gff_fasta)[1]
            output_fasta = os.path.join(output_path, record.id + gff_fasta_extension)
            SeqIO.write(record, output_fasta, 'fasta')

        for region in regions:
            region_pf_path = os.path.join(output_path, region + '.pf')
            with open(region_pf_path, 'w', encoding='utf-8') as element_file:
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                element_file.write(';; ' + region + '\n')
                element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
                for feature in gff_database.features_of_type(tuple(gff_database.featuretypes())):
                    if feature.featuretype == 'gene':
                        if feature.chrom == region:
                            if len(feature.id) > 40:
                                logger.critical('Critical warning: gene ID ' + feature.id + ' of ' + input_path + 'is too long (more than 40 characters), this will cause errors in Pathway Tools.')
                            element_file.write('ID\t' + feature.id + '\n')
                            element_file.write('NAME\t' + feature.id + '\n')
                            element_file.write('STARTBASE\t' + str(feature.start) + '\n')
                            element_file.write('ENDBASE\t' + str(feature.stop) + '\n')
                            for child in gff_database.children(feature.id):
                                if 'product' in child.attributes:
                                    for product in child.attributes['product']:
                                        element_file.write('FUNCTION\t' + product + '\n')
                                if 'ec_number' in child.attributes:
                                    for ec in child.attributes['ec_number']:
                                        element_file.write('EC\t' + ec + '\n')
                                if 'db_xref' in child.attributes:
                                    for dbxref in child.attributes['db_xref']:
                                        if ':' in dbxref:
                                            element_file.write('DBLINK\t' + dbxref + '\n')
                                if 'Dbxref' in child.attributes:
                                    for dbxref in child.attributes['Dbxref']:
                                        if ':' in dbxref:
                                            element_file.write('DBLINK\t' + dbxref + '\n')
                                if child.featuretype == 'CDS':
                                    element_file.write('PRODUCT-TYPE\tP' + '\n')
                                    element_file.write('PRODUCT-ID\tprot ' + feature.id + '\n')
                                elif child.featuretype == 'tRNA':
                                    element_file.write('PRODUCT-TYPE\tTRNA' + '\n')
                                    element_file.write('PRODUCT-ID\ttrna ' + feature.id + '\n')
                                elif child.featuretype == 'rRNA':
                                    element_file.write('PRODUCT-TYPE\tRRNA' + '\n')
                                    element_file.write('PRODUCT-ID\rrna ' + feature.id + '\n')
                                elif child.featuretype == 'pseudogene':
                                    element_file.write('PRODUCT-TYPE\tPSEUDO' + '\n')

                            element_file.write('//\n\n')

    elif any([species_file.endswith('.pf') for species_file in os.listdir(input_path)]):
        write_taxon_id_file(input_name, taxon_id, output_folder)
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        for species_file in os.listdir(input_path):
            species_file_path = os.path.join(input_path, species_file)
            if species_file.endswith('.pf'):
                if os.path.exists(species_file_path.replace('.pf', '.fasta')):
                    fasta_file = species_file.replace('.pf', '.fasta')
                elif os.path.exists(species_file_path.replace('.pf', '.fsa')):

                    fasta_file = species_file.replace('.pf', '.fsa')
                else:
                    logger.critical('Error: No fasta file (.fasta or .fsa) for PathoLogic file ' + species_file_path)
                    return None
                output_species_path = os.path.join(output_path, species_file)
                input_fasta_path = os.path.join(input_path, fasta_file)
                output_fasta_path = os.path.join(output_path, fasta_file)
                shutil.copy2(species_file_path, output_species_path)
                shutil.copy2(input_fasta_path, output_fasta_path)

    return True