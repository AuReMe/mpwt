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
Check input folders and input files.

Create PathoLogic input files:
-organism-params.dat
-genetic-elements.dat
-flat_files_creation.lisp
"""

import csv
import logging
import os

from Bio import SeqIO
from gffutils.iterators import DataIterator
from mpwt import utils

logger = logging.getLogger(__name__)


def compare_input_ids_to_ptools_ids(compare_ids, ptools_run_ids, set_operation):
    """ Compare species IDs in input folder with IDs present in PGDB folder.

    Args:
        compare_ids (list): species IDs (folder and GBK/GFF file name)
        ptools_run_ids (list): PGDB IDs inside ptools-local
        set_operation (str): difference or intersection between the two inputs list
    Returns:
        list: ID intersecting or different between the two lists
    """
    # To compare them, lower case the species IDs, then run a difference or an intersection between set.
    lower_case_compare_ids = list(map(lambda x:x.lower(), compare_ids))

    lower_compare_ids = dict(zip(lower_case_compare_ids, compare_ids))

    # Difference to obtain all IDs which are not in PGDB folder (which need to be run on PathoLogic).
    if set_operation == 'difference':
        compare_ids_ptools = set(lower_case_compare_ids) - set(ptools_run_ids)

    # Intersection to obtain all IDs which are already in PGDB folder.
    # Which need only to create BioPAX/flat files and move them in output folder.
    elif set_operation == 'intersection':
        compare_ids_ptools = set(ptools_run_ids).intersection(set(lower_case_compare_ids))

    new_compare_ids = [lower_compare_ids[compare_id] for compare_id in compare_ids_ptools]

    return new_compare_ids


def check_input_and_existing_pgdb(run_ids, input_folder, output_folder, number_cpu_to_use):
    """ Check input structure and data in output folder and ptools-local.

    Args:
        run_ids (list): species IDs (folder and GBK/GFF file name)
        input_folder (str): pathname to the input folder
        output_folder (str): pathname to the output folder
        number_cpu_to_use (int): number of CPU to use for multiprocessing
    Returns:
        list: input IDs for PathoLogic and BioPAX/flat files creation
        list: input IDs for BioPAX/flat files creation
    """
    # Check if there are files/folders inside the input folder.
    # And do not use hidden folder/file (beginning with '.').
    species_folders = [species_folder for species_folder in os.listdir(input_folder) if not species_folder.startswith('.')]
    if len(species_folders) == 0:
        logger.critical("No folder containing genbank/gff file. In {0} you must have sub-folders containing Genbank/GFF file.".format(input_folder))
        return None, None

    # Remove Pathologic taxon ID file.
    if 'taxon_id.tsv' in species_folders:
        species_folders.remove('taxon_id.tsv')

    # Check if there is a Genbank, a GFF or a PathoLogic file inside each subfolder.
    check_species_folders = []
    for species_folder in species_folders:
        species_input_files = []
        species_folder_path = os.path.join(input_folder, species_folder)
        for species_file in os.listdir(species_folder_path):
            species_filename, species_file_extension = os.path.splitext(species_file)
            if species_file_extension in ['.gbk', '.gbff', '.gff']:
                if species_filename == species_folder:
                    check_species_folders.append(species_folder)
                    species_input_files.append(species_file_extension)
            if any(input_extension in species_file for input_extension in ['.pf']):
                check_species_folders.append(species_folder)
                species_input_files.append(species_file_extension)
        species_input_files = list(set(species_input_files))
        if len(species_input_files) > 1:
            logger.critical('Multiple input files for {0}, there must be only one type of files among: GenBank, GFF or multiple PF files'.format(species_folder))
            return None, None
        elif len(species_input_files) == 0:
            logger.critical('Missing input file for {0}. A GenBank file, GFF file or multiple PF files are required.'.format(species_folder))
            return None, None

    check_species_folders = list(set(check_species_folders))

    missing_input_files = list(set(run_ids) - set(check_species_folders))
    if len(check_species_folders) == 0:
        logger.critical('Missing Genbank/GFF/PF file for: {0} \nCheck for input files (.gbk/.gbff/.gff/.pf)'.format(','.join(missing_input_files)))
        return None, None

    # Check the structure of the input folder.
    invalid_characters = ['.', '/']
    for species_folder in check_species_folders:
        species_folder_path = os.path.join(input_folder, species_folder)
        if os.path.isfile(species_folder_path):
            logger.critical('Error: file inside the input_folder ({0}) instead of a subfolder. Check that you have a structure file of input_folder/species_1/species1.gbk and not input_folder/species_1.gbk.'.format(species_folder_path))
            return None, None
        elif os.path.isdir(species_folder_path):
            if any(char in invalid_characters for char in species_folder):
                logger.critical('Error: . or / in genbank/gff name {0} \nGenbank name is used as an ID in Pathway Tools and Pathway Tools does not create PGDB with . in ID.'.format(species_folder))
                return None, None

    # Take run_ids and remove folder with error (with the intersection with check_species_folders) and if there is already present output.
    clean_run_ids = set(run_ids).intersection(set(check_species_folders))

    if output_folder:
        if os.path.exists(output_folder):
            if os.path.isdir(output_folder):
                already_present_outputs = [output_pgdb for output_pgdb in os.listdir(output_folder)]
                new_run_ids = clean_run_ids - set(already_present_outputs)
                new_run_ids = list(new_run_ids)
                for pgdb in already_present_outputs:
                    if pgdb in clean_run_ids:
                        logger.warning("! PGDB {0} already in output folder {1}, no inference will be launched on this species.".format(pgdb, output_folder))

                if len(new_run_ids) == 0:
                    logger.info("All PGDBs are already present in the output folder. Remove them if you want a new inference.")
                    return None, None
            else:
                logger.info(output_folder + " is not a valid output folder.")
                return None, None
        else:
            new_run_ids = list(clean_run_ids)

    else:
        new_run_ids = list(clean_run_ids)

    # Check for PGDB in ptools-local to see if PGDB are already present but they haven't been exported.
    already_present_pgdbs = [pgdb_species_folder[:-3] for pgdb_species_folder in utils.list_pgdb()]

    # Check the already finished PGDBs.
    if already_present_pgdbs:
        pathologic_builds = compare_input_ids_to_ptools_ids(new_run_ids, already_present_pgdbs, 'intersection')

        # Check for unfinished build of PGDB using their pathologic.log file.
        logger.info("Check and delete unfinished builds of Pathway Tools.")
        unfinished_builds = []
        finished_builds = []
        for pathologic_build in pathologic_builds:
            pathologic_build_lower = pathologic_build.lower()
            pathologic_file = os.path.join(*[input_folder, pathologic_build, 'pathologic.log'])
            if os.path.exists(pathologic_file):
                with open(pathologic_file, 'r') as pathologic_log:
                    pathologic_string = pathologic_log.read()
                    if 'Done' in pathologic_string:
                        finished_builds.append(pathologic_build_lower)
                    else:
                        unfinished_builds.append(pathologic_build_lower)

        # Delete the unfinished PGDBs.
        if unfinished_builds:
            utils.remove_pgdbs([unfinished_build + 'cyc' for unfinished_build in unfinished_builds], number_cpu_to_use)

        already_present_pgdbs = list(set(already_present_pgdbs) - set(unfinished_builds))

        run_patho_flat_ids = compare_input_ids_to_ptools_ids(new_run_ids, already_present_pgdbs, 'difference')
        run_flat_ids = compare_input_ids_to_ptools_ids(new_run_ids, already_present_pgdbs, 'intersection')

        for run_flat_id in run_flat_ids:
            logger.info("! PGDB {0} already in ptools-local, no PathoLogic inference will be launched on this species.".format(run_flat_id))
        return run_patho_flat_ids, run_flat_ids

    return new_run_ids, None


def create_flat_creation_script(pgdb_id, lisp_pathname):
    """ Create a lisp script allowing flat files creation.

    Args:
        pgdb_id (str): ID of a PGDB
        lisp_pathname (str): pathname to the output list script
    Returns:
        bool: True if lisp_pathname has been created
    """
    with open(lisp_pathname, 'w', encoding='utf-8') as lisp_file:
        lisp_file.write(';; Most files should begin by specifying that they will be interpreted')
        lisp_file.write('\n')
        lisp_file.write(';; within the Lisp package called ecocyc.')
        lisp_file.write('\n')
        lisp_file.write('(in-package :ecocyc)')
        lisp_file.write('\n')
        lisp_file.write(';; Select the organism PGDB as the current PGDB')
        lisp_file.write('\n')
        lisp_file.write("(select-organism :org-id '" + pgdb_id + ')')
        lisp_file.write('\n')
        lisp_file.write(';; Create attribute-values files without the progression pop-up')
        lisp_file.write('\n')
        lisp_file.write('(let ((*progress-noter-enabled?* NIL))')
        lisp_file.write('\n')
        lisp_file.write('        (create-flat-files-for-current-kb))')

    return os.path.isfile(lisp_pathname)


def extract_taxon_id(run_folder, pgdb_id, taxon_id, taxon_file):
    """ Extract taxon ID and other informations from taxon_id.tsv file.
    Other informatiosn are:
        - circular (cirularity of genome)
        - element_type (chromosome, plasmid, contig, mitochondria or chloroplast)
        - codon_table

    Args:
        run_folder (str): ID of a species of the input folder
        pgdb_id (str): ID of a PGDB
        taxon_id (str): Taxon ID for the corresponding species
        taxon_file (bool): Boolean indicating if a taxon_file must be used
    Returns:
        taxon_error (bool): Error status (True: error, False: no error).
        taxon_id (str): Taxon ID for the corresponding species
        taxon_datas (dict): Name of element file (or 'one_input' if only one file)
    """
    input_folder =  os.path.abspath(os.path.join(run_folder ,os.pardir))
    taxon_datas = {}
    taxon_id_found = None

    known_element_types = [':CHRSM', ':PLASMID', ':MT', ':PT', ':CONTIG']
    known_codon_table = ['0', '1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23']
    known_species = []

    taxon_id_path = os.path.join(input_folder, 'taxon_id.tsv')
    if not os.path.exists(taxon_id_path):
        logger.critical('Missing taxon_id.tsv file in {0}.'.format(input_folder))
        return True, None, None

    with open(taxon_id_path) as taxon_id_file:
        taxon_id_reader = csv.DictReader(taxon_id_file, delimiter='\t')
        for data in taxon_id_reader:
            if 'species' not in data:
                logger.critical('Missing "species" header in taxon_id.tsv file {0}.'.format(taxon_id_path))
                return True, None, None
            species = data['species']
            known_species.append(species)
            if pgdb_id == species:
                if 'taxon_id' in data:
                    if data['taxon_id'] != '':
                        if not taxon_id_found:
                            taxon_id = data['taxon_id']
                            taxon_id_found = True
                            logger.info('taxon_id.tsv: find taxon ID {0} for {1}'.format(taxon_id, pgdb_id))
                    else:
                        logger.critical('Missing taxon ID for {0} in {1}.'.format(pgdb_id, taxon_id_path))
                        return True, None, None
                else:
                    logger.critical('Missing taxon ID for {0} in {1}.'.format(pgdb_id, taxon_id_path))
                    return True, None, None

                if 'circular' in data:
                    if data['circular'] != '':
                        if data['circular'] == 'Y' or data['circular'] == 'N':
                            circular = data['circular']
                        else:
                            logger.critical('taxon_id.tsv: wrong circular for {0}, {1} instead of Y or N'.format(pgdb_id, data['circular']))
                            return True, None, None
                    else:
                        circular = None
                else:
                    circular = None

                if 'element_type' in data:
                    if data['element_type'] != '':
                        if data['element_type'] in known_element_types:
                            element_type = data['element_type']
                        else:
                            logger.critical('taxon_id.tsv: wrong element_type for {0}, {1} instead of {2}'.format(pgdb_id, data['element_type'], ', '.join(known_element_types)))
                            return True, None, None
                    else:
                        element_type = None
                else:
                    element_type = None

                if 'codon_table' in data:
                    if data['codon_table'] != '':
                        if data['codon_table'] in known_codon_table:
                            codon_table = data['codon_table']
                        else:
                            logger.critical('taxon_id.tsv: wrong codon_table for {0}, {1} instead of {2}'.format(pgdb_id, data['codon_table'], ', '.join(known_codon_table)))
                            return True, None, None
                    else:
                        codon_table = None
                else:
                    codon_table = None

                if 'corresponding_file' in data:
                    if data['corresponding_file'] != '':
                        corresponding_file = data['corresponding_file']
                        taxon_datas[corresponding_file] = {}
                        if circular is not None:
                            taxon_datas[corresponding_file]['circular'] = circular
                        if element_type is not None:
                            taxon_datas[corresponding_file]['element_type'] = element_type
                        if codon_table is not None:
                            taxon_datas[corresponding_file]['codon_table'] = codon_table
                    else:
                        if circular is not None:
                            taxon_datas['circular'] = circular
                        if element_type is not None:
                            taxon_datas['element_type'] = element_type
                        if codon_table is not None:
                            taxon_datas['codon_table'] = codon_table
                else:
                    if circular is not None:
                        taxon_datas['circular'] = circular
                    if element_type is not None:
                        taxon_datas['element_type'] = element_type
                    if codon_table is not None:
                        taxon_datas['codon_table'] = codon_table

                if 'reference_pgdb' in data:
                    if data['reference_pgdb'] != '':
                        taxon_datas['reference_pgdbs'] = data['reference_pgdb'].split(',')

    if pgdb_id not in known_species and taxon_id == '':
        logger.critical('Missing pgdb ID for {0} in {1}.'.format(pgdb_id, taxon_id_path))
        return True, None, None

    if taxon_file and pgdb_id not in known_species:
        logger.critical('Missing pgdb ID for {0} in {1}.'.format(pgdb_id, taxon_id_path))
        return True, None, None

    return False, taxon_id, taxon_datas


def create_flats_and_lisp(run_folder, taxon_file):
    """
    Read Genbank/GFF/PF files and create Pathway Tools needed file.
    Create also a lisp file to create flat files from Pathway tools results.
    The name of the PGDB created by Pathway Tools will be the name of the species with '_' instead of space.

    Create organism-params.dat:
    ID  pgdb_id
    STORAGE FILE
    NCBI-TAXON-ID   taxon_id
    NAME    species_name

    Create genetic-elements.dats:
    NAME    
    ANNOT-FILE  gbk_name
    //

    Create flat_files_creation.lisp:
    (in-package :ecocyc)
    (select-organism :org-id 'pgdb_id)
    (create-flat-files-for-current-kb)

    Args:
        run_folder (str): ID of a species of the input folder
        taxon_file (bool): Boolean indicating if a taxon_file must be used
    Returns:
        list: boolean list, True if all files have been created
   """
    # Look for a Genbank/GFF files in the run folder.
    # PGDB ID corresponds to the name of the species folder.
    pgdb_id = os.path.basename(run_folder)
    gbk_name = pgdb_id + ".gbk"
    gbk_pathname = os.path.join(run_folder, gbk_name)
    gbff_name = pgdb_id + ".gbff"
    gbff_pathname = os.path.join(run_folder, gbff_name)
    gff_name = pgdb_id + ".gff"
    gff_pathname = os.path.join(run_folder, gff_name)

    organism_dat = os.path.join(run_folder, 'organism-params.dat')
    genetic_dat = os.path.join(run_folder, 'genetic-elements.dat')
    lisp_pathname = os.path.join(run_folder, 'flat_files_creation.lisp')

    fasta_extensions = ['.fasta', '.fsa']

    taxon_id = ""
    taxon_error = False
    species_name = ""
    taxon_datas = {}

    if os.path.isfile(gbk_pathname) or os.path.isfile(gbff_pathname):
        if os.path.isfile(gbk_pathname):
            input_name = gbk_name
            input_path = gbk_pathname
        else:
            input_name = gbff_name
            input_path = gbff_pathname
        # Take the species name and the taxon id from the genbank file.
        with open(input_path, "r") as gbk:
            # Take the first record of the genbank (first contig/chromosome) to retrieve the species name.
            try:
                first_seq_record = next(SeqIO.parse(gbk, "genbank"))
            except StopIteration:
                logger.critical('Issue with the genbank {0}, it can be empty or malformatted.'.format(input_path))
                return None

            try:
                species_name = first_seq_record.annotations['organism']
            except KeyError:
                logger.critical('No organism in the Genbank {0} In the SOURCE you must have: ORGANISM  Species name'.format(pgdb_id))
                return None

            # Take the source feature of the first record.
            # This feature contains the taxon ID in the db_xref qualifier.
            src_features = [feature for feature in first_seq_record.features if feature.type == "source"]
            for src_feature in src_features:
                if 'db_xref' in src_feature.qualifiers:
                    src_dbxref_qualifiers = src_feature.qualifiers['db_xref']
                    for src_dbxref_qualifier in src_dbxref_qualifiers:
                        if 'taxon:' in src_dbxref_qualifier:
                            taxon_id = src_dbxref_qualifier.replace('taxon:', '')
                if not taxon_id:
                    logger.info('No taxon ID in the Genbank {0} In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.'.format(gbk_pathname))
                    logger.info('Try to look in the taxon_id.tsv file')
                    taxon_error, taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id, taxon_file)
            if taxon_file:
                taxon_error, taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id, taxon_file)

    elif os.path.isfile(gff_pathname):
        input_name = gff_name
        # Check if there is a fasta file.
        gff_fasta = None
        for fasta_extension in fasta_extensions:
            fasta_input_name = input_name.replace('.gff', fasta_extension)
            fasta_path = os.path.join(run_folder, fasta_input_name)
            if os.path.exists(fasta_path):
                gff_fasta = fasta_input_name
        if not gff_fasta:
            logger.critical('No fasta file (.fasta or .fsa) with the GFF of {0}'.format(pgdb_id))
            return None

        # Instead of parsing and creating a database from the GFF, parse the file and extract the first region feature.
        try:
            region_feature = [feature for feature in DataIterator(gff_pathname) if feature.featuretype == 'region'][0]
        except IndexError:
            logger.critical('No region feature in the GFF file of {0}, GFF file must have region features.'.format(pgdb_id))
            return None

        try:
            region_feature.attributes['Dbxref']
        except KeyError:
            logger.critical('No Dbxref in GFF file of {0} GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))

        for dbxref in region_feature.attributes['Dbxref']:
            if 'taxon' in dbxref:
                taxon_id = dbxref.split('taxon:')[1]
        if not taxon_id or taxon_file:
            if not taxon_id:
                logger.info('Missing "taxon:" in GFF file of {0} GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))
                logger.info('Try to look in the taxon_id.tsv file')
            taxon_error, taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id, taxon_file)

    # Look for PF files.
    elif all([True for species_file in os.listdir(run_folder) if '.pf' in species_file or '.fasta' in species_file or '.fsa' in species_file]):
        for species_file in os.listdir(run_folder):
            if '.pf' in species_file:
                # Check if there is a fasta file.
                pf_fasta = None
                for fasta_extension in fasta_extensions:
                    fasta_species_name = species_file.replace('.pf', fasta_extension)
                    fasta_path = os.path.join(run_folder, fasta_species_name)
                    if os.path.exists(fasta_path):
                        pf_fasta = fasta_species_name
                if not pf_fasta:
                    logger.critical('No fasta file (.fasta or .fsa) with the Pathologic file of {0}, this could lead to warnings in Pathway Tools.'.format(pgdb_id))

        taxon_error, taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id, taxon_file)

    if taxon_error == True:
        logger.critical('Issue with taxon ID of {0}.'.format(run_folder))
        return None

    # Create the organism-params dat file.
    with open(organism_dat, 'w', encoding='utf-8') as organism_file:
        organism_writer = csv.writer(organism_file, delimiter='\t', lineterminator='\n')
        organism_writer.writerow(['ID', pgdb_id])
        organism_writer.writerow(['STORAGE', "FILE"])
        organism_writer.writerow(['NCBI-TAXON-ID', taxon_id])
        organism_writer.writerow(['NAME', species_name])
        if 'reference_pgdbs' in taxon_datas:
            for reference_pgdb in taxon_datas['reference_pgdbs']:
                organism_writer.writerow(['REF-ORGID', reference_pgdb])

    # Create the genetic-elements dat file.
    with open(genetic_dat, 'w', encoding='utf-8') as genetic_file:
        if os.path.isfile(gff_pathname) or os.path.isfile(gbk_pathname) or os.path.isfile(gbff_pathname):
            genetic_writer = csv.writer(genetic_file, delimiter='\t', lineterminator='\n')
            genetic_writer.writerow(['NAME', ''])
            genetic_writer.writerow(['ANNOT-FILE', input_name])
            if os.path.isfile(gff_pathname):
                genetic_writer.writerow(['SEQ-FILE', gff_fasta])
            if 'circular' in taxon_datas:
                circular = taxon_datas['circular']
                genetic_writer.writerow(['CIRCULAR?', circular])
            if 'element_type' in taxon_datas:
                element_type = taxon_datas['element_type']
                genetic_writer.writerow(['TYPE', element_type])
            if 'codon_table' in taxon_datas:
                codon_table = taxon_datas['codon_table']
                genetic_writer.writerow(['CODON-TABLE', codon_table])
            genetic_writer.writerow(['//'])
        elif all([True for species_file in os.listdir(run_folder) if '.pf' in species_file or '.fasta' in species_file or '.fsa' in species_file]):
            genetic_writer = csv.writer(genetic_file, delimiter='\t', lineterminator='\n')
            for species_file in os.listdir(run_folder):
                    if '.pf' in species_file:
                        species_file_name = os.path.splitext(species_file)[0]
                        genetic_writer.writerow(['NAME', species_file.replace('.pf', '')])
                        genetic_writer.writerow(['ID', species_file.replace('.pf', '')])
                        genetic_writer.writerow(['ANNOT-FILE', species_file])
                        fasta_path = os.path.join(run_folder, species_file.replace('.pf', '.fasta'))
                        fsa_path = os.path.join(run_folder, species_file.replace('.pf', '.fsa'))
                        if os.path.exists(fasta_path):
                            genetic_writer.writerow(['SEQ-FILE', species_file.replace('.pf', '.fasta')])
                        elif os.path.exists(fsa_path):
                            genetic_writer.writerow(['SEQ-FILE', species_file.replace('.pf', '.fsa')])

                        if species_file_name in taxon_datas:
                            if 'circular' in taxon_datas[species_file_name]:
                                circular = taxon_datas[species_file_name]['circular']
                                genetic_writer.writerow(['CIRCULAR?', circular])
                            if 'element_type' in taxon_datas[species_file_name]:
                                element_type = taxon_datas[species_file_name]['element_type']
                                genetic_writer.writerow(['TYPE', element_type])
                            if 'codon_table' in taxon_datas[species_file_name]:
                                codon_table = taxon_datas[species_file_name]['codon_table']
                                genetic_writer.writerow(['CODON-TABLE', codon_table])
                        else:
                            if 'circular' in taxon_datas:
                                circular = taxon_datas['circular']
                                genetic_writer.writerow(['CIRCULAR?', circular])
                            if 'element_type' in taxon_datas:
                                element_type = taxon_datas['element_type']
                                genetic_writer.writerow(['TYPE', element_type])
                            if 'codon_table' in taxon_datas:
                                codon_table = taxon_datas['codon_table']
                                genetic_writer.writerow(['CODON-TABLE', codon_table])
                        genetic_writer.writerow(['//'])
    # Create the lisp script.
    check_lisp_file = create_flat_creation_script(pgdb_id, lisp_pathname)

    return all([os.path.isfile(organism_dat), os.path.isfile(genetic_dat), check_lisp_file])


def read_taxon_id(run_folder):
    """
    Search for Taxon ID in genbank or GFF files.
    For GenBank file searc for ''taxon:' key in 'db_xref' qualifier.
    For GFF file search for 'taxon' in dbxref feature.

    Args:
        run_folder (str): path to the input folder
    """
    taxon_ids = {}

    for input_folder in os.listdir(run_folder):
        input_folder_path = os.path.join(run_folder, input_folder)
        for input_file in os.listdir(input_folder_path):
            if '.gbk' in input_file:
                gbk_pathname = os.path.join(input_folder_path, input_file)
                # Take the species name and the taxon id from the genbank file.
                with open(gbk_pathname, "r") as gbk:
                    # Take the first record of the genbank (first contig/chromosome) to retrieve the species name.
                    first_seq_record = next(SeqIO.parse(gbk, "genbank"))
                    # Take the source feature of the first record.
                    # This feature contains the taxon ID in the db_xref qualifier.
                    src_features = [feature for feature in first_seq_record.features if feature.type == "source"]
                    for src_feature in src_features:
                        try:
                            src_dbxref_qualifiers = src_feature.qualifiers['db_xref']
                            for src_dbxref_qualifier in src_dbxref_qualifiers:
                                if 'taxon:' in src_dbxref_qualifier:
                                    taxon_id = src_dbxref_qualifier.replace('taxon:', '')
                        except KeyError:
                            logger.info('No taxon ID in the Genbank {0} In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.'.format(gbk_pathname))

            elif '.gff' in input_file:
                gff_pathname = os.path.join(input_folder_path, input_file)

                # Instead of parsing and creating a database from the GFF, parse the file and extract the first region feature.
                try:
                    region_feature = [feature for feature in DataIterator(gff_pathname) if feature.featuretype == 'region'][0]
                except IndexError:
                    raise IndexError('No region feature in the GFF file of {0}, GFF file must have region features.'.format(input_folder))

                try:
                    region_feature.attributes['Dbxref']
                except KeyError:
                    raise KeyError('No Dbxref in GFF file of {0} GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(input_folder))

                for dbxref in region_feature.attributes['Dbxref']:
                    if 'taxon' in dbxref:
                        taxon_id = dbxref.split('taxon:')[1]

            elif '.pf' in input_file:
                logger.info('No taxon ID associated to a PathoLogic Format. {0} will have a missing taxon_id'.format(input_folder))
                taxon_id = "missing"
        taxon_ids[input_folder] = taxon_id

    return taxon_ids


def pwt_input_files(run_folder, taxon_file):
    """
    Check if files needed by Pathway Tools are available, if not create them.
    Check if there is a pathologic.log from a previous run. If yes, delete it.

    Args:
        run_folder (str): path to the input folder
        taxon_file (str): path to the taxon_id.tsv file
    """
    required_files = set(['organism-params.dat', 'genetic-elements.dat', 'flat_files_creation.lisp'])
    files_in = set(next(os.walk(run_folder))[2])

    species_folder = os.path.basename(run_folder)
    pathologic_log = os.path.join(run_folder, 'pathologic.log')

    if 'pathologic.log' in files_in:
        os.remove(pathologic_log)

    error_found = False
    missing_string = ''
    if required_files.issubset(files_in):
        missing_string = 'No missing files'
    else:
        missing_string = 'Missing {0}'.format('; '.join(required_files.difference(files_in))) + '. Inputs file created for {0}'.format(species_folder)
        check_datas_lisp = create_flats_and_lisp(run_folder, taxon_file)
        if check_datas_lisp is None:
            logger.critical('|Input Check|{0}| Error with the creation of input files of {1}.'.format(species_folder, run_folder))
            error_found = True
            return error_found

    logger.info('|Input Check|{0}| {1}'.format(species_folder, missing_string))

    return error_found


def create_only_flat_lisp(pgdbs_folder_path, tmp_folder):
    """
    Create a lisp script file for each PGDB in the ptools-local folder.
    Return a generator with the PGDB IDs.

    Args:
        pgdbs_folder_path (str): pathname to species PGDB in ptools-local
        tmp_folder (str): termporary folder where lisp script will be stored
    Returns:
        generator: generator with the PGDB IDs.
    """
    for species_pgdb in os.listdir(pgdbs_folder_path):
        species_dir_path = os.path.join(pgdbs_folder_path, species_pgdb)
        if os.path.isdir(species_dir_path):
            pgdb_id = species_pgdb[:-3]
            pgdb_pathname = os.path.join(tmp_folder, pgdb_id)
            tmp_pgdb_path = os.path.join(tmp_folder, pgdb_id)
            if not os.path.exists(tmp_pgdb_path):
                os.mkdir(tmp_pgdb_path)
            lisp_pathname = os.path.join(pgdb_pathname, 'flat_files_creation.lisp')
            check_lisp_file = create_flat_creation_script(pgdb_id, lisp_pathname)
            if not check_lisp_file:
                raise Exception('Error with the creation of the lisp script for {0}'.format(species_pgdb))

            yield pgdb_id
