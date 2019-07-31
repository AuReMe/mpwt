"""
Create PathoLogic input files:
-organism-params.dat
-genetic-elements.dats
-dat_creation.lisp
"""

import csv
import logging
import os
import sys

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
    # Which need only to create BioPAX/dat files and move them in output folder.
    elif set_operation == 'intersection':
        compare_ids_ptools = set(ptools_run_ids).intersection(set(lower_case_compare_ids))

    new_compare_ids = [lower_compare_ids[compare_id] for compare_id in compare_ids_ptools]

    return new_compare_ids


def check_input_and_existing_pgdb(run_ids, input_folder, output_folder):
    """ Check input structure and data in output folder and ptools-local.

    Args:
        run_ids (list): species IDs (folder and GBK/GFF file name)
        input_folder (str): pathname to the input folder
        output_folder (str): pathname to the output folder
    Returns:
        list: input IDs for PathoLogic and BioPAX/dat creation
        list: input IDs for BioPAX/dat creation
    """
    # Check if there are files/folders inside the input folder.
    # And remove hidden folder/file (beginning with '.').
    species_folders = [species_folder for species_folder in os.listdir(input_folder) if not species_folder.startswith('.')]
    if len(species_folders) == 0:
        logger.critical("No folder containing genbank/gff file. In {0} you must have sub-folders containing Genbank/GFF file.".format(input_folder))
        return None, None

    # Remove Pathologic taxon ID file.
    if 'taxon_id.tsv' in species_folders:
        species_folders.remove('taxon_id.tsv')

    # Check if there is a Genbank, a GFF or a PathoLogic file inside each subfolder.
    input_extensions = ['.gbk', '.gff', '.pf']
    species_folders = list(set([species_folder for species_folder in species_folders
                                    for species_file in os.listdir(input_folder+'/'+species_folder)
                                        if any(input_extension in species_file for input_extension in input_extensions)]))
    missing_input_files = list(set(run_ids) - set(species_folders))
    if len(species_folders) == 0:
        logger.critical('Missing Genbank/GFF file for: {0} \nCheck if you have a Genbank file and if it ends with .gbk or .gff'.format(' '.join(missing_input_files)))
        return None, None

    # Check the structure of the input folder.
    invalid_characters = ['.', '/']
    for species_folder in species_folders:
        if os.path.isfile(input_folder+'/'+species_folder):
            logger.critical('Error: file inside the input_folder ({0}) instead of a subfolder. Check that you have a structure file of input_folder/species_1/species1.gbk and not input_folder/species_1.gbk.'.format(input_folder+'/'+species_folder))
            return None, None
        elif os.path.isdir(input_folder+'/'+species_folder):
            if any(char in invalid_characters for char in species_folder):
                logger.critical('Error: . or / in genbank/gff name {0} \nGenbank name is used as an ID in Pathway Tools and Pathway Tools does not create PGDB with . in ID.'.format(species_folder))
                return None, None

    # Take run_ids and remove folder with error (with the intersection with species_folders) and if there is already present output.
    clean_run_ids = set(run_ids).intersection(set(species_folders))
    if output_folder:
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
        new_run_ids = []
        for species_folder in species_folders:
            new_run_ids.append(species_folder)

    # Check for PGDB in ptools-local to see if PGDB are already present but they haven't been exported.
    already_present_pgdbs = [pgdb_species_folder[:-3] for pgdb_species_folder in utils.list_pgdb()]
    if already_present_pgdbs != []:
        run_patho_dat_ids = compare_input_ids_to_ptools_ids(new_run_ids, already_present_pgdbs, 'difference')
        run_dat_ids = compare_input_ids_to_ptools_ids(new_run_ids, already_present_pgdbs, 'intersection')
        for run_dat_id in run_dat_ids:
            logger.info("! PGDB {0} already in ptools-local, no PathoLogic inference will be launched on this species.".format(run_dat_id))
        return run_patho_dat_ids, run_dat_ids

    return new_run_ids, None


def create_dat_creation_script(pgdb_id, lisp_pathname):
    """ Create a lisp script allowing dat extraction.

    Args:
        pgdb_id (str): ID of a PGDB
        lisp_pathname (str): pathname to the output list script
    Returns:
        bool: True if lisp_pathname has been created
    """
    with open(lisp_pathname, 'w') as lisp_file:
        lisp_file.write('(in-package :ecocyc)')
        lisp_file.write('\n')
        lisp_file.write("(select-organism :org-id '" + pgdb_id + ')')
        lisp_file.write('\n')
        lisp_file.write('(let ((*progress-noter-enabled?* NIL))')
        lisp_file.write('\n')
        lisp_file.write('        (create-flat-files-for-current-kb))')

    return os.path.isfile(lisp_pathname)


def extract_taxon_id(run_folder, pgdb_id, taxon_id):
    """ Extract taxon ID from taxon_id.tsv file.

    Args:
        run_folder (str): ID of a species of the input folder
        pgdb_id (str): ID of a PGDB
    Returns:
        taxon_id (str): Taxon ID for the corresponding species
        taxon_datas (dict): Name of element file (or 'one_input' if only one file)
    """
    input_folder =  os.path.abspath(os.path.join(run_folder ,os.pardir))
    taxon_datas = {}
    taxon_id_found = None

    known_element_types = [':CHRSM', ':PLASMID', ':MT', ':PT', ':CONTIG']
    known_codon_table = ['0', '1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23']
    with open(input_folder + '/taxon_id.tsv') as pf_taxon_id:
        taxon_id_reader = csv.DictReader(pf_taxon_id, delimiter='\t')
        for data in taxon_id_reader:
            species = data['species']
            if pgdb_id == species:
                if 'taxon_id' in data:
                    if data['taxon_id'] != '':
                        if not taxon_id_found:
                            taxon_id = data['taxon_id']
                            taxon_id_found = True
                            logger.info('taxon_id.tsv: find taxon ID {0} for {1}'.format(taxon_id, pgdb_id))
                    else:
                        raise Exception('Missing taxon ID for {0} in {1}.'.format(pgdb_id, input_folder + '/taxon_id.tsv'))
                else:
                    raise Exception('Missing taxon ID for {0} in {1}.'.format(pgdb_id, input_folder + '/taxon_id.tsv'))

                if 'circular' in data:
                    if data['circular'] != '':
                        if data['circular'] == 'Y' or data['circular'] == 'N':
                            circular = data['circular']
                        else:
                            raise Exception('taxon_id.tsv: wrong circular for {0}, {1} instead of Y or N'.format(pgdb_id, circular))
                    else:
                        circular = None
                else:
                    circular = None

                if 'element_type' in data:
                    if data['element_type'] != '':
                        if data['element_type'] in known_element_types:
                            element_type = data['element_type']
                        else:
                            raise Exception('taxon_id.tsv: wrong circular for {0}, {1} instead of {2}'.format(pgdb_id, data['element_type'], ', '.join(known_element_types)))
                    else:
                        element_type = None
                else:
                    element_type = None

                if 'codon_table' in data:
                    if data['codon_table'] != '':
                        if data['codon_table'] in known_codon_table:
                            codon_table = data['codon_table']
                        else:
                            raise Exception('taxon_id.tsv: wrong circular for {0}, {1} instead of {2}'.format(pgdb_id, data['codon_table'], ', '.join(known_codon_table)))
                    else:
                        codon_table = None
                else:
                    codon_table = None

                if 'corresponding_file' in data:
                    if data['corresponding_file'] != '':
                        corresponding_file = data['corresponding_file']
                        taxon_datas[corresponding_file] = [circular, element_type, codon_table]
                    else:
                        taxon_datas['circular'] = circular
                        taxon_datas['element_type'] = element_type
                        taxon_datas['codon_table'] = codon_table
                else:
                        taxon_datas['circular'] = circular
                        taxon_datas['element_type'] = element_type
                        taxon_datas['codon_table'] = codon_table

    return taxon_id, taxon_datas


def create_dats_and_lisp(run_folder, taxon_file):
    """
    Read Genbank/GFF/PF files and create Pathway Tools needed file.
    Create also a lisp file to create dat files from Pathway tools results.
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

    Create dat_creation.lisp:
    (in-package :ecocyc)
    (select-organism :org-id 'pgdb_id)
    (create-flat-files-for-current-kb)

    Args:
        run_folder (str): ID of a species of the input folder
    Returns:
        list: boolean list, True if all files have been created
   """
    # Look for a Genbank/GFF files in the run folder.
    # PGDB ID corresponds to the name of the species folder.
    pgdb_id = run_folder.split('/')[-2]
    gbk_name = pgdb_id + ".gbk"
    gbk_pathname = run_folder + gbk_name
    gff_name = pgdb_id + ".gff"
    gff_pathname = run_folder + gff_name

    organism_dat = run_folder + 'organism-params.dat'
    genetic_dat = run_folder + 'genetic-elements.dat'

    taxon_id = ""
    species_name = ""
    taxon_datas = {}

    if os.path.isfile(gbk_pathname):
        input_name = gbk_name
        # Take the species name and the taxon id from the genbank file.
        with open(gbk_pathname, "r") as gbk:
            # Take the first record of the genbank (first contig/chromosome) to retrieve the species name.
            first_seq_record = next(SeqIO.parse(gbk, "genbank"))
            try:
                species_name = first_seq_record.annotations['organism']
            except KeyError:
                raise KeyError('No organism in the Genbank {0} In the SOURCE you must have: ORGANISM  Species name'.format(pgdb_id))

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
                    logger.info('Try to look in the taxon_id.tsv file')
                    taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id)
            if taxon_file:
                taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id)

    elif os.path.isfile(gff_pathname):
        input_name = gff_name
        # Check if there is a fasta file.
        try:
            with open(run_folder + input_name.replace('.gff', '.fasta'), 'r'):
                gff_fasta = input_name.replace('.gff', '.fasta')
        except FileNotFoundError:
            raise FileNotFoundError('No fasta file with the GFF of {0}'.format(pgdb_id))

        # Instead of parsing and creating a database from the GFF, parse the file and extract the first region feature.
        region_feature = [feature for feature in DataIterator(gff_pathname) if feature.featuretype == 'region'][0]
        try:
            region_feature.attributes['Dbxref']
        except KeyError:
            raise KeyError('No Dbxref in GFF file of {0} GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))

        for dbxref in region_feature.attributes['Dbxref']:
            if 'taxon' in dbxref:
                taxon_id = dbxref.split('taxon:')[1]
        if not taxon_id or taxon_file:
            if not taxon_id:
                logger.info('Missing "taxon:" in GFF file of {0} GFF file must have a ;Dbxref=taxon:taxonid; in the region feature.'.format(pgdb_id))
                logger.info('Try to look in the taxon_id.tsv file')
            taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id)

    # Look for PF files.
    elif all([True for species_file in os.listdir(run_folder) if '.pf' in species_file or '.fasta' in species_file]):
        for species_file in os.listdir(run_folder):
            if '.pf' in species_file:
                # Check if there is a fasta file.
                try:
                    pf_fasta = open(run_folder + species_file.replace('.pf', '.fasta'), 'r')
                    pf_fasta.close()
                except FileNotFoundError:
                    raise FileNotFoundError('No fasta file with the Pathologic file of {0}'.format(pgdb_id))

        taxon_id, taxon_datas = extract_taxon_id(run_folder, pgdb_id, taxon_id)

    lisp_pathname = run_folder + "dat_creation.lisp"

    # Create the organism-params dat file.
    with open(organism_dat, 'w') as organism_file:
        organism_writer = csv.writer(organism_file, delimiter='\t', lineterminator='\n')
        organism_writer.writerow(['ID', pgdb_id])
        organism_writer.writerow(['STORAGE', "FILE"])
        organism_writer.writerow(['NCBI-TAXON-ID', taxon_id])
        organism_writer.writerow(['NAME', species_name])

    # Create the genetic-elements dat file.
    with open(genetic_dat, 'w') as genetic_file:
        if os.path.isfile(gff_pathname) or os.path.isfile(gbk_pathname):
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
        elif all([True for species_file in os.listdir(run_folder) if '.pf' in species_file or '.fasta' in species_file]):
            genetic_writer = csv.writer(genetic_file, delimiter='\t', lineterminator='\n')
            for species_file in os.listdir(run_folder):
                    if '.pf' in species_file:
                        species_file_name = os.path.splitext(species_file)[0]
                        if species_file_name in taxon_datas:
                            circular = taxon_datas[species_file_name][0]
                            element_type = taxon_datas[species_file_name][1]
                            codon_table = taxon_datas[species_file_name][2]
                        else:
                            if 'circular' in taxon_datas:
                                circular = taxon_datas['circular']
                            if 'element_type' in taxon_datas:
                                element_type = taxon_datas['element_type']
                            if 'codon_table' in taxon_datas:
                                codon_table = taxon_datas['codon_table']
                        genetic_writer.writerow(['NAME', species_file.replace('.pf', '')])
                        genetic_writer.writerow(['ID', species_file.replace('.pf', '')])
                        genetic_writer.writerow(['ANNOT-FILE', species_file])
                        genetic_writer.writerow(['SEQ-FILE', species_file.replace('.pf', '.fasta')])
                        if circular:
                            genetic_writer.writerow(['CIRCULAR?', circular])
                        if element_type:
                            genetic_writer.writerow(['TYPE', element_type])
                        if codon_table:
                            genetic_writer.writerow(['CODON-TABLE', codon_table])
                        genetic_writer.writerow(['//'])
    # Create the lisp script.
    check_lisp_file = create_dat_creation_script(pgdb_id, lisp_pathname)

    return all([os.path.isfile(organism_dat), os.path.isfile(genetic_dat), check_lisp_file])


def pwt_input_files(multiprocess_input):
    """
    Check if files needed by Pathway Tools are available, if not create them.
    Check if there is a pathologic.log from a previous run. If yes, delete it.

    Args:
        multiprocess_input (dict): multiprocess dictionary input
    """
    run_folder = multiprocess_input['species_input_folder_path']
    taxon_file = multiprocess_input['taxon_file']

    required_files = set(['organism-params.dat', 'genetic-elements.dat', 'dat_creation.lisp'])
    files_in = set(next(os.walk(run_folder))[2])

    species_folder = run_folder.split('/')[-2]

    if 'pathologic.log' in files_in:
        os.remove(run_folder + 'pathologic.log')

    missing_string = ''
    if required_files.issubset(files_in):
        missing_string = 'no missing files'
    else:
        missing_string = 'missing {0}'.format('; '.join(required_files.difference(files_in))) + '. Inputs file created for {0}'.format(run_folder.split('/')[-2])
        check_datas_lisp = create_dats_and_lisp(run_folder, taxon_file)
        if not check_datas_lisp:
            raise Exception('Error with the creation of input files of {0}'.format(run_folder))

    logger.info('Checking inputs for {0}: {1}.'.format(species_folder, missing_string))


def create_mpwt_input(run_ids, input_folder, pgdbs_folder_path,
                      patho_hole_filler=None, dat_extraction=None, output_folder=None,
                      size_reduction=None, only_dat_creation=None, taxon_file=None):
    """
    Create input list for all multiprocess function, containing one lsit for each input subfolder.
    All arguments are also stored.

    Args:
        run_ids (list): input species IDs
        input_folder (str): pathname to input folder
        pgdbs_folder_path (str): pathname to species PGDB in ptools-local
        patho_hole_filler (bool): PathoLogic Hole Filler argument
        dat_extraction (bool): BioPAX/attribute-values file extraction argument
        output_folder (str): pathname to output folder
        size_reduction (bool): ptools-local PGDB deletion after processing argument
        only_dat_creation (bool): only create BioPAX/attribute values argument
    Returns:
        dictionary: contain all these data for multiprocessing
    """
    multiprocess_inputs = []
    for run_id in run_ids:
        multiprocess_input = {}
        input_folder_path = input_folder + '/' + run_id + '/'
        species_pgdb_folder = pgdbs_folder_path + run_id.lower() + 'cyc/'
        pgdb_id_folders = (run_id, species_pgdb_folder)
        if only_dat_creation:
            multiprocess_input['pgdb_folders'] = retrieve_complete_id(pgdb_id_folders)
        else:
            multiprocess_input['pgdb_folders'] = pgdb_id_folders
        multiprocess_input['species_input_folder_path'] = input_folder_path
        multiprocess_input['patho_hole_filler'] = patho_hole_filler
        multiprocess_input['dat_extraction'] = dat_extraction
        multiprocess_input['output_folder'] = output_folder
        multiprocess_input['size_reduction'] = size_reduction
        multiprocess_input['taxon_file'] = taxon_file
        multiprocess_inputs.append(multiprocess_input)

    return multiprocess_inputs


def create_only_dat_lisp(pgdbs_folder_path, tmp_folder):
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
        if os.path.isdir(pgdbs_folder_path + species_pgdb):
            pgdb_id = species_pgdb[:-3]
            pgdb_pathname = tmp_folder + pgdb_id + '/'
            os.mkdir(tmp_folder + pgdb_id)
            lisp_pathname = pgdb_pathname + "dat_creation.lisp"
            check_lisp_file = create_dat_creation_script(pgdb_id, lisp_pathname)
            if not check_lisp_file:
                raise Exception('Error with the creation of the lisp script for {0}'.format(species_pgdb))

            yield pgdb_id


def retrieve_complete_id(pgdb_id_folder):
    """
    Retrieve the ID of the PGDB from the genetic-elements.dat file.

    Args:
        pgdb_id_folder (tuple): second tuple argument is the pathname to the PGDB
    Returns:
        tuple: (new PGDB ID (according to input file), pathname to PGDB folder)
    """
    with open(pgdb_id_folder[1] + '/1.0/input/genetic-elements.dat') as organism_file:
        for line in organism_file:
            if 'ANNOT-FILE' in line and ';;' not in line:
                pgdb_id_complete = line.split('\t')[1].replace('.gff','').replace('.gbk','').strip()

    return (pgdb_id_complete, pgdb_id_folder[1])
