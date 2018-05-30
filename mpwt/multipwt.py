#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From genbank file this script will create Pathway-Tools input data, then run Pathway-Tools's PathoLogic on them and at last it will generate dat files for AuReMe.
The script takes a folder name as argument. The folder structure expected is:
Folder_input
├── Folder_for_species_1
│   └── Genbank_species_1
├── Folder_for_species_2
│   └── Genbank_species_2
├── Folder_for_species_3
│   └── Genbank_species_3
│
"""

import argparse
import csv
import datetime
import getpass
import os
import shutil
import subprocess
import sys

from Bio import SeqIO
from multiprocessing import Pool, cpu_count

def run():
    from mpwt.cleaning_pwt import cleaning, cleaning_input

    parser = argparse.ArgumentParser(usage="python pathway_tools_multiprocess.py -f FOLDER")
    parser.add_argument("-f", "--folder", dest = "folder", metavar = "FOLDER", help = "Folder containing sub-folders with Genbank file.")
    parser.add_argument("-o", "--output", dest = "output", metavar = "FOLDER", help = "Output folder path. Will create a output folder in this folder.", default=None)
    parser.add_argument("-d", "--dat", dest = "extract_dat", help = "Will extract only dat files from Pathway-Tools results.", action='store_true', default=None)
    parser.add_argument("-r", "--reduce", dest = "reduce_size", help = "Will delete files in ptools-local to reduce size of results.", action='store_true', default=None)
    parser.add_argument("clean", nargs='?', help = "Arguments to clean ptools-local folder.")

    parser_args = parser.parse_args(sys.argv[1:])

    input_folder = parser_args.folder
    output_folder = parser_args.output
    dat_extraction = parser_args.extract_dat
    size_reduction = parser_args.reduce_size

    if parser_args.clean:
        print('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')
        cleaning()
        if input_folder:
            cleaning_input(input_folder, output_folder)
        if len(sys.argv) == 2:
            sys.exit()

    multiprocess_pwt(input_folder, output_folder, dat_extraction,size_reduction)

def multiprocess_pwt(input_folder,output_folder=None,dat_extraction=None,size_reduction=None):
    # Run folder contains sub-folders containing GBK file
    run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
    if output_folder is not None:
        if os.path.exists(output_folder) == False:
            print('No output directory, it will be created.')
            os.mkdir(output_folder)
        run_ids = check_existing_pgdb(run_ids, output_folder)
    genbank_paths = [input_folder + "/" + run_id + "/" for run_id in run_ids]
    if len(genbank_paths) == 0:
        sys.exit("No folder containing genbank file. In " + input_folder + " you must have sub-folders containing Genbank file.")
    p = Pool(processes=cpu_count())
    print('~~~~~~~~~~Creation of input data from Genbank~~~~~~~~~~')
    for genbank_path in genbank_paths:
        pwt_run(genbank_path)
    print('~~~~~~~~~~Inference on the data~~~~~~~~~~')
    p.map(run_pwt, genbank_paths)
    print('~~~~~~~~~~Check inference~~~~~~~~~~')
    check_pwt(genbank_paths)
    print('~~~~~~~~~~Creation of the PGDB-METADATA.ocelot file~~~~~~~~~~')
    pgdb_folders = {}
    for genbank_path in genbank_paths:
        pgdb_id_folder = create_metadata(genbank_path)
        pgdb_folders[genbank_path] = pgdb_id_folder
    print('~~~~~~~~~~Creation of the .dat files~~~~~~~~~~')
    p.map(run_pwt_dat, genbank_paths)
    print('~~~~~~~~~~End of the Pathway-Tools Inference~~~~~~~~~~')
    print('~~~~~~~~~~Moving result files~~~~~~~~~~')
    for genbank_path in pgdb_folders:
        move_pgdb(genbank_path, pgdb_folders[genbank_path], output_folder, dat_extraction,size_reduction)
    print('~~~~~~~~~~The script have finished! Thank you for using it.')

def check_existing_pgdb(run_ids, output_folder):
    """
    Check output folder for already existing PGDB, don't create them.
    """
    already_present_pgdbs = [output_pgdb for output_pgdb in os.listdir(output_folder)]

    new_run_ids = set(run_ids) - set(already_present_pgdbs)

    new_run_ids = list(new_run_ids)

    if len(new_run_ids) == 0:
        sys.exit("All PGDBs are already present in the output folder. Remove them if you want a new inference.")

    return new_run_ids

def pwt_run(run_folder):
    """
    Check if files needed by Pathway-Tools are available, if not create them.
    """
    required_files = set(['organism-params.dat','genetic-elements.dat','script.lisp'])
    files_in = set(next(os.walk(run_folder))[2])
    print("Checking for pathwaytools inputs:")
    if required_files.issubset(files_in):
        print("OK")
    else:
        print("%s missing" %"; ".join(required_files.difference(files_in)))
        create_dats_and_lisp(run_folder)

def create_dats_and_lisp(run_folder):
    """
    Read Genbank file and create Pathway Tools needed file.
    Create also a lisp file to create dat files from Pathway tools results.
    The name of the PGDB created by Pathway Tools will be the name of the species with '_' instead of space.

    Genbank example:

    LOCUS       scaffold1         XXXXXX bp    DNA     linear   INV DD-MMM-YYYY
    DEFINITION  My species genbank.
    ACCESSION   scaffold1
    VERSION     scaffold1
    KEYWORDS    Key words.
    SOURCE      Source
    ORGANISM  Species name
                Taxonomy; Of; My; Species; With;
                The; Genus.
    FEATURES             Location/Qualifiers
        source          1..XXXXXX
                        /scaffold="scaffold1"
                        /db_xref="taxon:taxonid"

    Create organism-params.dat:
    ID  myDBName
    STORAGE FILE
    NCBI-TAXON-ID   taxon_id
    NAME    species_name

    Create genetic-elements.dats:
    NAME    
    ANNOT-FILE  gbk_name
    //

    Create script.lisp:
    (in-package :ecocyc)
    (select-organism :org-id 'myDBName)
    (create-flat-files-for-current-kb)
    """
    # Look for a Genbank file in the run folder.
    for input_file in os.listdir(run_folder):
        if input_file.endswith(".gbk") or input_file.endswith(".gb") or input_file.endswith(".gbff"):
            gbk_file = run_folder + input_file
            gbk_name = input_file

    # Check if a Genbank file have been found.
    try:
        gbk_file
    except NameError:
        raise NameError('Missing Genbank file. Check if you have a Genbank file and if it ends with .gbk or .gb or .gbff')

    organism_dat = run_folder + 'organism-params.dat'
    genetic_dat = run_folder + 'genetic-elements.dat'

    taxon_id = ""
    species_name = ""

    # Take the species name and the taxon id from the genbank file.
    with open(gbk_file, "rU") as gbk:
        # Take the first record of the genbank (first contig/chromosome) to retrieve the species name.
        first_seq_record = next(SeqIO.parse(gbk, "genbank"))
        try:
            species_name = first_seq_record.annotations['organism']
        except KeyError:
            raise KeyError('No organism in the Genbank. In the SOURCE you must have: ORGANISM  Species name')      

        # Take the source feature of the first record.
        # This feature contains the taxon ID in the db_xref qualifier.
        src_feature = [f for f in first_seq_record.features if f.type == "source"][0]
        try:
            taxon_id = src_feature.qualifiers['db_xref'][0].replace('taxon:', '')
        except KeyError:
            raise KeyError('No taxon ID in the Genbank. In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.')               

    lisp_file = run_folder + "script.lisp"

    # The name of the PGDB will be the name of the species.
    myDBName = gbk_name.split('.')[0]

    # Create the organism-params dat file.
    with open(organism_dat, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['ID', myDBName])
        writer.writerow(['STORAGE', "FILE"])
        writer.writerow(['NCBI-TAXON-ID', taxon_id])
        writer.writerow(['NAME', species_name])

    # Create the genetic-elements dat file.
    with open(genetic_dat, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['NAME', ''])
        writer.writerow(['ANNOT-FILE', gbk_name])
        writer.writerow(['//'])

    # Create the lisp script.
    with open(lisp_file, 'w') as file:
        file.write("(in-package :ecocyc)")
        file.write('\n')
        file.write("(select-organism :org-id '" + myDBName + ")")
        file.write('\n')
        file.write("(create-flat-files-for-current-kb)")

def check_pwt(genbank_paths):
    """
    Check PathoLogic's log.
    """
    failed_inferences = []
    passed_inferences = []

    with open('log_error.txt', 'w') as output_file:
        with open('resume_inference.tsv', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
            writer.writerow(['species', 'gene_number', 'protein_number', 'pathway_number', 'reaction_number', 'compound_number'])
            for genbank_path in genbank_paths:
                patho_log = genbank_path + '/pathologic.log'

                output_file.write('------------ Species: ')
                output_file.write(genbank_path)
                output_file.write('\n')

                fatal_error_index = None

                with open(patho_log, 'r') as input_file:
                    for index, line in enumerate(input_file):
                        if 'fatal error' in line:
                            fatal_error_index = index
                            output_file.write(line)
                            writer.writerow([genbank_path, 'ERROR', '', '', '', ''])
                            failed_inferences.append(genbank_path)
                        if fatal_error_index is not None:
                            if index > fatal_error_index:
                                output_file.write(line)
                        if 'Build done.' in  line:
                            output_file.write(line)
                            resume_inference_line = next(input_file)
                            output_file.write(resume_inference_line)
                            gene_number = int(resume_inference_line.split('PGDB contains ')[1].split(' genes')[0])
                            protein_number = int(resume_inference_line.split('genes, ')[1].split(' proteins')[0])
                            pathway_number = int(resume_inference_line.split('proteins, ')[1].split(' base pathways')[0])
                            reaction_number = int(resume_inference_line.split('base pathways, ')[1].split(' reactions')[0])
                            compound_number = int(resume_inference_line.split('reactions, ')[1].split(' compounds')[0])
                            writer.writerow([genbank_path, gene_number, protein_number, pathway_number, reaction_number, compound_number])

                            passed_inferences.append(genbank_path)

                output_file.write('------------\n\n')

    with open('log_error.txt','r') as contents:
        save = contents.read()
    with open('log_error.txt', 'w') as output_file:
            output_file.write('Inference statistics:\n')
            if len(passed_inferences) > 0:
                print('\n' + str(len(passed_inferences)) + ' builds have passed!\n')   
                output_file.write('Build done: ' + str(len(passed_inferences)))
                output_file.write('\tSpecies: ' + ', '.join(passed_inferences)+'\n')
            if len(failed_inferences) > 0:
                print('WARNING: ' + str(len(failed_inferences)) + ' builds have failed! See the log for more information.\n')
                output_file.write('Build failed: ' + str(len(failed_inferences)))
                output_file.write('\tSpecies: ' + ', '.join(failed_inferences)+'\n\n')
            output_file.write(save)
    
    if len(failed_inferences) > 0:
        sys.exit("Stop the inference.")

    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', 'log_error.txt'])
    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', 'resume_inference.tsv'])

def create_metadata(run_folder):
    """
    Create the PGDB-METADATA.ocelot file required for the creation of the dat files by Pathway-Tools (with run_pwt_dat).
    Because if we want the multiprocessing to work with Pathway-Tools we have to disable the metadata saving (because multiple process can't open and write this file at the same time).
    This will stop the creation of this file so after the creation of PGDB in multiprocess (run_pwt) the script will create these files.

    Create PGDB-METADAT.ocelot and add data for each species:
    (:OCELOT-KB :NAME PGDB-METADATA :FORMAT :V1-SEXPR :PACKAGE "ECOCYC" :WRITE-DATE
    hour_date :USER username :PROCESS-ID process_id)

    ( myDBName NIL (
    (ORGID myDBName )
    (OCELOT-GFP::PARENTS |PGDBs|)
    (ROOT-PATHNAME "file_path")
    (FULL-SPECIES-NAME "species_name")
    (SPECIES-NAME "species_name")
    (NCBI-TAX-ID "taxon_id")
    (RANK |species|)
    (PGDB-LABEL "species_name") )
    NIL)

    """
    for input_file in os.listdir(run_folder):
        if input_file.endswith(".gbk") or input_file.endswith(".gb") or input_file.endswith(".gbff"):
            gbk_file = run_folder + '/' + input_file
            gbk_name = input_file

    taxon_id = ""
    species_name = ""

    # Take the species name and the taxon id from the genbank file.
    with open(gbk_file, "rU") as gbk:
        # Take the first record of the genbank (first contig/chromosome) to retrieve the species name.
        first_seq_record = next(SeqIO.parse(gbk, "genbank"))
        try:
            species_name = first_seq_record.annotations['organism']
        except KeyError:
            raise KeyError('No organism in the Genbank. In the SOURCE you must have: ORGANISM  Species name')      

        # Take the source feature of the first record.
        # This feature contains the taxon ID in the db_xref qualifier.
        src_feature = [f for f in first_seq_record.features if f.type == "source"][0]
        try:
            taxon_id = src_feature.qualifiers['db_xref'][0].replace('taxon:', '')
        except KeyError:
            raise KeyError('No taxon ID in the Genbank. In the FEATURES source you must have: /db_xref="taxon:taxonid" Where taxonid is the Id of your organism. You can find it on the NCBI.')               

    # The name of the PGDB will be the name of the species.
    myDBName = gbk_name.split('.')[0]

    pathway_tools_str = subprocess.check_output('type pathway-tools', shell=True)
    pathway_tools_path = pathway_tools_str.decode('UTF-8').split('is ')[1].strip('\n')
    pathway_tools_file = open(pathway_tools_path, 'r')   
    ptools_local_str = [line for line in pathway_tools_file if 'PTOOLS_LOCAL_PATH' in line][0]
    ptools_local_path = ptools_local_str.split(';')[0].split('=')[1].replace('"', '').strip(' ') + '/ptools-local'
    file_path = ptools_local_path.replace('\n', '') +'/pgdbs/user/'

    pathway_tools_file.close()

    if os.path.exists(file_path+'PGDB-METADATA.ocelot') == False:
        current_time = datetime.datetime.now()
        hour_date = current_time.strftime("%H:%M:%S on %a %b %d, %Y")
        process_id = str(os.getpid())
        username = getpass.getuser()
        with open(file_path+'PGDB-METADATA.ocelot', 'w') as metadata_file:
            metadata_file.write('''(:OCELOT-KB :NAME PGDB-METADATA :FORMAT :V1-SEXPR :PACKAGE "ECOCYC" :WRITE-DATE
 "''' + hour_date + '''" :USER "''' + username + '''" :PROCESS-ID ''' + process_id + ''')\n''')

    with open(file_path+'PGDB-METADATA.ocelot', 'a') as metadata_file:
         metadata_file.write('(' + myDBName +' NIL (\n')
         metadata_file.write('(ORGID '+ myDBName +')\n')
         metadata_file.write('(OCELOT-GFP::PARENTS |PGDBs|)\n')
         metadata_file.write('(ROOT-PATHNAME "' + file_path + '")\n')
         metadata_file.write('(FULL-SPECIES-NAME "' + species_name + '")\n')
         metadata_file.write('(SPECIES-NAME "' + species_name + '")\n')
         metadata_file.write('(NCBI-TAX-ID "' + taxon_id + '")\n')
         metadata_file.write('(RANK |species|)\n')
         metadata_file.write('(PGDB-LABEL "' + species_name +'") )\n')
         metadata_file.write('NIL)')
         metadata_file.write('\n\n')

    # Replace all / by _ to ensure that there is no error with the path with myDBName.
    pgdb_folder = file_path + myDBName.replace('/', '_').lower() + 'cyc/'
    pgdb_id_folder = (myDBName, pgdb_folder)

    return pgdb_id_folder

def run_pwt(genbank_path):
    """
    Create PGDB using files created during 'create_dats_and_lisp'.
    """
    cmd_pwt = "pathway-tools -no-web-cel-overview -no-cel-overview -disable-metadata-saving -nologfile -patho %s" %genbank_path
    print(cmd_pwt)
    subprocess.call(cmd_pwt, shell=True)

def run_pwt_dat(genbank_path):
    """
    Create dat file using a lisp script created during 'create_dats_and_lisp'.
    Add an input to the subprocess call to close the Navigator Window opening proposition ('Enter name of X-window server to connect to (of the form HOST:N.M):').
    If this proposition is not closed the script can't continue.
    """
    cmd_dat = "pathway-tools -no-web-cel-overview -no-cel-overview -disable-metadata-saving -nologfile -load %s/script.lisp" %genbank_path
    print(cmd_dat)
    p = subprocess.Popen(cmd_dat, shell=True, stdin=subprocess.PIPE)
    p.communicate(input=b'none')

def move_pgdb(genbank_path, pgdb_folder, output_folder, dat_extraction, size_reduction):
    """
    Move the result files inside the shared folder containing the input data.
    """
    pgdb_folder_dbname = pgdb_folder[0]
    pgdb_folder_path = pgdb_folder[1]

    output_dat_path = pgdb_folder_path + '/1.0/data'

    if output_folder is None:
        output_species = genbank_path + '/output/'
    else:
        output_species = output_folder + '/' + pgdb_folder_dbname +'/'

    if dat_extraction == True:
        pgdb_folder_path = output_dat_path

    if size_reduction is True:
        for pgdb_file in os.listdir(pgdb_folder_path):
            if os.path.exists(output_species) == False:
                os.mkdir(output_species)
            if dat_extraction == True:
                if '.dat' in pgdb_file:
                    shutil.move(pgdb_folder_path+'/'+pgdb_file, output_species+pgdb_file)
            elif dat_extraction is None:
                shutil.move(pgdb_folder_path+'/'+pgdb_file, output_species+pgdb_file)
        shutil.rmtree(pgdb_folder_path)
    else:
        shutil.copytree(pgdb_folder_path, output_species)
        if dat_extraction == True:
            for pgdb_file in os.listdir(output_species):
                if '.dat' not in pgdb_file:
                    os.remove(output_species+'/'+pgdb_file)

    # Give access to the file for user outside the container.
    subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_species])
    if output_folder is not None:
        subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_folder])

if __name__ == '__main__':
    run()
