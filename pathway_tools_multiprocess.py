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
import docopt
import getpass
import os
import shutil
import subprocess
import sys

from Bio import SeqIO
from multiprocessing import Pool, cpu_count

parser = argparse.ArgumentParser(usage="python pathway_tools_multiprocess.py -f FOLDER")
parser.add_argument("-f", "--folder", dest = "folder", metavar = "FOLDER", help = "Folder containing sub-folders with Genbank file.")
parser.add_argument("-o", "--output", dest = "output", metavar = "FOLDER", help = "Output folder path. Will create a output folder in this folder.",default=None)

parser_args = parser.parse_args(sys.argv[1:])

def main():
    folder = parser_args.folder
    output_folder = parser_args.output
    # Run folder contains sub-folders containing GBK file
    run_ids = [folder_id for folder_id in next(os.walk(folder))[1]]
    genbank_paths = [folder + "/" + run_id + "/" for run_id in run_ids]
    p = Pool(processes=cpu_count())
    print('~~~~~~~~~~Creation of input data from Genbank~~~~~~~~~~')
    for genbank_path in genbank_paths:
        pwt_run(genbank_path)
    print('~~~~~~~~~~Inference on the data~~~~~~~~~~')
    p.map(run_pwt, genbank_paths)
    print('~~~~~~~~~~Creation of the PGDB-METADATA.ocelot file~~~~~~~~~~')
    pgdb_folders = {}
    for genbank_path in genbank_paths:
        pgdb_folder, species_name = create_metadata(genbank_path)
        pgdb_folders[genbank_path] = (pgdb_folder, species_name)
    print('~~~~~~~~~~Creation of the .dat files~~~~~~~~~~')
    p.map(run_pwt_dat, genbank_paths)
    print('~~~~~~~~~~End of the Pathway-Tools Inference~~~~~~~~~~')
    print('~~~~~~~~~~Moving result files~~~~~~~~~~')
    for genbank_path in pgdb_folders:
        move(genbank_path, pgdb_folders[genbank_path], output_folder)
    print('~~~~~~~~~~The script have finished! Thank you for using it.')
    #[None for _ in resultats]

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
    for file in os.listdir(run_folder):
        if file.endswith(".gbk") or file.endswith(".gb") or file.endswith(".gbff"):
            gbk_file = run_folder + file
            gbk_name = file

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
    myDBName = species_name.replace(' ', '_').replace('/', '_')

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
    for file in os.listdir(run_folder):
        if file.endswith(".gbk") or file.endswith(".gb") or file.endswith(".gbff"):
            gbk_file = run_folder + '/' + file

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
    myDBName = species_name.replace(' ', '_')

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
    species_name = species_name.replace('/', '_')

    return pgdb_folder, species_name

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

def move(genbank_path, pgdb_folder, output_folder):
    """
    Move the result files inside the shared folder containing the input data.
    """
    if output_folder is None:
        output_folder_path = genbank_path + 'output/'
        if not os.path.exists(output_folder_path):
            os.makedirs(output_folder_path)
        # Give access to the file for user outside the container.
        subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_folder_path])

    else:
        output_folder_path = output_folder + '/' + pgdb_folder[1] + '_output/'
        if not os.path.exists(output_folder_path):
            os.makedirs(output_folder_path)
        # Give access to the file for user outside the container.
        subprocess.call(['chmod', '-R', 'u=rwX,g=rwX,o=rwX', output_folder])

    shutil.move(pgdb_folder[0], output_folder_path)
if __name__ == "__main__":
    main()
