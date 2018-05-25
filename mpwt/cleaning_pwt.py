#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

from mpwt.multipwt import check_existing_pgdb

def cleaning():
    """
    Clean Pathway-Tools PGDB's folder.
    The script will delete folders and files in ptools-local/pgdbs/user.

    """
    pathway_tools_str = subprocess.check_output('type pathway-tools', shell=True)
    pathway_tools_path = pathway_tools_str.decode('UTF-8').split('is ')[1].strip('\n')

    pathway_tools_file = open(pathway_tools_path, 'r')   
    ptools_local_str = [line for line in pathway_tools_file if 'PTOOLS_LOCAL_PATH' in line][0]
    ptools_local_path = ptools_local_str.split(';')[0].split('=')[1].replace('"', '').strip(' ') + '/ptools-local'
    file_path = ptools_local_path.replace('\n', '') +'/pgdbs/user/'
    pathway_tools_file.close()

    pgdb_metadata_path = file_path + 'PGDB-METADATA.ocelot'
    if os.path.isfile(pgdb_metadata_path):
        os.remove(pgdb_metadata_path)
        print('PGDB-METADATA.ocelot has been removed.')

    pgdb_counter_path = file_path + 'PGDB-counter.dat'
    if os.path.isfile(pgdb_counter_path):
        os.remove(pgdb_counter_path)
        print('PGDB-counter.dat has been removed.')

    for pgdb_folder in os.listdir(file_path):
        pgdb_folder_path = file_path + pgdb_folder
        shutil.rmtree(pgdb_folder_path)
        print(pgdb_folder + ' has been removed.')

def cleaning_input(folder, output_folder=None):
    run_ids = [folder_id for folder_id in next(os.walk(folder))[1]]
    if output_folder is not None:
        run_ids = check_existing_pgdb(run_ids, output_folder)
    genbank_paths = [folder + "/" + run_id + "/" for run_id in run_ids]
    for genbank_path in genbank_paths:
        os.remove(genbank_path + 'script.lisp')
        os.remove(genbank_path + 'pathologic.log')
        os.remove(genbank_path + 'genetic-elements.dat')
        os.remove(genbank_path + 'organism-params.dat')
        print('Remove temporary datas.')