#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

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
