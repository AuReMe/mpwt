#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

from mpwt.multipwt import check_existing_pgdb, ptools_path

def cleaning(verbose=None):
    """
    Clean Pathway-Tools PGDB's folder.
    The script will delete folders and files in ptools-local/pgdbs/user.
    """
    ptools_local_path = ptools_path()
    file_path = ptools_local_path.replace('\n', '') +'/pgdbs/user/'

    pgdb_metadata_path = file_path + 'PGDB-METADATA.ocelot'
    if os.path.isfile(pgdb_metadata_path):
        os.remove(pgdb_metadata_path)
        if verbose:
            print('PGDB-METADATA.ocelot has been removed.')

    pgdb_counter_path = file_path + 'PGDB-counter.dat'
    if os.path.isfile(pgdb_counter_path):
        os.remove(pgdb_counter_path)
        if verbose:
            print('PGDB-counter.dat has been removed.')

    for pgdb_folder in os.listdir(file_path):
        pgdb_folder_path = file_path + pgdb_folder
        shutil.rmtree(pgdb_folder_path)
        if verbose:
            print(pgdb_folder + ' has been removed.')

def delete_pgdb(pgdb_name):
    """
    Remove a specific PGDB.
    """
    ptools_local_path = ptools_path()
    pgdb_path = ptools_local_path.replace('\n', '') +'/pgdbs/user/' + pgdb_name

    shutil.rmtree(pgdb_path)

    print(pgdb_name + ' (at ' + pgdb_path + ') has been removed.')

def cleaning_input(input_folder, output_folder=None, verbose=None):
    """
    Remove script.lisp, pathologic.log, genetic-elements.dat and organism-params.dat in a genbank folder.
    """
    run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]
    if output_folder:
        if os.path.exists(output_folder) == False:
            if verbose:
                print('No output directory, it will be created.')
            os.mkdir(output_folder)
        run_ids = check_existing_pgdb(run_ids, input_folder, output_folder)
    genbank_paths = [input_folder + "/" + run_id + "/" for run_id in run_ids]
    for genbank_path in genbank_paths:
        lisp_script = genbank_path + 'script.lisp'
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
        if verbose:
            species = genbank_path.split('/')[-2]
            print('Remove ' + species + ' temporary datas.')
