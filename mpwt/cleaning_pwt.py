#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

from mpwt.multipwt import check_input_and_existing_pgdb, find_ptools_path
from multiprocessing import Pool, cpu_count


def delete_pgdb(pgdb_name):
    """
    Remove a specific PGDB.
    """
    ptools_local_path = find_ptools_path()
    pgdb_path = ptools_local_path.replace('\n', '') +'/pgdbs/user/' + pgdb_name
    if os.path.isdir(pgdb_path):
        shutil.rmtree(pgdb_path)
        print('{0} (at {1}) has been removed.'.format(pgdb_name, pgdb_path))
    else:
        print(pgdb_path + " not a folder.")


def remove_pgbds(to_delete_pgdbs, number_cpu=None):
    """
    Delete all PGDB inside to_delete_pgdbs using multiprocessing.
    Check if there is a Pool and if not spawn one.
    """
    # Use the number of cpu given by the user or all the cpu available.
    if number_cpu:
        number_cpu_to_use = int(number_cpu)
    else:
        number_cpu_to_use = cpu_count()
    mpwt_pool = Pool(processes=number_cpu_to_use)

    mpwt_pool.map(delete_pgdb, to_delete_pgdbs)

    mpwt_pool.close()
    mpwt_pool.join()

def cleaning(number_cpu=None, verbose=None):
    """
    Clean Pathway-Tools PGDB's folder.
    The script will delete folders and files in ptools-local/pgdbs/user.
    """
    ptools_local_path = find_ptools_path()
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

    # Extract all pgdbs inside ptools-local. Then delete them.
    all_pgdbs = os.listdir(file_path)
    remove_pgbds(all_pgdbs)


def cleaning_input(input_folder, output_folder=None, verbose=None):
    """
    Remove dat_creation.lisp, pathologic.log, genetic-elements.dat and organism-params.dat in a genbank folder.
    """
    run_ids = [folder_id for folder_id in next(os.walk(input_folder))[1]]

    if output_folder:
        if os.path.exists(output_folder) == False:
            output_folder = None
        run_ids = check_input_and_existing_pgdb(run_ids, input_folder, output_folder, verbose)
        if not run_ids:
            return

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
            if verbose:
                species = genbank_path.split('/')[-2]
                print('Remove ' + species + ' temporary datas.')
