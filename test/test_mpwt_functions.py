#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test mpwt functions which don't need a Pathway-Tools environment.
"""

import mpwt

def test_create_dats_and_lisp():
    mpwt.multipwt.create_dats_and_lisp('test/tca_cycle_ecoli/')

    genetic_pathname = 'test/tca_cycle_ecoli/genetic-elements.dat'
    organism_pathname = 'test/tca_cycle_ecoli/organism-params.dat'
    lisp_pathname = 'test/tca_cycle_ecoli/script.lisp'

    genetic_string_expected = 'NAME\t\nANNOT-FILE\ttca_cycle_ecoli.gbk\n//\n'
    organism_string_expected = 'ID\ttca_cycle_ecoli\nSTORAGE\tFILE\nNCBI-TAXON-ID\t511145\nNAME\tEscherichia coli str. K-12 substr. MG1655\n'
    lisp_string_expected = '''(in-package :ecocyc)\n(select-organism :org-id 'tca_cycle_ecoli)\n(create-flat-files-for-current-kb)'''

    with open(genetic_pathname, 'r') as genetic_file:
        genetic_string_found = genetic_file.read()
        assert genetic_string_found == genetic_string_expected

    with open(organism_pathname, 'r') as organism_file:
        organism_string_found = organism_file.read()
        assert organism_string_found == organism_string_expected

    with open(lisp_pathname, 'r') as lisp_file:
        lisp_string_found = lisp_file.read()
        assert lisp_string_found == lisp_string_expected

    mpwt.cleaning_input('test')


def test_check_existing_pgdb():
    mpwt.delete_pgdb('fatty_acid_beta_oxydation_icyc')
    mpwt.delete_pgdb('tca_cycle_ecolicyc')
    run_ids = mpwt.multipwt.check_existing_pgdb(['fatty_acid_beta_oxydation_I', 'tca_cycle_ecoli'], 'test', None)
    assert sorted(run_ids) == sorted(['fatty_acid_beta_oxydation_I', 'tca_cycle_ecoli'])

    run_ids = mpwt.multipwt.check_existing_pgdb(['wrong.gbk.id.gbk'], 'test_wrong_id', None)
    assert run_ids == None


def test_ptools_path():
    ptools_path = mpwt.multipwt.ptools_path()
    assert ptools_path == '/root/ptools-local'


def test_extract_pgdb_pathname():
    pgdb_id_folder = mpwt.multipwt.extract_pgdb_pathname('test/fatty_acid_beta_oxydation_I/')
    assert pgdb_id_folder == ('fatty_acid_beta_oxydation_I', '/root/ptools-local/pgdbs/user/fatty_acid_beta_oxydation_icyc/')
