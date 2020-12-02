#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test mpwt functions which don't need a Pathway-Tools environment.
"""

import mpwt
import os


def test_create_dats_and_lisp():
    input_path = os.path.join(*['test', 'fatty_acid_beta_oxydation_I'])
    mpwt.pathologic_input.create_flats_and_lisp(input_path, False)

    genetic_pathname = os.path.join(*['test', 'fatty_acid_beta_oxydation_I', 'genetic-elements.dat'])
    organism_pathname = os.path.join(*['test', 'fatty_acid_beta_oxydation_I', 'organism-params.dat'])
    lisp_pathname = os.path.join(*['test', 'fatty_acid_beta_oxydation_I' , 'flat_files_creation.lisp'])

    genetic_string_expected = 'NAME\t\nANNOT-FILE\tfatty_acid_beta_oxydation_I.gbk\n//\n'
    organism_string_expected = 'ID\tfatty_acid_beta_oxydation_I\nSTORAGE\tFILE\nNCBI-TAXON-ID\t511145\nNAME\tEscherichia coli str. K-12 substr. MG1655\n'
    lisp_string_expected = ''';; Most files should begin by specifying that they will be interpreted\n;; within the Lisp package called ecocyc.\n(in-package :ecocyc)\n;; Select the organism PGDB as the current PGDB\n(select-organism :org-id 'fatty_acid_beta_oxydation_I)\n;; Create attribute-values files without the progression pop-up\n(let ((*progress-noter-enabled?* NIL))\n        (create-flat-files-for-current-kb))'''

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


def test_check_input_and_existing_pgdb():
    mpwt.remove_pgdbs('fatty_acid_beta_oxydation_icyc,fatty_acid_beta_oxydation_i_gffcyc,fatty_acid_beta_oxydation_i_pfcyc')

    run_patho_ids, run_flat_ids = mpwt.pathologic_input.check_input_and_existing_pgdb(['fatty_acid_beta_oxydation_I', 'fatty_acid_beta_oxydation_I_pf', 'fatty_acid_beta_oxydation_I_gff'], 'test', None, 1)
    assert sorted(run_patho_ids) == sorted(['fatty_acid_beta_oxydation_I', 'fatty_acid_beta_oxydation_I_pf', 'fatty_acid_beta_oxydation_I_gff'])
    assert not run_flat_ids

    run_patho_ids, run_flat_ids = mpwt.pathologic_input.check_input_and_existing_pgdb(['wrong.gbk.id.gbk'], 'test_wrong_id', None, 1)
    assert not run_patho_ids
    assert not run_flat_ids
