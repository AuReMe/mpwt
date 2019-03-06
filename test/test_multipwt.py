#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test mpwt on genbank files containing E. coli genes implied in the TCA cycle and in the Fatty Acid Beta oxydation.
Need an environment with Pathway-Tools installed.
"""

import mpwt
import os
import shutil
import subprocess


def reaction_extraction(pathway_file_pathname):
    expected_reactions = []

    with open(pathway_file_pathname, 'r') as pathway_file:
        for line in pathway_file:
            if 'REACTION-LIST' in line:
                reaction = line.split(' - ')[1].strip()
                expected_reactions.append(reaction)

    return expected_reactions


def tca_reactions():
    return ["RXN-14971", "MALATE-DEH-RXN", "ISOCITDEH-RXN", "MALATE-DEHYDROGENASE-ACCEPTOR-RXN",
                    "ACONITATEDEHYDR-RXN", "CITSYN-RXN", "ACONITATEHYDR-RXN", "2OXOGLUTARATEDEH-RXN", "SUCCCOASYN-RXN", "FUMHYDR-RXN"]


def fabo_reactions():
    return ["ACYLCOASYN-RXN", "ACYLCOADEHYDROG-RXN", "ENOYL-COA-DELTA-ISOM-RXN", "ENOYL-COA-HYDRAT-RXN",
                    "OHBUTYRYL-COA-EPIM-RXN", "OHACYL-COA-DEHYDROG-RXN", "KETOACYLCOATHIOL-RXN"]


def test_multiprocess_pwt_import():
    """
    Test mpwt when called in a python script.
    """
    mpwt.remove_pgbds('fatty_acid_beta_oxydation_icyc,tca_cycle_ecolicyc')
    mpwt.cleaning_input('test')
    mpwt.multiprocess_pwt('test', 'test_output', patho_inference=True, dat_creation=True, dat_extraction=True, size_reduction=False, verbose=True)

    pathway_tca_pathname = "test_output/tca_cycle_ecoli/pathways.dat"
    expected_tca_reactions = reaction_extraction(pathway_tca_pathname)
    assert set(tca_reactions()).issubset(set(expected_tca_reactions))

    pathway_fabo_pathname = "test_output/fatty_acid_beta_oxydation_I/pathways.dat"
    expected_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    mpwt.cleaning_input('test')
    shutil.rmtree('test_output')
    shutil.rmtree('__pycache__')


def test_multiprocess_pwt_call():
    """
    Test mpwt when called in terminal.
    """
    subprocess.call(['mpwt', '--clean'])
    subprocess.call(['mpwt', '-f', 'test', '--clean'])
    subprocess.call(['mpwt', '-f', 'test', '--patho'])

    subprocess.call(['mpwt', '-o', 'test_output', '--dat', '--md'])

    pgdbs = mpwt.list_pgdb()
    assert sorted(pgdbs) == ['fatty_acid_beta_oxydation_icyc', 'tca_cycle_ecolicyc']

    pathway_tca_pathname = "test_output/tca_cycle_ecoli/pathways.dat"
    expected_tca_reactions = reaction_extraction(pathway_tca_pathname)
    assert set(tca_reactions()).issubset(set(expected_tca_reactions))

    pathway_fabo_pathname = "test_output/fatty_acid_beta_oxydation_I/pathways.dat"
    expected_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    subprocess.call(['mpwt', '-f', 'test', '--clean'])
    shutil.rmtree('test_output')
 