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


def test_multiprocess_pwt():
    mpwt.cleaning()
    mpwt.cleaning_input('test')
    mpwt.multiprocess_pwt('test', 'test_output', dat_extraction=True, size_reduction=False, verbose=True)

    pathway_tca_pathname = "test_output/tca_cycle_ecoli/pathways.dat"

    tca_reactions = ["RXN-14971", "MALATE-DEH-RXN", "ISOCITDEH-RXN", "MALATE-DEHYDROGENASE-ACCEPTOR-RXN",
                    "ACONITATEDEHYDR-RXN", "CITSYN-RXN", "ACONITATEHYDR-RXN", "2OXOGLUTARATEDEH-RXN", "SUCCCOASYN-RXN", "FUMHYDR-RXN"]
    expected_tca_reactions = []
    with open(pathway_tca_pathname, 'r') as pathway_file:
        for line in pathway_file:
            if 'REACTION-LIST' in line:
                reaction = line.split(' - ')[1].strip()
                expected_tca_reactions.append(reaction)

    assert set(tca_reactions).issubset(set(expected_tca_reactions))

    pathway_fabo_pathname = "test_output/fatty_acid_beta_oxydation_I/pathways.dat"

    fabo_reactions = ["ACYLCOASYN-RXN", "ACYLCOADEHYDROG-RXN", "ENOYL-COA-DELTA-ISOM-RXN", "ENOYL-COA-HYDRAT-RXN",
                    "OHBUTYRYL-COA-EPIM-RXN", "OHACYL-COA-DEHYDROG-RXN", "KETOACYLCOATHIOL-RXN"]
    expected_fabo_reactions = []
    with open(pathway_fabo_pathname, 'r') as pathway_file:
        for line in pathway_file:
            if 'REACTION-LIST' in line:
                reaction = line.split(' - ')[1].strip()
                expected_fabo_reactions.append(reaction)

    assert set(fabo_reactions).issubset(set(expected_fabo_reactions))

    os.remove('log_error.txt')
    mpwt.cleaning_input('test')
    os.remove('resume_inference.tsv')
    shutil.rmtree('test_output')
    shutil.rmtree('__pycache__')
