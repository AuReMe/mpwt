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
import logging

logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
logger = logging.getLogger(__name__)
logging.getLogger("mpwt").setLevel(logging.DEBUG)


def reaction_extraction(pathway_file_pathname):
    expected_reactions = []

    with open(pathway_file_pathname, 'r') as pathway_file:
        for line in pathway_file:
            if 'REACTION-LIST' in line:
                reaction = line.split(' - ')[1].strip()
                expected_reactions.append(reaction)

    return expected_reactions


def fabo_reactions():
    return ["ACYLCOASYN-RXN", "ACYLCOADEHYDROG-RXN", "ENOYL-COA-DELTA-ISOM-RXN", "ENOYL-COA-HYDRAT-RXN",
                    "OHBUTYRYL-COA-EPIM-RXN", "OHACYL-COA-DEHYDROG-RXN", "KETOACYLCOATHIOL-RXN"]


def test_multiprocess_pwt_import():
    """
    Test mpwt when called in a python script.
    """
    mpwt.remove_pgdbs('fatty_acid_beta_oxydation_icyc,fatty_acid_beta_oxydation_i_gffcyc,fatty_acid_beta_oxydation_i_pfcyc')
    mpwt.cleaning_input('test')

    mpwt.create_pathologic_file('test', 'test_pf')
    mpwt.multiprocess_pwt('test_pf', 'test_output', patho_inference=True, flat_creation=True, dat_extraction=True, size_reduction=False, number_cpu=3, verbose=True)

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I_gff', 'pathways.dat'])
    expected_tca_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_tca_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I', 'pathways.dat'])
    expected_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I_pf', 'pathways.dat'])
    expected_pf_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_pf_fabo_reactions))

    mpwt.cleaning_input('test')
    shutil.rmtree('test_pf')
    shutil.rmtree('test_output')
    shutil.rmtree('__pycache__')


def test_multiprocess_pwt_call():
    """
    Test mpwt when called in terminal.
    """
    subprocess.call(['mpwt', '--delete', 'fatty_acid_beta_oxydation_icyc,fatty_acid_beta_oxydation_i_gffcyc,fatty_acid_beta_oxydation_i_pfcyc'])
    subprocess.call(['mpwt', '-f', 'test', '--clean'])
    subprocess.call(['mpwt', '-f', 'test', '--patho'])

    subprocess.call(['mpwt', '-o', 'test_output', '--flat', '--md', '--mx', '--mc', '--mo', '--cpu', '3'])

    pgdbs = mpwt.list_pgdb()
    assert set(['fatty_acid_beta_oxydation_i_gffcyc', 'fatty_acid_beta_oxydation_i_pfcyc', 'fatty_acid_beta_oxydation_icyc']).issubset(set(pgdbs))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_i_gff', 'pathways.dat'])
    expected_tca_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_tca_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_i', 'pathways.dat'])
    expected_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_i_pf', 'pathways.dat'])
    expected_pf_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_pf_fabo_reactions))

    subprocess.call(['mpwt', '-f', 'test', '--clean'])
    shutil.rmtree('test_output')
