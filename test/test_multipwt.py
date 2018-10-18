#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test mpwt on a genbank file containing E. coli genes implied in the TCA cycle.
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

def test_multiprocess_pwt():
    mpwt.multiprocess_pwt('test', 'test_output', dat_extraction=True, size_reduction=False, verbose=True)

    pathway_pathname = "test_output/tca_cycle_ecoli/pathways.dat"

    tca_reactions = ["RXN-14971", "MALATE-DEH-RXN", "ISOCITDEH-RXN", "MALATE-DEHYDROGENASE-ACCEPTOR-RXN",
                    "ACONITATEDEHYDR-RXN", "CITSYN-RXN", "ACONITATEHYDR-RXN", "2OXOGLUTARATEDEH-RXN", "SUCCCOASYN-RXN", "FUMHYDR-RXN"]
    reactions = []
    with open(pathway_pathname, 'r') as pathway_file:
        for line in pathway_file:
            if 'REACTION-LIST' in line:
                reaction = line.split(' - ')[1].strip()

                reactions.append(reaction)

    assert set(tca_reactions).issubset(set(reactions))

test_create_dats_and_lisp()
test_multiprocess_pwt()
