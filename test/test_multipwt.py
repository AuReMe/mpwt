#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test mpwt on a genbank file containing E. coli genes implied in the TCA cycle.
"""

import mpwt

def test_multiprocess_pwt():
    mpwt.multiprocess_pwt('test', 'test_output', dat_extraction=True, size_reduction=False, verbose=True)

    pathway_pathname = "test_output/tca_cycle_ecoli/pathways.dat"

    tca_reactions = ["RXN-14971", "MALATE-DEH-RXN", "ISOCITDEH-RXN", "MALATE-DEHYDROGENASE-ACCEPTOR-RXN",
                    "ACONITATEDEHYDR-RXN", "CITSYN-RXN", "ACONITATEHYDR-RXN", "2OXOGLUTARATEDEH-RXN", "SUCCCOASYN-RXN", "FUMHYDR-RXN"]
    reactions = []
    with open(pathway_pathname) as pathway_file:
        for line in pathway_file:
            if 'REACTION-LIST' in line:
                reaction = line.split(' - ')[1].strip()

                reactions.append(reaction)

    assert set(tca_reactions).issubset(set(reactions))

test_multiprocess_pwt()