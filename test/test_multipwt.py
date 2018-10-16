#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test for mpwt data using gebank file from Buchnera aphidicola str. APS (Acyrthosiphon pisum).
"""

import gzip
import mpwt
import os
import shutil
import urllib.request

def test_on_genome():
    os.makedirs('test_data/baphidicola')

    gbk_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/605/GCA_000009605.1_ASM960v1/GCA_000009605.1_ASM960v1_genomic.gbff.gz'
    compressed_file_pathname = 'test_data/baphidicola/baphidicola.gbk.gz'
    gbk_file_pathname = 'test_data/baphidicola/baphidicola.gbk'
    run_folder = 'test_data'
    # Retrieve Buchnera aphidicola genome.
    urllib.request.urlretrieve(gbk_url, compressed_file_pathname)

    # Extract Buchnera aphidicola genome.
    with gzip.open(compressed_file_pathname, 'rb') as compressed_file:
        with open(gbk_file_pathname, 'wb') as gbk_file:
            shutil.copyfileobj(compressed_file, gbk_file)

    os.remove(compressed_file_pathname)

    # Run mpwt on genome.
    mpwt.multiprocess_pwt(run_folder, output_folder=None, dat_extraction=True, verbose=True)