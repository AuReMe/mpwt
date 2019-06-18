#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From Genbank/GFF/PF files this script will create Pathway Tools input data, then run Pathway Tools's PathoLogic on them. It can also generate dat files.
The script takes a folder name as argument.

usage:
    mpwt -f=DIR [-o=DIR] [--patho] [--hf] [--dat] [--md] [--cpu=INT] [-r] [-v] [--clean] [--log=FOLDER]
    mpwt --dat [-f=DIR] [-o=DIR] [--md] [--cpu=INT] [-v]
    mpwt -o=DIR [--md] [--cpu=INT] [-v]
    mpwt --clean [--cpu=INT] [-v]
    mpwt --delete=STR [--cpu=INT]
    mpwt --list

options:
    -h --help     Show help.
    -f=DIR     Working folder containing sub-folders with Genbank/GFF/PF files.
    -o=DIR    Output folder path. Will create a output folder in this folder.
    --patho    Will run an inference of Pathologic on the input files.
    --hf    Use with --patho. Run the Hole Filler using Blast.
    --dat    Will create BioPAX/attribute-value dat files from PGDB.
    --md    Move only the dat files into the output folder.
    --clean    Clean ptools-local folder, before any other operations.
    --delete=STR    Give a PGDB name and it will delete it (if multiple separe them with a ',', example: ecolicyc,athalianacyc).
    -r    Will delete files in ptools-local and compress results files to reduce results size (use it with -o).
    --cpu=INT     Number of cpu to use for the multiprocessing (default=1).
    --log=FOLDER     Create PathoLogic log files inside the given folder (use it with --patho).
    --list     List all PGDBs inside the ptools-local folder.
    -v     Verbose.

"""

import docopt
import logging
import os
import sys

from mpwt import utils
from mpwt.mpwt_workflow import multiprocess_pwt
from multiprocessing import Pool

logging.basicConfig(format='%(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def run_mpwt():
    """
    Function used with a mpwt call in the terminal.
    """
    args = docopt.docopt(__doc__)

    argument_number = len(sys.argv[1:])

    input_folder = args['-f']
    output_folder = args['-o']
    patho_inference = args['--patho']
    patho_hole_filler = args['--hf']
    dat_creation = args['--dat']
    move_dat = args['--md']
    size_reduction = args['-r']
    number_cpu = args['--cpu']
    patho_log = args['--log']
    pgdb_to_deletes = args['--delete']
    pgdb_list = args['--list']
    verbose = args['-v']

    if pgdb_list:
        pgdbs = utils.list_pgdb()
        if pgdbs == []:
            logger.info('No PGDB inside ptools-local.')
        else:
            logger.info(str(len(pgdbs)) + ' PGDB inside ptools-local:\n' + '\t'.join(pgdbs))
        return

    # Delete PGDB if use of --delete argument.
    # Use a set to remove redudant PGDB.
    if pgdb_to_deletes:
        utils.remove_pgbds(list(set(pgdb_to_deletes.split(','))), number_cpu)
        return

    if args['--clean']:
        if verbose:
            logger.info('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')
        utils.cleaning(number_cpu, verbose)
        if input_folder:
            utils.cleaning_input(input_folder, verbose)
        if argument_number == 1 or (argument_number == 2 and verbose) or (argument_number == 4 and verbose and number_cpu):
            sys.exit()

    multiprocess_pwt(input_folder,
                    output_folder,
                    patho_inference,
                    patho_hole_filler,
                    dat_creation,
                    move_dat,
                    size_reduction,
                    number_cpu,
                    patho_log,
                    verbose)


if __name__ == "__main__":
    run_mpwt()
