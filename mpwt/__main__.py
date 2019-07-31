#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From Genbank/GFF/PF files this script will create Pathway Tools input data, then run Pathway Tools's PathoLogic on them. It can also generate dat files.
The script takes a folder name as argument.

usage:
    mpwt -f=DIR [-o=DIR] [--patho] [--hf] [--dat] [--md] [--cpu=INT] [-r] [-v] [--clean] [--log=FOLDER] [--ignore-error] [--taxon-file]
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
    --cpu=INT     Number of cpu to use for the multiprocessing (default=1). [default: 1]
    --log=FOLDER     Create PathoLogic log files inside the given folder (use it with --patho).
    --list     List all PGDBs inside the ptools-local folder.
    --ignore-error     Ignore errors (PathoLogic and dat creation) and continue for successful builds.
    --taxon-file     For the use of the taxon_id.tsv file to find the taxon ID.
    -v     Verbose.

"""

import docopt
import logging
import os
import sys

from mpwt import utils
from mpwt.mpwt_workflow import multiprocess_pwt
from multiprocessing import Pool

logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
logger = logging.getLogger('mpwt')
logger.setLevel(logging.CRITICAL)


def run_mpwt():
    """
    Function used with a mpwt call in the terminal.
    """
    args = docopt.docopt(__doc__)

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
    ignore_error = args['--ignore-error']
    taxon_file = args['--taxon-file']
    verbose = args['-v']

    if verbose:
        logger.setLevel(logging.DEBUG)

    if pgdb_list:
        pgdbs = utils.list_pgdb()
        if pgdbs == []:
            logger.critical('No PGDB inside ptools-local.')
        else:
            logger.critical(str(len(pgdbs)) + ' PGDB inside ptools-local:\n' + '\t'.join(pgdbs))
        return

    # Delete PGDB if use of --delete argument.
    # Use a set to remove redudant PGDB.
    if pgdb_to_deletes:
        utils.remove_pgdbs(list(set(pgdb_to_deletes.split(','))), number_cpu)
        return

    if args['--clean']:
        if verbose:
            logger.info('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')

        if input_folder:
            utils.cleaning_input(input_folder, verbose)
            input_pgdb_to_deletes = [species.lower() + 'cyc' for species in os.listdir(input_folder) if not species.startswith('.') and species != 'taxon_id.tsv']
            utils.remove_pgdbs(input_pgdb_to_deletes, number_cpu)
        else:
            utils.cleaning(number_cpu, verbose)
        if not patho_inference and not dat_creation and not move_dat and not output_folder:
            sys.exit()

    multiprocess_pwt(input_folder=input_folder,
                    output_folder=output_folder,
                    patho_inference=patho_inference,
                    patho_hole_filler=patho_hole_filler,
                    dat_creation=dat_creation,
                    dat_extraction=move_dat,
                    size_reduction=size_reduction,
                    number_cpu=number_cpu,
                    patho_log=patho_log,
                    ignore_error=ignore_error,
                    taxon_file=taxon_file,
                    verbose=verbose)


if __name__ == "__main__":
    run_mpwt()
