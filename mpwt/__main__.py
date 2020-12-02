#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
From Genbank/GFF/PF files this script will create Pathway Tools input data, then run Pathway Tools's PathoLogic on them. It can also generate flat files.
The script takes a folder name as argument.

usage:
    mpwt -f=FOLDER [-o=FOLDER] [--patho] [--hf] [--op] [--tp] [--nc] [--flat] [--md] [--mx] [--mo] [--mc] [-p=FLOAT] [--cpu=INT] [-r] [-v] [--clean] [--log=FOLDER] [--ignore-error] [--taxon-file]
    mpwt --flat [-f=FOLDER] [-o=FOLDER] [--md] [--mx] [--mo] [--mc] [--cpu=INT] [-v]
    mpwt -o=FOLDER [--md] [--mx] [--mo] [--mc] [--cpu=INT] [-v]
    mpwt --clean [--cpu=INT] [-v]
    mpwt --delete=STR [--cpu=INT]
    mpwt --list
    mpwt --version
    mpwt topf -f=FOLDER -o=FOLDER [--cpu=INT] [--clean]

options:
    -h --help     Show help.
    -f=FOLDER     Working folder containing sub-folders with Genbank/GFF/PF files.
    -o=FOLDER    Output folder path. Will create a output folder in this folder.
    --patho    Will run an inference of Pathologic on the input files.
    --hf    Use with --patho. Run the Hole Filler using Blast.
    --op    Use with --patho. Run the Operon predictor of Pathway-Tools.
    --tp    Use with --patho. Run the Transport Inference Parser of Pathway-Tools.
    --nc    Use with --patho. Turn off loading of Pubmed entries.
    -p=FLOAT   Use with --patho. Modify PathoLogic pathway prediction score.
    --flat    Will create BioPAX/attribute-value flat files from PGDB.
    --md    Move the dat files into the output folder.
    --mx    Move the metabolic-reactions.xml file into the output folder.
    --mo    Move owl files into the output folder.
    --mc    Move tabular files into the output folder.
    --clean    Clean ptools-local folder, before any other operations.
    --delete=STR    Give a PGDB name and it will delete it (if multiple separe them with a ',', example: ecolicyc,athalianacyc).
    -r    Will delete files in ptools-local and compress results files to reduce results size (use it with -o).
    --cpu=INT     Number of cpu to use for the multiprocessing (default=1). [default: 1]
    --log=FOLDER     Create PathoLogic log files inside the given folder (use it with --patho).
    --list     List all PGDBs inside the ptools-local folder.
    --ignore-error     Ignore errors (PathoLogic and flat-files creation) and continue for successful builds.
    --taxon-file     For the use of the taxon_id.tsv file to find the taxon ID.
    -v     Verbose.
    --version     Version
    topf     Will convert Genbank and/or GFF files into PathoLogic Format file.

"""

import docopt
import logging
import os
import sys
import pkg_resources

from mpwt import utils, to_pathologic
from mpwt.mpwt_workflow import multiprocess_pwt
from multiprocessing import Pool

logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
logger = logging.getLogger(__name__)
logging.getLogger("mpwt").setLevel(logging.CRITICAL)

VERSION = pkg_resources.get_distribution("mpwt").version
LICENSE = """Copyright (C) AuReMe
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
Mpwt is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
"""

def run_mpwt():
    """
    Function used with a mpwt call in the terminal.
    """
    args = docopt.docopt(__doc__)

    input_folder = args['-f']
    output_folder = args['-o']
    patho_inference = args['--patho']
    patho_hole_filler = args['--hf']
    patho_operon_predictor = args['--op']
    patho_transporter_inference = args['--tp']
    no_download_articles = args['--nc']
    flat_creation = args['--flat']
    move_dat = args['--md']
    move_xml = args['--mx']
    move_owl = args['--mo']
    move_col = args['--mc']
    size_reduction = args['-r']
    number_cpu = args['--cpu']
    patho_log = args['--log']
    clean_arg = args['--clean']
    pgdb_to_deletes = args['--delete']
    pgdb_list = args['--list']
    ignore_error = args['--ignore-error']
    taxon_file = args['--taxon-file']
    pathway_score = args['-p']
    verbose = args['-v']
    topf = args['topf']
    version = args['--version']

    if version:
        print("Mpwt v" + VERSION  + "\n" + LICENSE)
        sys.exit()

    if verbose:
        logging.getLogger("mpwt").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    if pgdb_list:
        pgdbs = utils.list_pgdb()
        if pgdbs == []:
            logger.critical('No PGDB inside ptools-local.')
        else:
            logger.critical(str(len(pgdbs)) + ' PGDB inside ptools-local:\n' + '\t'.join(pgdbs))
        sys.exit()

    # Delete PGDB if use of --delete argument.
    # Use a set to remove redudant PGDB.
    if pgdb_to_deletes:
        utils.remove_pgdbs(list(set(pgdb_to_deletes.split(','))), number_cpu)
        return

    if clean_arg:
        if verbose:
            logger.info('~~~~~~~~~~Remove local PGDB~~~~~~~~~~')

        if input_folder:
            utils.cleaning_input(input_folder, verbose)
            input_pgdb_to_deletes = [species.lower() + 'cyc' for species in os.listdir(input_folder) if not species.startswith('.') and species != 'taxon_id.tsv']
            utils.remove_pgdbs(input_pgdb_to_deletes, number_cpu)
        else:
            utils.cleaning(number_cpu, verbose)
        if not patho_inference and not flat_creation and not move_dat and not output_folder:
            sys.exit()

    if topf:
        if input_folder and output_folder:
            to_pathologic.create_pathologic_file(input_folder, output_folder, number_cpu)
        sys.exit()

    multiprocess_pwt(input_folder=input_folder,
                    output_folder=output_folder,
                    patho_inference=patho_inference,
                    patho_hole_filler=patho_hole_filler,
                    patho_operon_predictor=patho_operon_predictor,
                    patho_transporter_inference=patho_transporter_inference,
                    no_download_articles=no_download_articles,
                    flat_creation=flat_creation,
                    dat_extraction=move_dat,
                    xml_extraction=move_xml,
                    owl_extraction=move_owl,
                    col_extraction=move_col,
                    size_reduction=size_reduction,
                    number_cpu=number_cpu,
                    patho_log=patho_log,
                    ignore_error=ignore_error,
                    pathway_score=pathway_score,
                    taxon_file=taxon_file,
                    verbose=verbose)


if __name__ == "__main__":
    run_mpwt()
