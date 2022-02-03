#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2018-2022 Arnaud Belcour - Inria Dyliss
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

"""
Description:
From Genbank/GFF/PF files this script will create Pathway Tools input data, then run Pathway Tools's PathoLogic on them. It can also generate flat files.
The script takes a folder name as argument.

usage:
    mpwt -f=FOLDER [-o=FOLDER] [--patho] [--hf] [--op] [--tp] [--nc] [--flat] [--md] [--mx] [--mo] [--mc] [-p=FLOAT] [--cpu=INT] [-r] [-v] [--clean] [--log=FOLDER] [--taxon-file]
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
    --taxon-file     For the use of the taxon_id.tsv file to find the taxon ID.
    --permission     Choose permission access to PGDB in ptools-local and output files, either 'all' or 'group' (by default it is user).
    -v     Verbose.
    --version     Version
    topf     Will convert Genbank and/or GFF files into PathoLogic Format file.

"""

import argparse
import logging
import os
import sys
import pkg_resources

from mpwt import utils, to_pathologic
from mpwt.mpwt_workflow import multiprocess_pwt

logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
logger = logging.getLogger(__name__)
logging.getLogger('mpwt').setLevel(logging.CRITICAL)

VERSION = pkg_resources.get_distribution('mpwt').version
LICENSE = """Copyright (C) 2018-2022 Arnaud Belcour - Inria Dyliss\n
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.\n

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.\n

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>\n
"""

def run_mpwt():
    """
    Function used with a mpwt call in the terminal.
    """
    parser = argparse.ArgumentParser(
        'mpwt',
        description='For specific help on each subcommand use: mpwt --help',
        allow_abbrev=False,
    )

    parser.add_argument(
        '-f',
        dest='input',
        required=False,
        help='Working folder containing sub-folders with Genbank/GFF/PF files.',
        metavar='INPUT_DIR')

    parser.add_argument(
        '-o',
        dest='output',
        required=False,
        help='Output folder path. Will create a output folder in this folder.',
        metavar='OUPUT_DIR')

    parser.add_argument(
        '--patho',
        dest='patho',
        help='Will run an inference of Pathologic on the input files.',
        required=False,
        action='store_true',
        default=False,
    )

    parser.add_argument(
        '--hf',
        dest='hf',
        help='Use with --patho. Run the Hole Filler using Blast.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--op',
        dest='op',
        help='Use with --patho. Run the Operon predictor of Pathway-Tools.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--tp',
        dest='tp',
        help='Use with --patho. Run the Transport Inference Parser of Pathway-Tools.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--nc',
        dest='nc',
        help='Use with --patho. Turn off loading of Pubmed entries.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '-p',
        dest='p',
        help='Use with --patho. Modify PathoLogic pathway prediction score. Must be a float between 0 and 1.',
        required=False,
    )

    parser.add_argument(
        '--flat',
        dest='flat',
        help='Will create BioPAX/attribute-value flat files from PGDB.',
        required=False,
        action='store_true',
        default=False,
    )

    parser.add_argument(
        '--md',
        dest='md',
        help='Move the dat files into the output folder.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--mx',
        dest='mx',
        help='Move the metabolic-reactions.xml file into the output folder.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--mo',
        dest='mo',
        help='Move owl files into the output folder.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--mc',
        dest='mc',
        help='Move tabular files into the output folder.',
        required=False,
        action='store_true',
        default=False,
    )

    parser.add_argument(
        '--clean',
        dest='clean',
        help='Clean ptools-local folder, before any other operations.',
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--delete',
        dest='delete',
        help="Give a PGDB name and mpwt will delete it (if multiple separe them with a ',', example: ecolicyc,athalianacyc).",
        required=False,
    )
    parser.add_argument(
        '-r',
        dest='r',
        help="Will delete files in ptools-local and compress results files to reduce results size (use it with -o).",
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--log',
        dest='log',
        help="Create PathoLogic log files inside the given folder (use it with --patho).",
        required=False,
    )
    parser.add_argument(
        '--list',
        dest='list',
        help="List all PGDBs inside the ptools-local folder.",
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--taxon-file',
        dest='taxon_file',
        help="For the use of the taxon_id.tsv file to find the taxon ID.",
        required=False,
    )
    parser.add_argument(
        '-v',
        dest='verbose',
        help="Verbose.",
        required=False,
        action='store_true',
        default=False,
    )
    parser.add_argument(
        'topf',
        help="Will convert Genbank and/or GFF files into PathoLogic Format file.",
        nargs='?',
    )

    parser.add_argument(
        '--version',
        dest='version',
        action='version',
        default=False,
        version='%(prog)s ' + VERSION + '\n' + LICENSE)

    parser.add_argument(
        '--cpu',
        help='Number of cpu to use for the multiprocessing (default=1). [default: 1]',
        required=False,
        type=int,
        default=1)
    parser.add_argument(
        '--permission',
        dest='permission',
        help="Choose permission access to PGDB in ptools-local and output files, either 'all' or 'group' (by default it is user).",
        required=False,
    )

    args = parser.parse_args()

    input_folder = args.input
    output_folder = args.output
    patho_inference = args.patho
    patho_hole_filler = args.hf
    patho_operon_predictor = args.op
    patho_transporter_inference = args.tp
    no_download_articles = args.nc
    flat_creation = args.flat
    move_dat = args.md
    move_xml = args.mx
    move_owl = args.mo
    move_col = args.mc
    size_reduction = args.r
    number_cpu = args.cpu
    patho_log = args.log
    clean_arg = args.clean
    pgdb_to_deletes = args.delete
    pgdb_list = args.list
    taxon_file = args.taxon_file
    pathway_score = args.p
    verbose = args.verbose
    topf = args.topf
    version = args.version
    permission = args.permission

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if version:
        print('Mpwt v' + VERSION  + '\n' + LICENSE)
        sys.exit()

    if verbose:
        logging.getLogger('mpwt').setLevel(logging.DEBUG)
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
        sys.exit()

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

    if topf is not None:
        if topf == 'topf':
            if input_folder and output_folder:
                to_pathologic.create_pathologic_file(input_folder, output_folder, number_cpu)
                sys.exit()
            else:
                sys.exit('topf argument needs input_folder (-f) and output_folder options (-o).')
        else:
            sys.exit(f'Wrong positional argument passed: {topf}, only "topf" is expected as a postional argument.')

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
                    pathway_score=pathway_score,
                    taxon_file=taxon_file,
                    verbose=verbose,
                    permission=permission)


if __name__ == '__main__':
    run_mpwt()
