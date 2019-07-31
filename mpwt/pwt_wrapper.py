"""
Wrapping Pathway Tools for PathoLogic and attribute-values dat files creation.
Move results files to an output folder.
"""

import logging
import os
import shutil
import subprocess

logger = logging.getLogger(__name__)


def pwt_error(species_input_folder_path, subprocess_returncode, subprocess_stdout, subprocess_stderr, cmd):
    """
    Print error messages when there is a subprocess error during PathoLogic run.

    Args:
        species_input_folder_path (str): pathname to the species input folder
        subprocess_returncode (int): code returns by subprocess
        subprocess_stdout (str): stdout returns by subprocess
        cmd (str): command used which generates the error
    """
    logger.critical('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')
    species_name = species_input_folder_path.split('/')[-2]
    logger.critical('Error for {0} with PathoLogic subprocess, return code: {1}'.format(species_name, str(subprocess_returncode)))
    if subprocess_stderr:
        logger.critical('An error occurred :' + subprocess_stderr.decode('utf-8'))

    logger.critical('=== Pathway Tools log ===')
    for line in subprocess_stdout:
        if line != '':
            logger.critical('\t' + line)

    # Look for error in pathologic.log.
    if '-patho' in cmd:
        fatal_error_index = None
        with open(species_input_folder_path + '/pathologic.log', 'r') as pathologic_log:
            for index, line in enumerate(pathologic_log):
                if line != '':
                    if 'fatal error' in line or 'Error' in line and not fatal_error_index:
                        fatal_error_index = index
                        logger.critical('=== Error in Pathologic.log ===')
                        logger.critical('\t' + 'Error from the pathologic.log file: {0}'.format(species_input_folder_path + '/pathologic.log'))
                        logger.critical('\t' + line)
                    if fatal_error_index:
                        if index > fatal_error_index:
                            logger.critical('\t' + line)

    logger.critical('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')


def run_pwt(multiprocess_input):
    """
    Create PGDB using files created during 'create_dats_and_lisp' ('organism-params.dat' and 'genetic-elements.dat').
    With verbose run check_output to retrieve the output of subprocess (and show when Pathway Tools has been killed).
    Otherwise send the output to the null device.
    Command used:
    pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho

    Args:
        multiprocess_input (dictionary): contains multiprocess input (mpwt argument: input folder, output folder, ...)
    Returns:
        boolean: True if there is an error during Pathway Tools run
    """
    species_input_folder_path = multiprocess_input['species_input_folder_path']
    patho_hole_filler = multiprocess_input['patho_hole_filler']

    cmd_options = ['-no-web-cel-overview', '-no-cel-overview', '-no-patch-download', '-disable-metadata-saving', '-nologfile']

    cmd_pwt = ['pathway-tools', *cmd_options, '-patho', species_input_folder_path]

    if patho_hole_filler:
        cmd_pwt.append('-hole-filler')

    logger.info(' '.join(cmd_pwt))

    error_status = None
    errors = ['Restart actions (select using :continue):', 'Error']
    patho_lines = []

    try:
        patho_subprocess = subprocess.Popen(cmd_pwt, stdout=subprocess.PIPE, universal_newlines="")
        # Check internal error of Pathway Tools (Error with Genbank).
        for patho_line in iter(patho_subprocess.stdout.readline, ""):
            patho_line = patho_line.decode('utf-8')
            if any(error in patho_line for error in errors):
                logger.info('Error possibly with the genbank file.')
                error_status = True
                patho_subprocess.kill()

            patho_lines.append(patho_line)
            patho_subprocess.poll()
            return_code = patho_subprocess.returncode
            # Check if Pathway Tools has been killed with returncode.
            # Also check if Pathway Tools has finished PathoLogic inference (returncode 0).
            if return_code or return_code == 0:
                if return_code == 0:
                    patho_subprocess.stdout.close()
                    return error_status
                elif return_code != 0:
                    raise subprocess.CalledProcessError(return_code, cmd_pwt)

    except subprocess.CalledProcessError as subprocess_error:
        # Check error with subprocess (when process is killed).
        pwt_error(species_input_folder_path, subprocess_error.returncode, patho_lines, patho_subprocess.stderr, cmd_pwt)
        error_status = True
    patho_subprocess.stdout.close()

    return error_status


def run_pwt_dat(multiprocess_input):
    """
    Create dat file using a lisp script created during 'create_dats_and_lisp'.
    Kill the subprocess when the command reach the Navigator Window opening proposition.
    If this proposition is not closed, the script can't continue.
    Command used:
    pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load

    Args:
        multiprocess_input (dictionary): contains multiprocess input (mpwt argument: input folder, output folder, ...)
    Returns:
        boolean: True if there is an error during lisp script execution
    """
    species_input_folder_path = multiprocess_input['species_input_folder_path']

    lisp_path = species_input_folder_path + 'dat_creation.lisp'
    cmd_options = ['-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_dat = ['pathway-tools', *cmd_options, '-load', lisp_path]

    logger.info(' '.join(cmd_dat))

    error_status = None
    dat_creation_ends = ['Opening Navigator window.', 'No protein-coding genes with sequence data found.  Cannot continue.']
    load_errors = ['Error']
    load_lines = []

    try:
        load_subprocess = subprocess.Popen(cmd_dat, stdout=subprocess.PIPE, universal_newlines="")
        for load_line in iter(load_subprocess.stdout.readline, ""):
            load_line = load_line.decode('utf-8')
            if any(dat_end in load_line for dat_end in dat_creation_ends):
                load_subprocess.stdout.close()
                load_subprocess.kill()
                return
            if any(error in load_line for error in load_errors):
                error_status = True
                load_subprocess.kill()

            load_lines.append(load_line)
            load_subprocess.poll()
            return_code = load_subprocess.returncode
            if return_code:
                raise subprocess.CalledProcessError(return_code, cmd_dat)

    except subprocess.CalledProcessError as subprocess_error:
        # Check error with subprocess (when process is killed).
        pwt_error(species_input_folder_path, subprocess_error.returncode, load_lines, load_subprocess.stderr, cmd_dat)
        error_status = True
    load_subprocess.stdout.close()

    return error_status


def run_move_pgdb(move_data):
    """
    Move the result files inside the shared folder containing the input data.
    pgdb_folder_dbname: ID of the species.
    pgdb_folder_path: path to the PGDB of the species (in ptools-local).

    Args:
        move_data (dictionary): contains multiprocess input (PGDB ID, ptools-local PGDB pathname, ...)
    """
    pgdb_folder_dbname = move_data['pgdb_folders'][0]
    pgdb_folder_path = move_data['pgdb_folders'][1]
    dat_extraction = move_data['dat_extraction']
    output_folder = move_data['output_folder']
    size_reduction = move_data['size_reduction']

    output_species = output_folder + '/' + pgdb_folder_dbname +'/'

    if dat_extraction:
        pgdb_tmp_folder_path = pgdb_folder_path + '/1.0/data'
    else:
        pgdb_tmp_folder_path = pgdb_folder_path

    if size_reduction:
        if dat_extraction:
            for pgdb_file in os.listdir(pgdb_tmp_folder_path):
                pgdb_file_pathname = pgdb_tmp_folder_path + '/' + pgdb_file
                if '.dat' not in pgdb_file:
                    os.remove(pgdb_file_pathname)
        shutil.make_archive(output_folder + '/' + pgdb_folder_dbname, 'zip', pgdb_tmp_folder_path)
        shutil.rmtree(pgdb_folder_path)
    else:
        shutil.copytree(pgdb_tmp_folder_path, output_species)
        if dat_extraction:
            for pgdb_file in os.listdir(output_species):
                if '.dat' not in pgdb_file:
                    os.remove(output_species+'/'+pgdb_file)
