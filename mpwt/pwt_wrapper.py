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
Wrapping Pathway Tools for PathoLogic and attribute-values flat files creation.
Move results files to an output folder.
"""
import chardet
import logging
import os
import shutil
import signal
import subprocess

logger = logging.getLogger(__name__)


def pwt_check_error(species_input_folder_path, subprocess_stdout, cmd, error_status, subprocess_returncode=None, subprocess_stderr=None):
    """
    Print error messages when there is a subprocess error during PathoLogic run.

    Args:
        species_input_folder_path (str): pathname to the species input folder
        subprocess_stdout (str): stdout returns by subprocess
        cmd (str): command used which generates the error
        error_status (bool): error status of Pathway Tools
        subprocess_returncode (int): code returns by subprocess
        subprocess_stderr (str): stderr returns by subprocess
    Returns:
        boolean: True if there is an error during Pathway Tools run
    """
    if error_status:
        logger.critical('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')
    species_name = os.path.basename(species_input_folder_path)
    if subprocess_returncode:
        logger.critical('Error for {0} with PathoLogic subprocess, return code: {1}'.format(species_name, str(subprocess_returncode)))
    if subprocess_stderr:
        logger.critical('An error occurred :' + subprocess_stderr.decode('utf-8'))

    # Look for error in pathologic.log.
    if '-patho' in cmd:
        pathologic_erros = ['fatal error', 'Error']
        patho_error_status = check_log(species_input_folder_path, 'pathologic.log', error_status, pathologic_erros)
    else:
        patho_error_status = None

    if '-load' in cmd:
        load_errors = ['Error', 'fatal error', 'No protein-coding genes with sequence data found.', 'Cannot continue.']
        flat_error_status = check_log(species_input_folder_path, 'flat_files_creation.log', error_status, load_errors)
    else:
        flat_error_status = None

    error_status = any([error_status, patho_error_status, flat_error_status])

    if error_status:
        logger.critical('=== Pathway Tools log ===')
        for line in subprocess_stdout:
            if line != '':
                logger.critical('\t' + line)
        logger.critical('!!!!!!!!!!!!!!!!!----------------------------------------!!!!!!!!!!!!!!!!!')

    return error_status


def check_log(species_input_folder_path, log_filename, error_status, log_errors):
    """
    Look for error and fatal error in pathologic.log after build.

    Args:
        species_input_folder_path (str): pathname to the species input folder
        log_filename (str): name of the log file
        error_status (bool): True if there is an error during Pathway Tools run
        log_errors (list): list of strings containing possible errors.
    Returns:
        boolean: True if there is an error during Pathway Tools run
    """
    fatal_error_index = None
    log_file_path = os.path.join(species_input_folder_path, log_filename)
    if log_filename == 'flat_files_creation.log':
        encoding = 'utf-8'
    else:
        encoding = 'ascii'
    with open(log_file_path, 'r', encoding=encoding) as log_file:
        for index, line in enumerate(log_file):
            if line != '':
                if not line.startswith(';;;'):
                    if any(error in line for error in log_errors) and not fatal_error_index:
                        fatal_error_index = index
                        logger.critical('=== Error in {0} for {1}==='.format(log_filename, species_input_folder_path))
                        logger.critical('\t' + 'Error from the {0} file: {1}'.format(log_filename, log_file_path))
                        logger.critical('\t' + line)
                        error_status = True
                    if fatal_error_index:
                        if index > fatal_error_index:
                            logger.critical('\t' + line)

    return error_status


def run_pwt(species_input_folder_path, patho_hole_filler, patho_operon_predictor, patho_transporter_inference):
    """
    Create PGDB using files created during 'create_flats_and_lisp' ('organism-params.dat' and 'genetic-elements.dat').
    With verbose run check_output to retrieve the output of subprocess (and show when Pathway Tools has been killed).
    Otherwise send the output to the null device.
    Command used:
    pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho

    Args:
        species_input_folder_path (str): path to input folder
        patho_hole_filler (bool): boolean to use or not PathoLogic Hole Filler
        patho_operon_predictor (bool): boolean to use or not PathoLogic Operon Predictor
        patho_transporter_inference (bool): boolean to use or not PathoLogic Transport Inference Parser
    Returns:
        boolean: True if there is an error during Pathway Tools run
    """
    cmd_options = ['-no-web-cel-overview', '-no-cel-overview', '-no-patch-download', '-disable-metadata-saving', '-nologfile']

    cmd_pwt = ['pathway-tools', *cmd_options, '-patho', species_input_folder_path]

    if patho_hole_filler:
        cmd_pwt.append('-hole-filler')

    if patho_operon_predictor:
        cmd_pwt.append('-operon-predictor')

    if patho_transporter_inference:
        cmd_pwt.append('-tip')

    species_name = os.path.basename(species_input_folder_path)
    logger.info('|PathoLogic|{}| '.format(species_name) + ' '.join(cmd_pwt))

    error_status = None
    # Errors are either a fatal error or opening the lisp listener.
    errors = ['fatal error', '[Current process: Initial Lisp Listener]']
    patho_lines = []

    # Name of the file containing the log from Pathway Tools terminal.
    pwt_log = os.path.join(species_input_folder_path, 'pwt_terminal.log')

    try:
        # Launch Pathway Tools PathoLogic.
        # Use start_new_session to group process ID to kill this process and its childs (with os.killpg).
        patho_subprocess = subprocess.Popen(cmd_pwt, stdout=subprocess.PIPE, start_new_session=True, universal_newlines="")
        with open(pwt_log, 'w', encoding='utf-8') as pwt_writer:
            for patho_line in iter(patho_subprocess.stdout.readline, b''):
                encoding = chardet.detect(patho_line)['encoding']
                patho_line = patho_line.decode(encoding, errors='replace')
                pwt_writer.write(patho_line)

                # An error occured, kill Pathway Tools.
                if any(error in patho_line for error in errors):
                    logger.info('Error possibly with the genbank file.')
                    error_status = True
                    patho_subprocess.kill()
                    os.killpg(os.getpgid(patho_subprocess.pid), signal.SIGKILL)

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
        error_status = True
        error_status = pwt_check_error(species_input_folder_path, patho_lines, cmd_pwt, error_status, subprocess_error.returncode, patho_subprocess.stderr)

    # Check pathologic.log for error.
    error_status = pwt_check_error(species_input_folder_path, patho_lines, cmd_pwt, error_status)

    patho_subprocess.stdout.close()

    return error_status


def run_pwt_flat(species_input_folder_path):
    """
    Create flat files using a lisp script created during 'create_flats_and_lisp'.
    Kill the subprocess when the command reach the Navigator Window opening proposition.
    If this proposition is not closed, the script can't continue.
    Command used:
    pathway-tools -no-patch-download -disable-metadata-saving -nologfile -load

    Args:
        species_input_folder_path (str): path to input folder
    Returns:
        boolean: True if there is an error during lisp script execution
    """
    lisp_path = os.path.join(species_input_folder_path, 'flat_files_creation.lisp')
    cmd_options = ['-no-patch-download', '-disable-metadata-saving', '-nologfile']
    cmd_flat = ['pathway-tools', *cmd_options, '-load', lisp_path]

    species_name = os.path.basename(species_input_folder_path)
    logger.info('|Flat files creation|{}| '.format(species_name) + ' '.join(cmd_flat))

    error_status = None
    flat_creation_ends = ['Opening Navigator window.']
    # Errors are either a fatal error or opening the lisp listener.
    load_errors = ['fatal error', '[Current process: Initial Lisp Listener]']
    load_lines = []

    # Name of the file containing the log from Pathway Tools terminal.
    flat_log = os.path.join(species_input_folder_path, 'flat_files_creation.log')

    try:
        # Launch Pathway Tools lisp command.
        # Use start_new_session to group process ID to kill this process and its childs (with os.killpg).
        load_subprocess = subprocess.Popen(cmd_flat, stdout=subprocess.PIPE, start_new_session=True, universal_newlines="")
        with open(flat_log, 'w', encoding='utf-8') as flat_file_writer:
            for load_line in iter(load_subprocess.stdout.readline, b''):
                encoding = chardet.detect(load_line)['encoding']
                load_line = load_line.decode(encoding, errors='replace')
                flat_file_writer.write(load_line)

                # Lisp command has finished, kill Pathway Toosl trying to open navigator.
                if any(flat_end in load_line for flat_end in flat_creation_ends):
                    load_subprocess.stdout.close()
                    load_subprocess.kill()
                    os.killpg(os.getpgid(load_subprocess.pid), signal.SIGKILL)
                    return

                # There is an error, kill lisp command and return the error.
                if any(error in load_line for error in load_errors):
                    if not load_line.startswith(';;;'):
                        error_status = True
                        load_lines.append(load_line)
                        load_subprocess.kill()
                        os.killpg(os.getpgid(load_subprocess.pid), signal.SIGKILL)

                load_lines.append(load_line)
                load_subprocess.poll()
                return_code = load_subprocess.returncode
                if return_code:
                    raise subprocess.CalledProcessError(return_code, cmd_flat)

        # Check for error.
        error_status = pwt_check_error(species_input_folder_path, load_lines, cmd_flat, error_status)

    except subprocess.CalledProcessError as subprocess_error:
        # Check error with subprocess (when process is killed).
        error_status = pwt_check_error(species_input_folder_path, load_lines, cmd_flat, error_status, subprocess_error.returncode, load_subprocess.stderr)

    load_subprocess.stdout.close()

    return error_status


def run_move_pgdb(pgdb_folder_dbname, pgdb_folder_path, output_folder, dat_extraction, size_reduction, xml_extraction, owl_extraction, col_extraction):
    """
    Move the result files inside the shared folder containing the input data.
    pgdb_folder_dbname: ID of the species.
    pgdb_folder_path: path to the PGDB of the species (in ptools-local).

    Args:
        pgdb_folder_dbname (str): species ID
        pgdb_folder_path (str): path to species PGDB folder
        output_folder (str): path to output folder
        dat_extraction (bool): to extract or not the attribute-values files (.dat files)
        size_reduction (bool): to compress or not the results
        xml_extraction (bool): to extract or not the metabolic-reactions.xml'
        owl_extraction (bool): to extract or not the owl files
        col_extraction (bool): to extract or not the tabular files (.col files)
    """
    logger.info('|Moving output files|{}| '.format(pgdb_folder_dbname))

    output_species = os.path.join(output_folder, pgdb_folder_dbname)

    keep_extensions = []
    if dat_extraction:
        keep_extensions.append('.dat')
    if xml_extraction:
        keep_extensions.append('.xml')
    if owl_extraction:
        keep_extensions.append('.owl')
    if col_extraction:
        keep_extensions.append('.col')

    if xml_extraction and len(keep_extensions) == 1:
        xml_file = os.path.join(*[pgdb_folder_path, '1.0', 'data', 'metabolic-reactions.xml'])
        if not os.path.exists(xml_file):
            logger.critical('Missing metabolic-reactions.xml for ' + pgdb_folder_dbname + '. To create it use mpwt --dat.')
            return
        shutil.copyfile(xml_file, output_species+'.xml')
        if size_reduction:
            shutil.rmtree(pgdb_folder_path)
        return

    if owl_extraction and len(keep_extensions) == 1:
        owl_level_2 = os.path.join(*[pgdb_folder_path, '1.0', 'data', 'biopax-level2.owl'])
        if not os.path.exists(owl_level_2):
            logger.critical('Missing biopax-level2.owl for ' + pgdb_folder_dbname + '. To create it use mpwt --dat.')
            return
        owl_level_3 = os.path.join(*[pgdb_folder_path, '1.0', 'data', 'biopax-level3.owl'])
        if not os.path.exists(owl_level_3):
            logger.critical('Missing biopax-level3.owl for ' + pgdb_folder_dbname + '. To create it use mpwt --dat.')
            return
        shutil.copyfile(owl_level_2, output_species+'-level2.owl')
        shutil.copyfile(owl_level_3, output_species+'-level3.owl')
        if size_reduction:
            shutil.rmtree(pgdb_folder_path)
        return

    if len(keep_extensions) > 0:
        pgdb_tmp_folder_path = os.path.join(*[pgdb_folder_path, '1.0', 'data'])
    else:
        pgdb_tmp_folder_path = pgdb_folder_path

    if not os.path.exists(pgdb_tmp_folder_path):
        logger.critical('Missing ' + pgdb_tmp_folder_path + ' folder.')
        return True

    # If size_reduction, mpwt will create a compressed version of the PGDB in output folder.
    # It will also delete the PGDB folder in ptools-local.
    if size_reduction:
        if len(keep_extensions) > 0:
            for pgdb_file in os.listdir(pgdb_tmp_folder_path):
                pgdb_file_pathname = os.path.join(pgdb_tmp_folder_path, pgdb_file)
                pgdb_file_extension = os.path.splitext(pgdb_file)[1]
                if pgdb_file_extension not in keep_extensions:
                    if os.path.isfile(pgdb_file_pathname):
                        os.remove(pgdb_file_pathname)
                    elif os.path.isdir(pgdb_file_pathname):
                        shutil.rmtree(pgdb_file_pathname)
        zip_input_path = os.path.join(output_folder, pgdb_folder_dbname)
        shutil.make_archive(zip_input_path, 'zip', pgdb_tmp_folder_path)
        shutil.rmtree(pgdb_folder_path)

    else:
        shutil.copytree(pgdb_tmp_folder_path, output_species)
        if len(keep_extensions) > 0:
            for pgdb_file in os.listdir(output_species):
                pgdb_file_pathname = os.path.join(output_species, pgdb_file)
                pgdb_file_extension = os.path.splitext(pgdb_file)[1]
                if pgdb_file_extension not in keep_extensions:
                    if os.path.isfile(pgdb_file_pathname):
                        os.remove(pgdb_file_pathname)
                    elif os.path.isdir(pgdb_file_pathname):
                        shutil.rmtree(pgdb_file_pathname)

    return False