#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test mpwt on genbank files containing E. coli genes implied in the TCA cycle and in the Fatty Acid Beta oxydation.
Used mock to avoid Pathway Tools.
"""
import mpwt
import os
import shutil

def fake_run_pwt(species_input_folder_path, patho_hole_filler, patho_operon_predictor, patho_transporter_inference,
            patho_complex_inference, run_flat_creation=None):
    output_folder = 'test_output'
    mpwt_expected = 'mpwt_expected'

    for organism_name in os.listdir(mpwt_expected):
        expected_folder_path = os.path.join(mpwt_expected, organism_name)
        output_organism_folder = os.path.join(output_folder, organism_name)
        if not os.path.exists(output_organism_folder):
            os.mkdir(output_organism_folder)
        for expected_file in os.listdir(expected_folder_path):
            output_organism_file = os.path.join(output_organism_folder, expected_file)
            if not os.path.exists(output_organism_file):
                shutil.copyfile(os.path.join(expected_folder_path, expected_file), output_organism_file)


def fake_run_pwt_flat(species_input_folder_path):
    output_folder = 'test_output'
    mpwt_expected = os.path.join('mpwt_expected')

    for organism_name in os.listdir(mpwt_expected):
        expected_folder_path = os.path.join(mpwt_expected, organism_name)
        output_organism_folder = os.path.join(output_folder, organism_name)
        if not os.path.exists(output_organism_folder):
            os.mkdir(output_organism_folder)
        for expected_file in os.listdir(expected_folder_path):
            output_organism_file = os.path.join(output_organism_folder, expected_file)
            if not os.path.exists(output_organism_file):
                shutil.copyfile(os.path.join(expected_folder_path, expected_file), output_organism_file)


def fake_run_move_pgdb(pgdb_folder_dbname, pgdb_folder_path, output_folder, dat_extraction, size_reduction, xml_extraction, owl_extraction, col_extraction):
    output_folder = 'test_output'
    mpwt_expected = os.path.join('mpwt_expected')

    for organism_name in os.listdir(mpwt_expected):
        expected_folder_path = os.path.join(mpwt_expected, organism_name)
        output_organism_folder = os.path.join(output_folder, organism_name)
        if not os.path.exists(output_organism_folder):
            os.mkdir(output_organism_folder)
        for expected_file in os.listdir(expected_folder_path):
            output_organism_file = os.path.join(output_organism_folder, expected_file)
            if not os.path.exists(output_organism_file):
                shutil.copyfile(os.path.join(expected_folder_path, expected_file), output_organism_file)


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


def test_mpwt_mocked(mocker):
    # Mock call to mpwt by creating expected files.
    mpwt_version = (28, 5)
    mocker.patch("mpwt.utils.get_ptools_version", return_value=mpwt_version)
    mocker.patch("mpwt.mpwt_workflow.run_pwt", wraps=fake_run_pwt)
    mocker.patch("mpwt.mpwt_workflow.run_pwt_flat", wraps=fake_run_pwt_flat)
    mocker.patch("mpwt.mpwt_workflow.run_move_pgdb", wraps=fake_run_move_pgdb)

    mpwt_input = 'test'
    output_folder = 'test_output'
    mpwt.mpwt_workflow.multiprocess_pwt(input_folder=mpwt_input, output_folder=output_folder, patho_inference=True, flat_creation=True, dat_extraction=True)

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I_gff', 'pathways.dat'])
    expected_tca_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_tca_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I', 'pathways.dat'])
    expected_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I_pf', 'pathways.dat'])
    expected_pf_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_pf_fabo_reactions))

    shutil.rmtree('test_output')


def test_mpwt_mocked_already_exisiting_output_folder(mocker):
    # Mock call to mpwt by creating expected files.
    mpwt_version = (28, 5)
    mocker.patch("mpwt.utils.get_ptools_version", return_value=mpwt_version)
    mocker.patch("mpwt.mpwt_workflow.run_pwt", wraps=fake_run_pwt)
    mocker.patch("mpwt.mpwt_workflow.run_pwt_flat", wraps=fake_run_pwt_flat)
    mocker.patch("mpwt.mpwt_workflow.run_move_pgdb", wraps=fake_run_move_pgdb)

    mpwt_input = 'test'
    output_folder = 'test_output'

    # Create output folders before run.
    os.mkdir(os.path.join('test_output'))
    os.mkdir(os.path.join('test_output', 'fatty_acid_beta_oxydation_I_gff'))
    os.mkdir(os.path.join('test_output', 'fatty_acid_beta_oxydation_I'))
    os.mkdir(os.path.join('test_output', 'fatty_acid_beta_oxydation_I_pf'))

    mpwt.mpwt_workflow.multiprocess_pwt(input_folder=mpwt_input, output_folder=output_folder, patho_inference=True, flat_creation=True, dat_extraction=True)

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I_gff', 'pathways.dat'])
    expected_tca_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_tca_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I', 'pathways.dat'])
    expected_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    pathway_fabo_pathname = os.path.join(*['test_output', 'fatty_acid_beta_oxydation_I_pf', 'pathways.dat'])
    expected_pf_fabo_reactions = reaction_extraction(pathway_fabo_pathname)
    assert set(fabo_reactions()).issubset(set(expected_pf_fabo_reactions))

    shutil.rmtree('test_output')