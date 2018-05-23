.. image:: https://img.shields.io/pypi/v/mpwt.svg
	:target: https://pypi.python.org/pypi/mpwt

Pathway-tools multiprocessing script
====================================

mpwt is a python package for running Pathway-Tools on multiple genomes using multiprocessing.

There is no guarantee that this script will work, it is a Work In Progress in early state.

.. contents:: Table of contents
   :backlinks: top
   :local:

Installation
------------

Requirements
~~~~~~~~~~~~

You must have an environment where Pathway-Tools is installed. Pathway-Tools can be obtained `here <http://bioinformatics.ai.sri.com/ptools/>`__.
For some versions you need to have Blast installed on you system, for further informations look at `this page <http://bioinformatics.ai.sri.com/ptools/installation-guide/released/blast.html>`__.

If your OS doesn't support Pathway-Tools, you can use a docker. If it's your case, look at `Pathway-Tools Multiprocessing Docker <https://github.com/ArnaudBelcour/pathway-tools-multiprocessing-docker>`__.
It is a dockerfile that will create container with Pathway-Tools, its dependancies and this package. You just need to give a Pathway-Tools installer as input.

Using pip
~~~~~~~~~

.. code:: sh

	pip install mpwt

Use
---

Input data
~~~~~~~~~~

The script takes a folder containing sub-folders as input. Each sub-folder contains a genbank file.

.. code-block:: text

    Folder_input
    ├── Folder_for_species_1
    │   └── Genbank_species_1
    ├── Folder_for_species_2
    │   └── Genbank_species_2
    ├── Folder_for_species_3
    │   └── Genbank_species_3
    │

Pathway-Tools will run on each genbank file. It will create an output folder containing all the result files from the PathoLogic inference for each species.

Example
~~~~~~~

.. code:: python

    import mpwt

    folder_input = "path/to/folder/input"
    folder_output = "path/to/folder/output"

    mpwt.multiprocess_pwt(folder_input, folder_output)

Useful functions
~~~~~~~~~~~~~~~~

1. multiprocess_pwt(folder_input, folder_output)

folder_input: folder containing sub-folders with Genbank file inside.

folder_output: output folder where all the result of Pathway-Tools will be moved. This argument is optional
If you don't enter an argument, results will be stored in a folder named output inside the sub-folders containg Genbank file.

2. cleaning()

Delete all the previous PGDB and the metadata files.

Errors
~~~~~~

If you encounter errors (and it is highly possible) there is some tips that can help you resolved them.

For error during PathoLogic inference, a log is created where you launch the command.
The log contains the summary of the build and the error for each species.
There is also a pathologic.log in each sub-folders.

For others errors, currently nothing is made to help you.
Maybe in the future.
