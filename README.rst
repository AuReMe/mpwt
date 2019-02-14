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

mpwt works only on Python 3 and it has been tested on Python 3.6.
It requires some python packages (biopython, docopt and gffutils) and Pathway-Tools. To avoid issues, Pathway-Tools version 22.5 is required.

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
Genbank files must have the same name as the folder in which they are located and also finished with a .gbk.

.. code-block:: text

    Folder_input
    ├── species_1
    │   └── species_1.gbk
    ├── species_2
    │   └── species_2.gbk
    ├── species_3
    │   └── species_3.gbk
    ..

Pathway-Tools will run on each genbank file.
It will create an output folder inside the folder containing all the result files from the PathoLogic inference for each species.
You can also choose another output folder.

Genbank file example:

.. code-block:: text

    LOCUS       scaffold1         XXXXXX bp    DNA     linear   INV DD-MMM-YYYY
    DEFINITION  My species genbank.
    ACCESSION   scaffold1
    VERSION     scaffold1
    KEYWORDS    Key words.
    SOURCE      Source
    ORGANISM  Species name
                Taxonomy; Of; My; Species; With;
                The; Genus.
    FEATURES             Location/Qualifiers
        source          1..XXXXXX
                        /scaffold="scaffold1"
                        /db_xref="taxon:taxonid"
        gene            START..STOP
                        /locus_tag="gene1"
        mRNA            START..STOP
                        /locus_tag="gene1"
        CDS             START..STOP
                        /locus_tag="gene1"
                        /db_xref="InterPro:IPRXXXXXX"
                        /EC_number="X.X.X.X"
                        /translation="AMINOAACIDSSEQUENCE"

Look at the `NCBI GBK format <http://www.insdc.org/files/feature_table.html#7.1.2>`__ for more informations.

GFF file example:

.. code-block:: text

    ##gff-version 3
    ##sequence-region region_id 1 XXXXXX
    region_id	Tools	region	1	XXXXXX	.	strand	.	ID=id0;Dbxref=taxon:id
    region_id	Tools	gene	START	STOP	.	strand	.	ID=geneID
    region_id	Tools	CDS	START	STOP	.	-	strand	ID=cds216;Parent=geneID;ec_number=x.x.x.x

Look at the `NCBI GFF format <https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/>`__ for more informations.

Pathway-Tools does not handle sequence in GFF files. This makes Pathway-Tools run faster compared to a run with a genbank file.
But PGDBs created with this format have missing sequences.

Input files created by mpwt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three input files are created by mpwt. Informations are extracted from the genbank file.
myDBName corresponds to the name of the folder and the genbank file.
taxonid corresponds to the taxonid in the db_xref of the source feature in the genbank.
species_name is extracted from the genbank file.

.. code-block:: text

    organism-params.dat:
    ID  myDBName
    STORAGE FILE
    NCBI-TAXON-ID   taxonid
    NAME    species_name

    genetic-elements.dats:
    NAME    
    ANNOT-FILE  gbk_pathname
    //

    dat_extraction.lisp:
    (in-package :ecocyc)
    (select-organism :org-id 'myDBName)
    (create-flat-files-for-current-kb)

Command Line Example
~~~~~~~~~~~~~~~~~~~~

mpwt is usable as a command line.

.. code:: sh

    mpwt -f path/to/folder/input [-o path/to/folder/output] [--patho] [--dat] [--cpu INT] [-r] [--clean] [--log path/to/folder/log] [-v]

Optional argument are identified by [].

-f input folder as described in Input data.

-o output folder containing PGDB data or dat files (see --dat arguments).

--patho will launch PathoLogic inference on input folder.

--dat will create dat files and only move them inside the output folder.

--cpu the number of cpu used for the multiprocessing.

-r delete files in ptools-local to reduce size of results.

--log folder where log files for PathoLogic inference will be store.

-v print some information about the processing of mpwt.

--delete delete a specific PGDB inside the ptools-local folder.

--clean clean ptools-local folder, before any other operations.

Possible uses of mpwt:

.. code:: sh

    mpwt -f path/to/folder/input --patho

Create PGDBs of studied organisms inside ptools-local.

.. code:: sh

    mpwt -f path/to/folder/input --patho --dat

Create PGDBs of studied organisms inside ptools-local and create dat files.

.. code:: sh

    mpwt -f path/to/folder/input --patho -o path/to/folder/output

Create PGDBs of studied organisms inside ptools-local.
Then extract the files inside the output folder.

.. code:: sh

    mpwt -f path/to/folder/input --patho --dat -o path/to/folder/output

Create PGDBs of studied organisms inside ptools-local and create dat files.
Then extract the dat files inside the output folder.

.. code:: sh

    mpwt --dat -o path/to/folder/output

Create dat files for the PGDB inside ptools-local.
And move them inside the output folder.

Python Example
~~~~~~~~~~~~~~

mpwt can be used in a python script with an import:

.. code:: python

    import mpwt

    folder_input = "path/to/folder/input"
    folder_output = "path/to/folder/output"

    mpwt.multiprocess_pwt(folder_input, folder_output, patho_inference=optional_boolean, dat_extraction=optional_boolean, size_reduction=optional_boolean, number_cpu=int, verbose=optional_boolean)

folder_input: folder containing sub-folders with Genbank file inside.

folder_output: output folder where all the result of Pathway-Tools will be moved. This argument is optional.
If you don't enter an argument, results will be inside the ptools-local folder.

patho_inference: True or nothing. If True, mpwt will launch PathoLogic inference.

dat_extraction: True or nothing. If True, mpwt will create dat files of the PGDBs.

size_reduction: True or nothing. If True, after moving the data to the output folder, mpwt will delete files in ptools-local. This to decrease the size of the results.

number_cpu: int or nothing. Number of cpu to use for the multiprocessing.

verbose: True or nothing. If true, mpwt will be verbose.

Useful functions
~~~~~~~~~~~~~~~~

1. multiprocess_pwt(folder_input, folder_output, patho_inference=optional_boolean, dat_extraction=optional_boolean, size_reduction=optional_boolean, number_cpu=int, verbose=optional_boolean)

Run the multiprocess Pathway-Tools on input folder.

2. cleaning()

Delete all the previous PGDB and the metadata files.

This can also be used with a command line argument:

.. code:: sh

    mpwt --clean

If you use clean and the argument -f input_folder, it will delete input files ('dat_extraction.lisp', 'pathologic.log', 'genetic-elements.dat' and 'organism-params.dat').

.. code:: sh

    mpwt --clean -f input_folder

2. delete_pgdb(pgdb_name)

With this command, it is possible to delete a specified db, where pgdb_name is the name of the PGDB (ending with 'cyc'). It can be multiple pgdbs, to do this, put all the pgdb IDs in a string separated by  a ','.

And as a command line:

.. code:: sh

    mpwt --delete mydbcyc1,mydbcyc2

4. ptools_path()

Return the path to ptools-local.

Errors
~~~~~~

If you encounter errors (and it is highly possible) there is some tips that can help you resolved them.

For error during PathoLogic inference, a log is created where you launch the command.
The log contains the summary of the build and the error for each species.
There is also a pathologic.log in each sub-folders.

If the build passed you have also the possibility to see the result of the inference with the file resume_inference.tsv.
For each species, it contains the number of genes/proteins/reactions/pathways/compounds in the metabolic network.

For others errors, currently nothing is made to help you.
Maybe in the future.
