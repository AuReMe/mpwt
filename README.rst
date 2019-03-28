.. image:: https://img.shields.io/pypi/v/mpwt.svg
	:target: https://pypi.python.org/pypi/mpwt

Pathway-Tools multiprocessing script
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

If your OS doesn't support Pathway-Tools, you can use a docker. If it's your case, look at `Pathway-Tools Multiprocessing Docker <https://github.com/ArnaudBelcour/mpwt-docker>`__.
It is a dockerfile that will create a container with Pathway-Tools, its dependencies and this package. You just need to give a Pathway-Tools installer as input.

You can also look at `Pathway-Tools Multiprocessing Singularity <https://github.com/ArnaudBelcour/mpwt-singularity>`__.
More manipulations are required compared to Docker but with this you can create a Singularity image.

Using pip
~~~~~~~~~

.. code:: sh

	pip install mpwt

Use
---

Input data
~~~~~~~~~~

The script takes a folder containing sub-folders as input. Each sub-folder contains a Genbank/GFF file.
Genbank files must have the same name as the folder in which they are located and also finished with a .gbk or a .gff.

.. code-block:: text

    Folder_input
    ├── species_1
    │   └── species_1.gbk
    ├── species_2
    │   └── species_2.gff
    │   └── species_2.fasta
    ├── species_3
    │   └── species_3.gbk
    ..

Pathway-Tools will run on each Genbank/GFF file. It will create the results in the ptools-local folder but you can also choose an output folder.

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
You can also look at the `example <http://bioinformatics.ai.sri.com/ptools/sample.gbff>`__ provided on Pathway-Tools site.

GFF file example:

.. code-block:: text

    ##gff-version 3
    ##sequence-region scaffold_1 1 XXXXXX
    scaffold_1	RefSeq	region	1	XXXXXXX	.	+	.	ID=region_id;Dbxref=taxon:XXXXXX
    scaffold_1	RefSeq	gene	START	STOP	.	-	.	ID=gene_id
    scaffold_1	RefSeq	CDS	START	STOP	.	-	0	ID=cds_id;Parent=gene_id

**Warning**: it seems that metabolic networks from GFF file have less reactions/pathways/compounds than metabolic networks from Genbank file.
Lack of some annotations (EC, GO) can be the reason explaining these differences.

Look at the `NCBI GFF format <https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/>`__ for more informations.

You have to provide a nucleotide sequence file associated with the GFF file containing the chromosome/scaffold/contig sequence.

.. code-block:: text

    >scaffold_1
    ATGATGCTGATACTGACTTAGCAT


Input files created by mpwt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three input files are created by mpwt. Informations are extracted from the Genbank/GFF file.
myDBName corresponds to the name of the folder and the Genbank/GFF file.
taxonid corresponds to the taxonid in the db_xref of the source feature in the Genbank/GFF.
species_name is extracted from the Genbank/GFF file.

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

    dat_creation.lisp:
    (in-package :ecocyc)
    (select-organism :org-id 'myDBName)
    (let ((*progress-noter-enabled?* NIL))
            (create-flat-files-for-current-kb))

Command Line and Python arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mpwt can be used as a command line.

.. code:: sh

    mpwt -f path/to/folder/input [-o path/to/folder/output] [--patho] [--hf] [--dat] [--md] [--cpu INT] [-r] [--clean] [--log path/to/folder/log] [-v]

Optional argument are identified by [].

mpwt can be used in a python script with an import:

.. code:: python

    import mpwt

    folder_input = "path/to/folder/input"
    folder_output = "path/to/folder/output"

    mpwt.multiprocess_pwt(folder_input,
    			  folder_output,
			  patho_inference=optional_boolean,
			  patho_hole_filler=optional_boolean,
			  dat_creation=optional_boolean,
			  dat_extraction=optional_boolean,
			  size_reduction=optional_boolean,
			  number_cpu=int,
			  patho_log=optional_folder_pathname,
			  verbose=optional_boolean)

Command line argument / Python argument: description

-f / folder_input(string: folder pathname): input folder as described in Input data.

-o / folder_output(string: folder pathname): output folder containing PGDB data or dat files (see --dat arguments).

--patho / patho_inference(boolean): will launch PathoLogic inference on input folder.

--hf /patho_hole_filler(boolean): (to use with --patho) will launch PathoLogic Hole Filler with Blast.

--dat / dat_creation(boolean): will create BioPAX/attribute-value dat files.

--md /dat_extraction(boolean): will move only the dat files inside the output folder.

--cpu / number_cpu(int): the number of cpu used for the multiprocessing.

-r / dat_extraction(boolean): delete files in ptools-local to reduce size of results.

--log / patho_log(string: folder pathname): folder where log files for PathoLogic inference will be store.

-v / verbose(boolean): print some information about the processing of mpwt.

--delete / mpwt.remove_pgdbs()(string: pgdb name): delete a specific PGDB inside the ptools-local folder.

--clean / mpwt.cleaning(): clean ptools-local folder, before any other operations.


Examples
~~~~~~~~

Possible uses of mpwt:

.. code:: sh

    mpwt -f path/to/folder/input --patho

.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
			  patho_inference=True)

Create PGDBs of studied organisms inside ptools-local.

.. code:: sh

    mpwt -f path/to/folder/input --patho --hf --log path/to/folder/log

.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
			  patho_inference=True,
			  patho_hole_filler=True,
			  patho_log='path/to/folder/log')

Create PGDBs of studied organisms inside ptools-local with the Hole-Filler.

.. code:: sh

    mpwt -f path/to/folder/input --patho --dat

.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
			  patho_inference=True,
                          dat_creation=True)

Create PGDBs of studied organisms inside ptools-local and create dat files.

.. code:: sh

    mpwt -f path/to/folder/input --patho -o path/to/folder/output

.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
                          folder_output='path/to/folder/output',
			  patho_inference=True)

Create PGDBs of studied organisms inside ptools-local.
Then move the files to the output folder.

.. code:: sh

    mpwt -f path/to/folder/input --patho --dat -o path/to/folder/output --md


.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
                          folder_output='path/to/folder/output',
			  patho_inference=True,
                          dat_creation=True,
			  dat_extraction=True)

Create PGDBs of studied organisms inside ptools-local and create dat files.
Then move the dat files to the output folder.

.. code:: sh

    mpwt --dat -o path/to/folder/output --md

.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(folder_output='path/to/folder/output',
                          dat_creation=True,
			  dat_extraction=True)

Create dat files for the PGDB inside ptools-local.
And move them to the output folder.

.. code:: sh

    mpwt -o path/to/folder/output

.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(folder_output='path/to/folder/output')

Move PGDB from ptools-local to the output folder.

.. code:: sh

    mpwt -o path/to/folder/output --md

.. code:: python

    import mpwt
    mpwt.multiprocess_pwt(folder_output='path/to/folder/output',
			  dat_extraction=True)

Move dat files from ptools-local to the output folder.


Useful functions
~~~~~~~~~~~~~~~~

1. multiprocess_pwt(folder_input, folder_output, patho_inference=optional_boolean, dat_creation=optional_boolean, dat_extraction=optional_boolean, size_reduction=optional_boolean, number_cpu=int, verbose=optional_boolean)

Run the multiprocess Pathway-Tools on input folder.

2. cleaning()

Delete all the previous PGDB and the metadata files.

This can also be used with a command line argument:

.. code:: sh

    mpwt --clean

If you use clean and the argument -f input_folder, it will delete input files ('dat_creation.lisp', 'pathologic.log', 'genetic-elements.dat' and 'organism-params.dat').

.. code:: sh

    mpwt --clean -f input_folder

2. remove_pgdbs(pgdb_name)

With this command, it is possible to delete a specified db, where pgdb_name is the name of the PGDB (ending with 'cyc'). It can be multiple pgdbs, to do this, put all the pgdb IDs in a string separated by  a ','.

And as a command line:

.. code:: sh

    mpwt --delete mydbcyc1,mydbcyc2

4. ptools_path()

Return the path to ptools-local.

5. list_pgdb()

Return a list containing all the PGDBs inside ptools-local folder. Can be used as a command with:

.. code:: sh

    mpwt --list

Errors
~~~~~~

If you encounter errors (and it is highly possible) there is some tips that can help you resolved them.

For error during PathoLogic inference, you can use the log arguments.
The log contains the summary of the build and the error for each species.
There is also a pathologic.log in each sub-folders.

If the build passed you have also the possibility to see the result of the inference with the file resume_inference.tsv.
For each species, it contains the number of genes/proteins/reactions/pathways/compounds in the metabolic network.

If Pathway-Tools crashed, mpwt can print some useful information in verbose mode.

Output
~~~~~~

If you did not use the output argument, results (PGDB with/without BioPAX/dat files) will be inside your ptools-local folder ready to be used with Pathway-Tools.
Have in mind that mpwt does not create the cellular overview and does not used the hole-filler. So if you want these results you should run them after.

If you used the output argument, there is two potential outputs depending on the use of the option --md/dat_extraction:

1. without this option, you will have a complete PGDB folder inside your results, for example:

.. code-block:: text

    Folder_output
    ├── species_1
    │   └── default-version
    │   └── 1.0
    │       └── data
    │           └── contains BioPAX/dat files if you used the --dat/dat_creation option.
    │       └── input
    │           └── species_1.gbk
    │           └── genetic-elements.dat
    │           └── organism-init.dat
    │           └── organism.dat
    │       └── kb
    │           └── species_1.ocelot
    │       └── reports
    │           └── contains Pathway-Tools reports.
    ├── species_2
    ..
    ├── species_3
    ..

2. with this option, you will only have the dat files, for example:

.. code-block:: text

    Folder_output
    ├── species_1
    │   └── classes.dat
    │   └── compounds.dat
    │   └── dnabindsites.dat
    │   └── enzrxns.dat
    │   └── genes.dat
    │   └── pathways.dat
    │   └── promoters.dat
    │   └── protein-features.dat
    │   └── proteins.dat
    │   └── protligandcplxes.dat
    │   └── pubs.dat
    │   └── reactions.dat
    │   └── regulation.dat
    │   └── regulons.dat
    │   └── rnas.dat
    │   └── species.dat
    │   └── terminators.dat
    │   └── transunits.dat
    │   └── ..
    ├── species_2
    ..
    ├── species_3
    ..

Release Notes
~~~~~~~~~~~~~

Changes between version are listed on the `release page <https://github.com/AuReMe/mpwt/releases>`__.
