.. image:: https://img.shields.io/pypi/v/mpwt.svg
	:target: https://pypi.python.org/pypi/mpwt

mpwt: Pathway Tools multiprocessing wrapper
===========================================

mpwt is a python package for running Pathway Tools on multiple genomes using multiprocessing.

There is no guarantee that this script will work, it is a Work In Progress in early state.

.. contents:: Table of contents
   :backlinks: top
   :local:

Installation
------------

Requirements
~~~~~~~~~~~~

mpwt works only on **Python 3** and it has been tested on Python 3.6.
It requires some python packages (`biopython <https://github.com/biopython/biopython>`__, `docopt <https://github.com/docopt/docopt>`__ and `gffutils <https://github.com/daler/gffutils>`__) and **Pathway Tools**. For the multiprocessing, mpwt uses the `multiprocessing library of Python 3 <https://docs.python.org/3/library/multiprocessing.html>`__.

You must have an environment where Pathway Tools is installed. Pathway Tools can be obtained `here <http://bioinformatics.ai.sri.com/ptools/>`__.

Pathway Tools needs **Blast**, so it must be install on your system. Depending on your system, Pathway Tools needs a file named **.ncbirc** to locate Blast, for more informations look at `this page <http://bioinformatics.ai.sri.com/ptools/installation-guide/released/blast.html>`__.

If your OS doesn't support Pathway Tools, you can use a docker. If it's your case, look at `Pathway Tools Multiprocessing Docker <https://github.com/ArnaudBelcour/mpwt-docker>`__.
It is a dockerfile that will create a container with Pathway Tools, its dependencies and this package. You just need to give a Pathway Tools installer as input.

You can also look at `Pathway Tools Multiprocessing Singularity <https://github.com/ArnaudBelcour/mpwt-singularity>`__.
More manipulations are required compared to Docker but with this you can create a Singularity image.

Using pip
~~~~~~~~~

.. code:: sh

	pip install mpwt

Use
---

Input data
~~~~~~~~~~

The script takes a folder containing sub-folders as input. Each sub-folder contains a Genbank/GFF file or multiple PathoLogic Format (PF) files.

.. code-block:: text

    Folder_input
    ├── species_1
    │   └── species_1.gbk
    ├── species_2
    │   └── species_2.gff
    │   └── species_2.fasta
    ├── species_3
    │   └── species_3.gbk
    ├── species_4
    │   └── scaffold_1.pf
    │   └── scaffold_1.fasta
    │   └── scaffold_2.pf
    │   └── scaffold_2.fasta
    taxon_id.tsv
    ..

Genbank files must have the same name as the folder in which they are located and also finished with a .gbk or a .gff.

For PF files, there is one file for each scaffold/contig and one corresponding fasta file.

Pathway Tools will run on each Genbank/GFF/PF files. It will create the results in the ptools-local folder but you can also choose an output folder.

Genbank
+++++++

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
                        /go_component="GO:XXXXXXX"
                        /EC_number="X.X.X.X"
                        /translation="AMINOAACIDSSEQUENCE"

Look at the `NCBI GBK format <http://www.insdc.org/files/feature_table.html#7.1.2>`__ for more informations.
You can also look at the `example <http://bioinformatics.ai.sri.com/ptools/sample.gbff>`__ provided on Pathway Tools site.

GFF
+++

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

PathoLogic Format
+++++++++++++++++

PF file example:

.. code-block:: text

    ;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; scaffold_1
    ;;;;;;;;;;;;;;;;;;;;;;;;;
    ID	gene_id
    NAME	gene_id
    STARTBASE	START
    ENDBASE	STOP
    FUNCTION	ORF
    PRODUCT-TYPE	P
    PRODUCT-ID	prot gene_id
    EC	X.X.X.X
    DBLINK	GO:XXXXXXX
    INTRON	START1-STOP1
    //

Look at the `Pathologic format <http://bioinformatics.ai.sri.com/ptools/tpal.pf>`__ for more informations.

You have to provide one nucleotide sequence for each pathologic containing one scaffold/contig.

.. code-block:: text

    >scaffold_1
    ATGATGCTGATACTGACTTAGCAT

Also to add the taxon ID we need the **taxon_id.tsv** (a tsv file with two values: the name of the folder containing the PF files and the taxon ID corresponding).

+------------+------------+
|species     |taxon_id    |
+============+============+
|species_4   |4           |
+------------+------------+

If you don't have taxon ID in your Genbank or GFF file, you can add one in this file for the corresponding species.

You can also add more informations for the genetic elements like **circularity of genome** (Y or N), **type of genetic element** (:CHRSM, :PLASMID, :MT (mitochondrial chromosome), :PT (chloroplast chromosome), or :CONTIG) or **codon table** (see the corresponding code below).

Example:

+------------+------------+------------+------------+------------+-------------------+
|species     |taxon_id    |  circular  |element_type| codon_table| corresponding_file|
+============+============+============+============+============+===================+
|species_1   |10          |    Y       | :CHRSM     |1           |                   |
+------------+------------+------------+------------+------------+-------------------+
|species_4   |4           |    N       | :CHRSM     |1           |  scaffold_1       |
+------------+------------+------------+------------+------------+-------------------+
|species_4   |4           |    N       | :MT        |1           |  scaffold_2       |
+------------+------------+------------+------------+------------+-------------------+

As you can see for **PF file** (species_4) you can use the column **corresponding_file** to add information for each PF files.

Genetic code for Pathway Tools:

+--------------------+-----------------------------------------------------------------------------------------------+
|Corresponding number|Genetic code                                                                                   |
+====================+===============================================================================================+
|0                   |Unspecified                                                                                    |
+--------------------+-----------------------------------------------------------------------------------------------+
|1                   | The Standard Code                                                                             |
+--------------------+-----------------------------------------------------------------------------------------------+
|2                   | The Vertebrate Mitochondrial Code                                                             |
+--------------------+-----------------------------------------------------------------------------------------------+
|3                   | The Yeast Mitochondrial Code                                                                  |
+--------------------+-----------------------------------------------------------------------------------------------+
|4                   | The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code  |
+--------------------+-----------------------------------------------------------------------------------------------+
|5                   |The Invertebrate Mitochondrial Code                                                            |
+--------------------+-----------------------------------------------------------------------------------------------+
|6                   | The Ciliate, Dasycladacean and Hexamita Nuclear Code                                          |
+--------------------+-----------------------------------------------------------------------------------------------+
|9                   | The Echinoderm and Flatworm Mitochondrial Code                                                |
+--------------------+-----------------------------------------------------------------------------------------------+
|10                  | The Euplotid Nuclear Code                                                                     |
+--------------------+-----------------------------------------------------------------------------------------------+
|11                  | The Bacterial, Archaeal and Plant Plastid Code                                                |
+--------------------+-----------------------------------------------------------------------------------------------+
|12                  | The Alternative Yeast Nuclear Code                                                            |
+--------------------+-----------------------------------------------------------------------------------------------+
|13                  |The Ascidian Mitochondrial Code                                                                |
+--------------------+-----------------------------------------------------------------------------------------------+
|14                  | The Alternative Flatworm Mitochondrial Code                                                   |
+--------------------+-----------------------------------------------------------------------------------------------+
|15                  |Blepharisma Nuclear Code                                                                       |
+--------------------+-----------------------------------------------------------------------------------------------+
|16                  | Chlorophycean Mitochondrial Code                                                              |
+--------------------+-----------------------------------------------------------------------------------------------+
|21                  | Trematode Mitochondrial Code                                                                  |
+--------------------+-----------------------------------------------------------------------------------------------+
|22                  |Scenedesmus obliquus Mitochondrial Code                                                        |
+--------------------+-----------------------------------------------------------------------------------------------+
|23                  | Thraustochytrium Mitochondrial Code                                                           |
+--------------------+-----------------------------------------------------------------------------------------------+

Input files created by mpwt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three input files are created by mpwt. Informations are extracted from the Genbank/GFF/PF file.
myDBName corresponds to the name of the folder and the Genbank/GFF/PF file.
taxonid corresponds to the taxonid in the db_xref of the source feature in the Genbank/GFF/PF.
The species_name is extracted from the Genbank/GFF/PF files.

.. code-block:: text

    **organism-params.dat**
    ID  myDBName
    STORAGE FILE
    NCBI-TAXON-ID   taxonid
    NAME    species_name

    **genetic-elements.dats**
    NAME    
    ANNOT-FILE  gbk_pathname
    //

    **dat_creation.lisp**
    (in-package :ecocyc)
    (select-organism :org-id 'myDBName)
    (let ((*progress-noter-enabled?* NIL))
            (create-flat-files-for-current-kb))

Command Line and Python arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mpwt can be used with the command line:

.. code:: sh

    mpwt -f path/to/folder/input [-o path/to/folder/output] [--patho] [--hf] [--op] [--nc] [--dat] [--md] [--cpu INT] [-r] [--clean] [--log path/to/folder/log] [--ignore-error] [-v]

Optional argument are identified by [].

mpwt can be used in a python script with an import:

.. code:: python

    import mpwt

    folder_input = "path/to/folder/input"
    folder_output = "path/to/folder/output"

    mpwt.multiprocess_pwt(input_folder=folder_input,
			  output_folder=folder_output,
			  patho_inference=optional_boolean,
			  patho_hole_filler=optional_boolean,
              patho_operon_predictor=optional_boolean,
              no_download_articles=optional_boolean,
			  dat_creation=optional_boolean,
			  dat_extraction=optional_boolean,
			  size_reduction=optional_boolean,
			  number_cpu=int,
			  patho_log=optional_folder_pathname,
			  ignore_error=optional_boolean,
			  taxon_file=optional_boolean,
			  verbose=optional_boolean)

+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
| Command line argument   | Python argument                                | description                                                             |
+=========================+================================================+=========================================================================+
|          -f             | input_folder(string: folder pathname)          | Input folder as described in Input data                                 |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          -o             | output_folder(string: folder pathname)         | Output folder containing PGDB data or dat files (see --dat arguments)   |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --patho        | patho_inference(boolean)                       | Launch PathoLogic inference on input folder                             |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --hf           | patho_hole_filler(boolean)                     | Launch PathoLogic Hole Filler with Blast                                |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --op           | patho_operon_predictor(boolean)                | Launch PathoLogic Operon Predictor                                      |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --nc           | no_download_articles(boolean)                  | Launch PathoLogic without loading PubMed citations                      |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --dat          | dat_creation(boolean)                          | Create BioPAX/attribute-value dat files                                 |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --md           | dat_extraction(boolean)                        | Move only the dat files inside the output folder                        |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --cpu          | number_cpu(int)                                | Number of cpu used for the multiprocessing                              |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          -r             | size_reduction(boolean)                        | Delete PGDB in ptools-local to reduce size and return compressed files  |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --log          | patho_log(string: folder pathname)             | Folder where log files for PathoLogic inference will be store           |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --delete       | mpwt.remove_pgdbs(string: pgdb name)           | Delete a specific PGDB                                                  |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          --clean        | mpwt.cleaning()                                | Delete all PGDBs in ptools-local folder or only PGDB from input folder  |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|     --ignore-error      | ignore_error(boolean)                          | Ignore errors and continue the workflow for successful build            |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|     --taxon-file        | taxon_file(boolean)                            | Force mpwt to use the taxon ID in the taxon_id.tsv file                 |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+
|          -v             | verbose(boolean)                               | Print some information about the processing of mpwt                     |
+-------------------------+------------------------------------------------+-------------------------------------------------------------------------+

There is also another argument:

.. code:: sh

    mpwt topf -f input_folder -o output_folder -c cpu_number

.. code:: python

    import mpwt
    mpwt.create_pathologic_file(input_folder, output_folder, cpu_number)

This argument reads the input data inside the input folder. Then it converts Genbank and GFF files into PathoLogic Format files. And if there is already PathoLogic files it copies them.

It can be used to avoid issues with parsing Genbank and GFF files. But it is an early Work in Progress.

Examples
~~~~~~~~

Possible uses of mpwt:

..

    .. code:: sh

        command line

    .. code:: python

        import mpwt
        python script

Create PGDBs of studied organisms inside ptools-local:

..

    .. code:: sh

        mpwt -f path/to/folder/input --patho

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
                patho_inference=True)

Convert Genbank and GFF files into PathoLogic files then create PGDBs of studied organisms inside ptools-local:

..

    .. code:: sh

        mpwt topf -f path/to/folder/input -o path/to/folder/pf
        mpwt -f path/to/folder/pf --patho

    .. code:: python

        import mpwt
        mpwt.create_pathologic_file(input_folder='path/to/folder/input', output_folder='path/to/folder/pf')
        mpwt.multiprocess_pwt(input_folder='path/to/folder/pf', patho_inference=True)

Create PGDBs of studied organisms inside ptools-local with Hole Filler, Operon Predictor and without loading PubMed citations (need Pathway Tools 23.5 or higher):

..

    .. code:: sh

        mpwt -f path/to/folder/input --patho --hf --op --nc --log path/to/folder/log

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
                patho_inference=True,
                patho_hole_filler=True,
                patho_operon_predictor=True,
                no_download_articles=True,
                patho_log='path/to/folder/log')

Create PGDBs of studied organisms inside ptools-local and create dat files:

..

    .. code:: sh

        mpwt -f path/to/folder/input --patho --dat

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
                patho_inference=True,
                            dat_creation=True)

Create PGDBs of studied organisms inside ptools-local.
Then move the files to the output folder.

..

    .. code:: sh

        mpwt -f path/to/folder/input --patho -o path/to/folder/output

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
                            output_folder='path/to/folder/output',
                patho_inference=True)

Create PGDBs of studied organisms inside ptools-local and create dat files.
Then move the dat files to the output folder.

..

    .. code:: sh

        mpwt -f path/to/folder/input --patho --dat -o path/to/folder/output --md


    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(input_folder='path/to/folder/input',
                            output_folder='path/to/folder/output',
                patho_inference=True,
                            dat_creation=True,
                dat_extraction=True)


Create dat files for the PGDB inside ptools-local.
And move them to the output folder.

..

    .. code:: sh

        mpwt --dat -o path/to/folder/output --md

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(output_folder='path/to/folder/output',
                            dat_creation=True,
                dat_extraction=True)

Move PGDB from ptools-local to the output folder:

..

    .. code:: sh

        mpwt -o path/to/folder/output

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(output_folder='path/to/folder/output')

Move dat files from ptools-local to the output folder:

..

    .. code:: sh

        mpwt -o path/to/folder/output --md

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(output_folder='path/to/folder/output',
                dat_extraction=True)


Useful functions
~~~~~~~~~~~~~~~~

- Run the multiprocess Pathway Tools on input folder

..

    .. code:: python

        import mpwt
        mpwt.multiprocess_pwt(input_folder,
                              output_folder,
                              patho_inference=optional_boolean,
                              dat_creation=optional_boolean,
                              dat_extraction=optional_boolean,
                              size_reduction=optional_boolean,
                              number_cpu=int,
                              verbose=optional_boolean)

- Delete all the previous PGDB and the metadata files

..

    .. code:: python

        import mpwt
        mpwt.cleaning()

    This can also be used with a command line argument:

    .. code:: sh

        mpwt --clean

    If you use clean and the argument -f input_folder, it will delete input files ('dat_creation.lisp', 'pathologic.log', 'genetic-elements.dat' and 'organism-params.dat') and the PGDB corresponding to the input folder.

    .. code:: sh

        mpwt -f input_folder --clean

    For example if you have:

    .. code-block:: text

        Folder_input
        ├── species_1
        │   └── species_1.gbk
        ├── species_2
        │   └── species_2.gff
        │   └── species_2.fasta
        ├── species_3
        │   └── species_3.gbk

    And you have in your ptools-local:

    .. code-block:: text

        ptools-local
        ├── pgdbs
            ├── user
                ├── species_1cyc
                │   └── ..
                ├── species_2cyc
                │   └── ..
                ├── species_3cyc
                │   └── ..
                ├── species_4cyc
                │   └── ..

    The command:

    .. code:: sh

        mpwt -f input_folder --clean

    will delete species_1cyc, species_2cyc and species_3cyc but not species_4cyc.

- Delete a specific PGDB

..

    With this command, it is possible to delete a specific PGDB, where pgdb_name is the name of the PGDB (ending with 'cyc'). It can be multiple pgdbs, to do this, put all the pgdb IDs in a string separated by  a ','.

    .. code:: python

        import mpwt
        mpwt.remove_pgdbs(pgdb_name)

    And as a command line:

    .. code:: sh

        mpwt --delete mydbcyc1,mydbcyc2

- Return the path of ptools-local

..

    .. code:: python

        import mpwt
        ptools_local_path = mpwt.find_ptools_path()


- Return a list containing all the PGDBs inside ptools-local folder

..

    .. code:: python

        import mpwt
        list_of_pgdbs = mpwt.list_pgdb()

    Can be used as a command with:

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

If Pathway Tools crashed, mpwt can print some useful information in verbose mode.
It will show the terminal in which Pathway Tools has crashed.
Also, if there is an error in pathologic.log, it will be shown after **=== Error in Pathologic.log ===**.

There is a `Pathway Tools forum <https://ask.pathwaytools.com/questions/>`__ where you can find informations on Pathway Tools errors.

You can also ignore PathoLogic errors by using the argument --ignore-error/ignore_error.
This option will ignore error and continue the mpwt workflow on the successful PathoLogic build.

Output
~~~~~~

If you did not use the output argument, results (PGDB with/without BioPAX/dat files) will be inside your ptools-local folder ready to be used with Pathway Tools.
Have in mind that mpwt does not create the cellular overview and does not used the hole-filler. So if you want these results you should run them after.

If you used the output argument, there is two potential outputs depending on the use of the option **--md/dat_extraction**:

- without --md/dat_extraction, you will have a complete PGDB folder inside your results, for example:

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
    │           └── contains Pathway Tools reports.
    ├── species_2
    ..
    ├── species_3
    ..

- with --md/dat_extraction, you will only have the dat files, for example:

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

- with the **-r /size_reduction** argument, you will have compressed zip files (and PGDBs inside ptools-local will be deleted):

.. code-block:: text

    Folder_output
    ├── species_1.zip
    ├── species_2.zip
    ├── species_3.zip
    ..

Release Notes
-------------

Changes between version are listed on the `release page <https://github.com/AuReMe/mpwt/releases>`__.

Citation
--------

Arnaud Belcour, Clémence Frioux, Meziane Aite, Anthony Bretaudeau, Anne Siegel (2019) Metage2Metabo: metabolic complementarity applied to genomes of large-scale microbiotas for the identification of keystone species. bioRxiv 803056; doi: `https://doi.org/10.1101/803056 <https://doi.org/10.1101/803056>`__.

Acknowledgements
----------------

`Mézaine Aite <https://github.com/mezianeAITE>`__ for his work on the first draft of this package.

`Clémence Frioux <https://github.com/cfrioux>`__ for her work and feedbacks.

Peter Karp, Suzanne Paley, Markus Krummenacker, Richard Billington and Anamika Kothari from the Bioinformatics Research Group of SRI International for their help on Pathway Tools and on Genbank format.

GenOuest bioinformatics (https://www.genouest.org/) core facility for providing the computing infrastructure to test this tool.

All the users that have tested this tool.
