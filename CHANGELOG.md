# Changelog

# mpwt 0.8.0 (2022-09-21)

## Add:

- new option Complex Inference `--cp` (issue #75).
- function get_ptools_version to extract Pathway Tools version.
- `-dump-flat-files-biopax` option with Pathway Tools 26.0 (issue #76).
- the possibility to use already present input files (issue #79).

## Fix:

- numerous issues with give_permission.
- issue in mpwt figure.
- issue if there are files in input folder (issue #80).

## Modify:

- import pwt_wrapper log creation.

# mpwt 0.7.2 (2022-04-15)

## Fix:

- an issue with log creation with Pathway Tools 26.0.

# mpwt 0.7.1 (2022-03-18)

## Add:

- support for specifying reference PGDB (with taxon_id file) to be used in PathoLogic.
- warning message when there is a missing pathologic.log file during log file creation (before mpwt crahses with a python error associated with `species_pathologic_informations`).

## Fix:

- an issue where mpwt does not stop if there is an error in pwt_input_files.
- input files not created if mpwt uses PGDB from ptools-local.

## Modify:

- update readme.

# mpwt 0.7.0 (2022-02-03)

This version should be compatible with Pathway Tools 25.5.

The behaviour of mpwt has been changed in this version. Before 0.7.0 mpwt will wait for all the PathoLogic processes to end before going to the next step. But if there was an error with one process, it will stop and no outputs were created even for successful processes.
Now, mpwt run the processes independently. This means that for example if you have 3 organisms (Org_A, Org_B and Org_C). If the process for Org_A fails but the processes for Org_B and Org_C succeed, then you will have an output folder containing results for Org_B and Org_C (whereas in older version no results were created).

## Add:

- mpwt can take as input PathoLogic files without fasta. It will output a warning message but it should be usable. This should make mpwt compatible with PathoLogic files created by EsMeCaTa.
- new option --permission/permission to choose the permission level for PGDB in ptools-local or output folder (useful when working with a container in a cluster).
- changelog file.
- citation in readme.

## Fix:

- issue with missing tmp path in run_move_pgdb.
- issue with topf (especially when using a GFF file some fields like dbxref were not correctly used).
- issue that prevents mpwt from killing Pathway Tools process that enters the Lisp Listener.
- issue with pathway prediction score not being correctly return to its previous value.
- String not Boolean and path incorrect in Readme (issue #70).
- Issue with Pathway Tools 25.5 when checking log (issue #71). The way of parsing the 'PGDB contains' line has been modified to be more robust. The parsing with str.split() has been replaced with a parsing with regex.
- numerous typos.

## Modify:

- refactor mpwt code to run each process independently from the other  (issue #68). Creation of the functions run_mpwt (to launch a complete run of mpwt on one organism) and independent_mpwt (to launch multiple run_mpwt on each inputs).
- refactor log creation. Instead of a function taking as input a list of dictionary (check_pwt) use 2 functions (check_mpwt_pathologic_runs and extract_pathologic). The first one takes as input a list of the paths to the input folder containing pathologic.log file and an output folder. Then it uses the second function (extract_pathologic) for each pathologic file to extraction the informations.
- replace docopt with argparse.
- update readme.
- year in license code.
- update mpwt workflow svg according to the changes in mpwt code. Put a white background (useful when using a GitHub dark theme).

# mpwt 0.6.3 (2021-04-02)

This version should be compatible with Pathway Tools 25.0.

## Fix:

- Issue with Pathway Tools 25.0 when checking log (issue #67).

# mpwt 0.6.2 (2021-04-01)

## Fix:

- Issue with encoding of flat_files_creation.log (issue #66).

## Modify:

- Use setup.cfg for version.

# mpwt 0.6.1 (2021-03-16)

## Modify:

- replace GPL license by LGPL license.

# mpwt 0.6.0 (2020-12-14)

**Warning**:
In this version, some option names have been changed:
- "--dat" ("dat_creation" in python code) has been replaced by "--flat" ("flat_creation" in python code) to better reflect the attribute-values flat files created by this option.

## Add:

- new options to extract other files to the ouput folder (issue #63):
    - "--mx" ("xml_creation") to extract XML files created by MetaFlux.
    - "--mo" ("owl_extraction") to extract owl files.
    - "--mc" ("col_extraction") to extract tabulated files (".col" extension).

## Fix:

- issue with "-r" (size_reduction) option: all files were extracted instead of only requested output files.

## Modify:

- replace "--dat" ("dat_creation") by "--flat" ("flat_creation") to better reflect data created by this command.
- use PGDB ID as output name when using mpwt without input_folder (like "mpwt -o output_folder").
- update readme.
  
# mpwt 0.5.9 (2020-11-20)

## Add:

- \_\_version\_\_ in mpwt init file to handle version (issue #59).

## Fix:

- error when creating dat files without input folder (issue #61).

## Modify:

- refactor how paths are handled (issue #60). This is a first step to the compatibility of mpwt on Windows. But this will need more tests and modification.

# mpwt 0.5.8 (2020-10-12)

mpwt has now a new dependency with **chardet** (linked to issue #56).

## Fix:

- issue with encoding when writing the terminal logs (issue #56).
- issue with Pathway Tools pop-up showing after attribute-value flat files creation (issue #57).

# mpwt 0.5.7 (2020-10-01)

## Add:

- an option (--tp/patho_transporter_inference) to use Transporter Inference Parser with PathoLogic (issue #53).
- a badge showing last Pathway Tools release compatible with mpwt.
- a picture trying to show how mpwt works.
- clean option when using topf.
- GitHub Actions to release mpwt on PyPI with GitHub Release.

## Fix:

- no join and close of Pool map if an error occurred (issue #52).
- issue when error can occur after PGDB creation but not take into account (issue #54).

## Modify:

- replace Pool map by starmap (issue #51).

# mpwt 0.5.6 (2020-07-28)

## Add:

- a check for unfinished builds. If during a previous run a build failed, mpwt will delete it and relaunch it (issue #49).
- a badge with a link to the preprint.

## Fix:

- numerous issues in topf (issue #30).
- numerous issues with taxon_id file when using circular, codon_table or element_type. 
- numerous issues with raise error messages not showing.
- typos.

# mpwt 0.5.5 (2020-05-27)

## Add:

- support for .fsa extension (issue #47).
- a license (issue #45).
- a --version argument to mpwt command-line (issue #44).
- pseudogene linked to mRNA in topf (issue #30).

## Fix:

- issue with encoding in mpwt_wrapper.
- issue when --md could be used without -o argument (issue #46).

Thanks to @cfrioux for her work on this release.

# mpwt 0.5.4 (2020-04-10)

## Add:

- support to .gbff extension for genbank files (issue #39).
- an option to use pathway prediction score of Pathway Tools (issue #42).
- a logger for mpwt (look at [For developer](https://github.com/AuReMe/mpwt#for-developer) for more informations).
- a check for ptools-init.dat file.

## Fix:

- issue when checking the input files (issue #37).
- issue when checking the existence of taxon_id.tsv (issue #38).
- issue if there is multiple files in the same input folder (issue #40).
- mpwt killing Pathway Tools process to create attribute-values files even if there is only a not fatal error (issue #41).
- numerous issues with to_pf (issue #30). Add DBLINK and a check for gene name length. Rewrite how CDS and RNAs are handled. But this option needs more testing and modification in 0.5.5.
- typos in Readme.

## Modify:

- pwt_wrapper: refactoring the code and add a new log (pwt_terminal.log) containing Pathway Tools log from the terminal.

## Check:

- no_download_articles/--nc option is currently not working.

# mpwt 0.5.3 (2020-01-09)

## Add:

- an option to create taxon_id.tsv files from GenBank and GFF (issue #35).
- an option to use operon predictor of Pathway Tools (issue #33).
- an option to avoid downloading PubMed entries using a parameter in ptools-init.dat (added in Pathway Tools 23.5). But it needs more testing and maybe modification in 0.5.4 (issue #34).
- a test for to_pf using GenBank and GFF from test data (issue #30).
- errors and warnings counts in log files.

## Fix:

- numerous issues with to_pf (issue #30). But this option needs more testing and modification in 0.5.4.
- incompatibility with Pathway Tools 23.5.
- typos in Readme.

# mpwt 0.5.2 (2019-10-17)

## Add:

- to_pf argument, allowing to convert GBK or GFF into PF files (issue #30). Need more polishing in 0.5.3.
- preprint  citation in Readme.

## Fix:

- error message when a GFF has no region feature (issue #31).

# mpwt 0.5.1 (2019-07-31)

## Add:

- taxon_id.tsv can be used to give other informations (element type, codon table, circularity of genome) to Pathway Tools (issue #28).

## Fix:

- test not working on Pathway Tools 23.0.
- error message if user do not give an integer for the number of CPUs.

## Modification:

- verbose and logging interactions.

# mpwt 0.5.0 (2019-07-02)

**Warning**:
In this version, the behaviours of some arguments have been changed.
- the -r/size_reduction argument will delete PGDB inside ptools-local and compressed the results in zip (issue #26). 
- the --clean argument (when used with -f/input_folder argument) will delete only PGDBs corresponding to species in the input_folder (issue #23).

## Add:

- support for Pathologic Format (PF) file (issue #19).
- taxon_id.tsv file to manage taxon_id for species. With the argument --taxon-id/taxon_file argument, you can force mpwt to use taxon from taxon_id.tsv for all type of files (Genbank, GFF or PF).
- a new argument --ignore-error/ignore_error to ignore PathoLogic failed builds (issue #21). With this argument, mpwt will continue to run even if PathoLogic have crashed for some species. It will perform the rest of the workflow for the successful build. 
- time measure of each steps of mpwt. They will be printed if you use -v/verbose or they can be accessed at the end of the log_error.txt file (created with --log).

## Fix:

- interaction between --clean and the other arguments.
- error if -r/size_reduction is used without output_folder argument.
- a typo error with the name of a function (remove_pgdbs).


## Modification:

- split multipwt.py in multiple modules to ease reading (issue #25).
- --clean behaviour when used with -f/input_folder (issue #23) to delete only input species PGDB.
- -r/size_reduction behaviour: it will delete the PGDB from ptools-local and it will move the results into a compressed zip file.
- add table and list in Readme.

# mpwt 0.4.2.4 (2019-06-07)

## Fix:

- wrong argument number given in the fix of the error with --clean when used with -v and --cpu (issue #20).

# mpwt 0.4.2.3 (2019-06-07)

## Fix:

- error with --clean when used with -v and --cpu (issue #20).
- error if PathoLogic build is aborted (issue #22).
- typos.

# mpwt 0.4.2.2 (2019-04-18)

## Fix:

- issue with BioPAX/dat creation (issue #17).
- issue with error report with BioPAX/dat creation.

## Modification:

- use by default 1 CPU, instead of using all CPUs available.

# mpwt 0.4.2.1 (2019-03-28)

## Fix:

- error with outdated argument in cleaning_input.
- several typos.

# mpwt 0.4.2 (2019-03-21)

## Add:

- logging instead of print (issue #12 ). 
- mpwt can create BioPAX/attribute-values files when PGDBs are already in ptools-local but not in output folder (issue #13  ). In the same run, mpwt can launch PathoLogic on a set of species having no PGDBs and then launch BioPAX/attribute-values files creation on this set of species and on other species with PGDBs but with no results in the output folder.
- fatal errors from pathologic.log are now printed (issue #14  ).

## Fix:

- number_cpu argument not used by cleaning function.

# mpwt 0.4.1 (2019-03-18)

## Add:

- option to use Pathway-Tools Hole-Filler (issue #9 ). 
- error message if Pathway-Tools is not installed or not in PATH (issue #10 ).
- fasta requirement for GFF file (issue #7 ).

## Fix:

- mpwt could take as input hidden files/folders (issue #8 ). Add test to check this case.
- issue with the name of the species output folder when giving an output folder but not an input folder.

## Modification:

- remove global variables.
- version of dependencies.

# mpwt 0.4.0 (2019-02-25)

First release on GitHub.

**Warning**:
In this version mpwt arguments have been changed to make the tool more flexible. Now to launch PathoLogic inference, you need the --patho argument. If you want to create BioPAX/attribute-values file, you need to use --dat and to move only the attribute-values dat files, you must use --md option. More informations in the readme.

## Add:

- mulitprocessing when creating input files, moving result files and deleting PGDBs (issue #1 ).

- better Pathway-Tools error handling (issue #2 ).

- more flexible tool (issue #3 ). Change how the argument work to make the tools more flexible.

- accept GFF file (issue #5 ). Only partial support.

- an argument to show the PGDBs in ptools-local (--list).


## Fix:

- issue with pop-up (issue #4 ). To avoid all the pop-ups, you need to use at least the version 22.5 of Pathway-Tools.

- when there is multiple db_xref qualifiers in the 'source' feature of the Genbank file, take the correct one.

- mpwt stops when there is PGDBs in ptools-local folder.
