# Pathway-tools multiprocessing script

Python script allowing to run pathway-Tools in multiprocess. To do this the script takes a folder containing sub-folders as input. Each sub-folder contains a genbank file.  

    Folder_input
    ├── Folder_for_species_1
    │   └── Genbank_species_1
    ├── Folder_for_species_2
    │   └── Genbank_species_2
    ├── Folder_for_species_3
    │   └── Genbank_species_3
    │

Pathway-Tools will run on each genbank file. It will create an output folder containing all the result files from the PathoLogic inference for each species.

Used in [Pathway-Tools Multiprocessing Docker](https://github.com/ArnaudBelcour/pathway-tools-multiprocessing-docker).
