# dataset_creation

This repository contains the code for my Master Thesis. The idea is to use common tools in Bioinformatics, like RDKit and CD-Hit, to process raw data, in hope of improving its quality for Drug-Target-Interaction prediction.

WORK IN PROGRESS

Part 1 data_cleaning:

As a prerequist, the user has to provide a .csv or .tsv file that has the following columns:
protein_IDs
protein_sequence
ligand_IDs
ligand_SMILE
interaction_value (like Kd, Ki, Kiba Scores)

Some databases have additional columns that specify sequences are multiprotein complexes are involved in binding. These should also be removed beforehand.


raw_input_to_raw_DTI_ready.py:
In the config file the input .csv/.tsv and the expected columns have to be specified, as well as the names of the output files.

The current version removes all missing data and throws out all rows that have ambiguous entries.

Four files are given as output:
A new .tsv file that includes only rows with usable raw data. This means data that still requires redundancy reduction for both the Drugs and the Targets.
A .fasta file with the protein sequences.
A file with the ligand IDs and SMILES.
Another .tsv file with the interaction matrix.


Part 2 double clustering:
		+ Removing sequences that are too similar and drawing representatives (using CD-Hit, work in progress) 
		+ Removing SMILES that have too similar fingerprints by calculating Jaccard distances and drawing representatives (work in progress)
		+ Averaging over interactions where multiple representative interactions take place
