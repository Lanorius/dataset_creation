# dataset_creation
This repo will hold part of my work on my master thesis. I am creating a dataset for Drug-Target-Interaction (DTI) problems.


WORK IN PROGRESS

data_cleaning:
process_inputs_to_clean.py is the main script
The current plan is for it to do the following steps
	- Take whichever version of bindingDB_all.tsv is the most current one
	- Perform a cleanup of the data (this is broadly defined, since I am currently researching which steps are resonable, but will also result in a large enough dataset for DTI problems)
	- Allow specificiations in the config file to choose how much sequence and/or compound similarity will be allowed
	
In the case that more sources of DTI data are found the script will be changed to work with them as well if the time allows it.
Idealy multipe sources of data will be combineable

DTI_set_creation:
cleaned_tsv_to_D_T_I.py takes the cleaned tsv and devides it into:
	-a .fasta file for the proteins
	-a file for the ligands
	-an interaction file 
The first two files can then be used to create embeddings or other representations for the final DTI prediction. 

building_pieces:
This folder holds scripts that I used to piece together the actual clean up script. 
Statistics on the data were also created using these scripts.
