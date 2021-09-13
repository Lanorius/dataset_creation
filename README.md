# dataset_creation
I am creating a dataset for Drug-Target-Interaction (DTI) problems.

WORK IN PROGRESS

data_cleaning:
process_inputs_to_clean.py is the main script
The current plan is for it to do the following steps
	- Take a version of bindingDB_all.tsv or similar data 
	- Perform a cleanup of the data 
	(this is broadly defined, since I am currently researching which steps are reasonable, but will also result in a large enough dataset for DTI problems)
		+ Removing sequences that are too similar (using CD-Hit, work in progress)
		+ Removing SMILES that have too similar fingerprints (work in progress)
	- Allow specification in the config file to specify how much similarity is allowed
	
In the case that more sources of DTI data are found the script will be changed to work with them as well if the time allows it.
Ideally, multiple sources of data will be combinable.


DTI_set_creation:
cleaned_tsv_to_D_T_I.py takes the cleaned tsv and divides it into:
	-a .fasta file for the proteins
	-a file for the ligands
	-an interaction file 
The first two files can then be used to create embeddings or other representations for the final DTI prediction. 

building_pieces:
This folder holds scripts that I used to piece together the actual clean up script. 
Statistics on the data were also created using these scripts.
