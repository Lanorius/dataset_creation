# My Master's Thesis: Drug-Target-Interaction (DTI) Dataset Creation


This repo compliments the DTI prediction algorithm in:
[Data Creator](https://github.com/Lanorius/binding-affinity-prediction): Predicits the datasets that are created here.

#### How to use the Data Cleaner
1. Clone this repo
2. Ensure the requirements are met
3. check in the src/config.ini file if all files are chosen correctly
4. run uisng "python main.py"

#### Databases

[BindingDB](https://www.bindingdb.org/bind/index.jsp): Database of measured binding affinity

<!--
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
-->


<!--		
# Prediction-Of-Binding-Affinity

Our aim was to predict the binding affinity of compounds and proteins by using a fully connected neural network.
To encode both types of molecules we used ChemVAE for the compounds and ProT5 for the proteins. 

![plot](./images/mymodel.png)


#### How to use our Model
1. Clone this repo
2. Ensure the requirements are met
3. check in the src/config.ini file if all files are chosen correctly
	*if you wish to run a special task, set the relevant task to true
	*if you wish to use overtraining, you have to set the parameters in the special params section
4. run uisng "python src/binding_prediction.py"

Several clustered and unclustered datasets are available in the data folder.

This algorithm currently only works with pKd scores. Leave the general section as it is. When the general setting is set to davis
the model works on both the pKd scores from the Davis set as well as those from the BindingDB.
The compound setting should not be changed as well, unless another way of encoding the compounds is implemented. 


#### Embeddings

[ChemVAE](https://github.com/aspuru-guzik-group/chemical_vae): Variational auto encoder that was used to create compound embeddings

[T5 Embeddings](https://github.com/agemagician/ProtTrans): Resource for creating protein sequence embeddings

#### Databases

[BindingDB](https://www.bindingdb.org/bind/index.jsp): Database of measured binding affinity
-->
