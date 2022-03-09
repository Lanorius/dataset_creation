# My Master's Thesis: Drug-Target-Interaction (DTI) Dataset Creation


This repo compliments the DTI prediction algorithm in:
[Data Creator](https://github.com/Lanorius/binding-affinity-prediction)

#### How to use the Data Cleaner
1. Clone this repo
2. Download the MayaChemTools collection, and place the folder next to the folder of this repo (or you can specify another location in the config file)
3. Install CD-Hit
4. Ensure the requirements are met, we found that the best way is to create an RDKit environment and install the other requirements to it
5. check in the src/config.ini file if all files are chosen correctly, including the raw data file
6. run by using "python main.py"

#### Databases

[BindingDB](https://www.bindingdb.org/bind/index.jsp): Database of measured binding affinity

#### RDKit

[RDKit](https://www.rdkit.org/): Open-Source Cheminformatics Software, install using the information in the following link https://www.rdkit.org/docs/Install.html

#### RDKit additional scripts

[MayaChemTools](http://www.mayachemtools.org/): collection of Perl and Python scripts, modules, and classes to support a variety of day-to-day computational discovery needs

#### CD-Hit

[CD-Hit](http://cd-hit.org): a widely used program for clustering and comparing protein or nucleotide sequences
