import pandas as pd

from ast import literal_eval

target_file="p1_targets.fasta"
drug_file="p1_drugs.csv"
interaction_file="p1_interaction.csv"



with open('ligands_can.txt') as l:
    ligands = l.read()
 
ligands = literal_eval(ligands)

with open('proteins.txt') as p:
    proteins = p.read()
 
proteins = literal_eval(proteins)

ligand_keys = list(ligands.keys())
proteins_keys = list(proteins.keys())


# interaction_file
affinity_values = pd.read_csv("drug-target_interaction_affinities_Kd__Davis_et_al.2011v1.txt", sep=" ", names=proteins_keys, index_col=None)
affinity_values['PubChem CID'] = ligand_keys
affinity_values = affinity_values.set_index('PubChem CID')

affinity_values.to_csv(interaction_file, sep="\t")


# drug_file
pd.DataFrame({'a': ligands.values(), 'b': ligands.keys()}).to_csv(drug_file, sep=" ", header=False, index=False)


# target_file
f = open(target_file, 'w')
for key in proteins_keys:
	f.write(">"+key+"\n")
	i = 0
	while i < len(proteins[key]):
		if i % 40 != 39:
			f.write(proteins[key][i])
			i+=1
		else:
			f.write(proteins[key][i]+"\n")
			i+=1
	if i % 40 != 0:
		f.write("\n")