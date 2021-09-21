from process_inputs_to_DTI import parse_config
import pandas as pd

'''
This file takes a cleaned DTI tsv created by input_to_clean_tsv.py and divides it into a .fasta file for the proteins,
a file with ligands, and an interaction file.

NEXT STEP: All files need redundancy reduction.
'''

# will be part of the config
path = '/media/lanorius/Elements/data/cleaned_frame.tsv'
sep = '\t'

# the following could be specified using the same config file that input_to_clean_tsv is using
protein_IDs = 'UniProt (SwissProt) Primary ID of Target Chain'
protein_sequence = 'BindingDB Target Chain  Sequence'
ligand_IDs = 'ChEMBL ID of Ligand'
ligand_SMILE = 'Ligand SMILES'  # maybe something else than the SMILE will be used in the future
interaction_value = 'Kd (nM)'

f = open('sequences.fasta', 'w')
d = open('drug_file.txt', 'w')

file = pd.read_csv(filepath_or_buffer=path, sep=sep, engine='python')
print(file)

for index, row in file.iterrows():

    f.write(">"+row[protein_IDs]+"\n")
    i = 0
    while i < len(row[protein_sequence]):
        if i % 40 != 39:
            f.write(row[protein_sequence][i])
            i += 1
        else:
            f.write(row[protein_sequence][i])
            f.write("\n")
            i += 1
    if i % 40 != 0:
        f.write("\n")

    d.write(row[ligand_IDs]+"\t"+row[ligand_SMILE]+"\n")

interactions = file.pivot_table(index=ligand_IDs, columns=protein_IDs, values=interaction_value, aggfunc='sum')
print(interactions)
interactions.to_csv('interactions.csv', sep='\t')

f.close()
d.close()
