import pandas as pd
import numpy as np


# this file serves the very simple purpose of repeating the clean up process of the BindingDB data as done in:
# https://link.springer.com/article/10.1007/s10822-021-00404-7
# As of today (16-08-21) Tanoori et al. might have done 

# Homo sapiens organisms were selected. Then, those lacking UniProt/PubChem IDs were removed.
# Among the remained interactions, those having multiple/unknown Kd values were also set aside to result in a subset of 11,689 interactions. 
# Since they were using the 2018 dataset repeating the same steps will now result in a larger final subset.

# I also removed all ligands without ChEMBL IDs 

# cols = ['BindingDB Reactant_set_id', 'Target Name Assigned by Curator or DataSource', 'Target Source Organism According to Curator or DataSource', 'Kd (nM)']
# cols = ['Ligand SMILES','PubChem CID','PubChem SID','BindingDB Target Chain  Sequence','Kd (nM)']
# cols = ['PubChem CID','PubChem SID']
cols = ['Target Name Assigned by Curator or DataSource','PubChem CID','UniProt (SwissProt) Primary ID of Target Chain','ChEMBL ID of Ligand','Target Source Organism According to Curator or DataSource','Ligand SMILES','BindingDB Target Chain  Sequence','Kd (nM)']
path = '/media/lanorius/Elements/data/BindingDB_All.tsv'

bindingdb_all = pd.read_csv(filepath_or_buffer=path, sep='\t', usecols = cols, error_bad_lines=False)
print(bindingdb_all.shape)
# bindingdb_all = pd.read_csv(filepath_or_buffer=path, sep='\t', error_bad_lines=False)

# cond = xor(bindingdb_all['PubChem CID'] == 'NaN', bindingdb_all['PubChem SID'] == 'NaN')
# bindingdb_all['que'] = np.where(cond, True)



# Homo sapiens organisms were selected.
bindingdb_all = bindingdb_all[bindingdb_all['Target Source Organism According to Curator or DataSource'] == "Homo sapiens"]

# Then, those lacking UniProt/PubChem IDs were removed.
bindingdb_all = bindingdb_all.dropna( how='any', subset=['PubChem CID'])
bindingdb_all = bindingdb_all.dropna( how='any', subset=['UniProt (SwissProt) Primary ID of Target Chain']) # this might need expanding

# I thought about additionally removing all rows with missing ChEMBL IDs fo the Ligands
# bindingdb_all = bindingdb_all.dropna( how='any', subset=['ChEMBL ID of Ligand'])

# Among the remained interactions, those having multiple/unknown Kd values were also set aside.
bindingdb_all['Kd (nM)'] = pd.to_numeric(bindingdb_all['Kd (nM)'],errors='coerce')
bindingdb_all = bindingdb_all.dropna( how='any', subset=['Kd (nM)'])


# bindingdb_all = bindingdb_all[bindingdb_all['Kd (nM)'].map(type) == str] # not needed since some have proper values cast as strings
# kicks out way too many rows because some cells hold multiple values, and some are saved as strings


# bindingdb_all = bindingdb_all[bindingdb_all['Kd (nM)'].map(float)] 
# needs adjustment to keep the strings with single values and only kicks out the rows with multiple values

print(bindingdb_all)

# bindingdb_all.to_csv('cleaned_bindingDB_all.tsv', sep='\t', header=True, index=True, mode='w')



'''
print(bindingdb_all)
# print(set([type(x) for x in bindingdb_all['Kd (nM)'].tolist()]))

# print([x for x in a if type(x) == float])


# bindingdb_all.to_csv('Tanoori_clean_BindingDB_All.tsv.tsv', sep = '\t', index=False)
'''


'''
a = [['10', '1.2'], ['15', 'NaN']]
df = pd.DataFrame(a, columns=['one', 'two'])


cond = xor(df['one'] == 'NaN', df['two'] == 'NaN')
df['que'] = np.where(cond, True)
'''
