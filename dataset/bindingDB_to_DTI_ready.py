from data_loading.process_inputs import parse_config
import pandas as pd
import ast   # to go from string to list while parsing calnames
import math  # to calculate better chunk sizes

# One idea would be to expect the user to do a very basic level of preprocessing
# This document will be much more powerful if we can expect the input to be the following:
#	-a csv or tsv file including columns with portein IDs, protein sequences, compound IDs, compound SMILES/Fingerprints, a parameter of affinity (Kd, Ki, ic50)
#	-for each file these columns have to be specified in the config file

# add a config file that specifies:
# The name of all input files
# The wanted parameters: sequence similarity, smile similarity, affinity values...
# The name of the output file
# Whatever else comes to mind with time


# Loading the data
name_of_file, specifications, params = parse_config()
# print(name_of_file['bindingDB_file'])
# print(params['sequence_similarity'])


cols_ = ['Target Name Assigned by Curator or DataSource', 'PubChem CID', 'UniProt (SwissProt) Primary ID of Target Chain', 'ChEMBL ID of Ligand', 'Target Source Organism According to Curator or DataSource', 'Ligand SMILES', 'BindingDB Target Chain  Sequence', 'Kd (nM)']
# belongs in the config, here for testing reasons

path = name_of_file['bindingDB_file']
sep = specifications['separator_one']
cols = ast.literal_eval(specifications['colnames_one'])

# print(cols == cols_)
print(type(cols_))
print(type(cols))

# chunksize = math.ceil(len(list(open(path)))/5)
# print(chunksize)
chunksize = 5*10**5
for chunk in pd.read_csv(filepath_or_buffer=path, sep=sep, chunksize=chunksize, usecols=cols, on_bad_lines='skip', engine='python'):
    chunk = chunk[chunk['Target Source Organism According to Curator or DataSource'] == "Homo sapiens"]
    print(chunk)

'''
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
'''


# bindingdb_all = bindingdb_all[bindingdb_all['Kd (nM)'].map(type) == str] # not needed since some have proper values cast as strings
# kicks out way too many rows because some cells hold multiple values, and some are saved as strings


# bindingdb_all = bindingdb_all[bindingdb_all['Kd (nM)'].map(float)] 
# needs adjustment to keep the strings with single values and only kicks out the rows with multiple values



# Preprocessing
# If you use multiple sources make sure each Drug-Target pair is unique
# Remove Sequences based on sequence similarity

# Remove SMILES based on (figure this out)

# Clustering?

# You could write an object for each source of data that you will have in the end

# The goal is to have a single dataset in the end that fullfills the given parameters

