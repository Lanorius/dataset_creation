from process_inputs_to_clean import parse_config
import pandas as pd
import ast   # to go from string to list while parsing calnames
import math  # to calculate better chunk sizes

'''
One idea would be to expect the user to do a very basic level of preprocessing
This document will be much more powerful if we can expect the input to be the following:
    -a csv or tsv file including columns with portein IDs, protein sequences, compound IDs, 
    compound SMILES/Fingerprints, a parameter of affinity (Kd, Ki, ic50)
    -for each file these columns have to be specified in the config file

add a config file that specifies:
The name of all input files
The wanted parameters: sequence similarity, smile similarity, affinity values...
The name of the output file
Whatever else comes to mind
'''

# Loading the data
name_of_file, name_of_output, specifications, params = parse_config()
# print(name_of_file['bindingDB_file'])
# print(params['sequence_similarity'])


path = name_of_file['bindingDB_file']
output = name_of_output['file_name']
sep = specifications['separator_one']
cols = ast.literal_eval(specifications['colnames_one'])

'''
Currently selecting:
Homo sapiens
Those lacking UniProt/PubChem IDs were removed
Those lacking ChEMBL IDs were removed
Those having multiple/unknown Kd values were removed
'''

cleaned_frame = pd.DataFrame(columns=cols)
# cleaned_frame will later be used to create a Drug, Target (fasta) and an interaction file

chunksize = math.ceil(len(list(open(path)))/5)  # allows handeling large sets of input data
for chunk in pd.read_csv(filepath_or_buffer=path, sep=sep, chunksize=chunksize, usecols=cols, on_bad_lines='skip',
                         engine='python'):
    chunk = chunk[chunk['Target Source Organism According to Curator or DataSource'] == "Homo sapiens"]
    chunk = chunk.dropna(how='any', subset=['PubChem CID'])
    chunk = chunk.dropna(how='any', subset=['UniProt (SwissProt) Primary ID of Target Chain'])  # might need expanding
    chunk = chunk.dropna(how='any', subset=['ChEMBL ID of Ligand'])
    chunk['Kd (nM)'] = pd.to_numeric(chunk['Kd (nM)'], errors='coerce')
    chunk = chunk.dropna(how='any', subset=['Kd (nM)'])
    cleaned_frame = pd.concat([cleaned_frame, chunk])

cleaned_frame.to_csv(output, sep='\t')

# Preprocessing
# If you use multiple sources make sure each Drug-Target pair is unique
# Remove Sequences based on sequence similarity

# Remove SMILES based on (figure this out)

# Clustering?

# You could write an object for each source of data that you will have in the end

# The goal is to have a single dataset in the end that fullfills the given parameters
