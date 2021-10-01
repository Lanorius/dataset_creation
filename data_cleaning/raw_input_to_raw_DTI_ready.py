from process_inputs import parse_config
import pandas as pd
import ast   # to go from string to list while parsing col-names
import math  # to calculate better chunk sizes
import os  # to remove temp_drugs.txt

'''
One idea would be to expect the user to do a very basic level of preprocessing
This document will be much more powerful if we can expect the input to be the following:
    -a csv or tsv file including columns with protein IDs, protein sequences, compound IDs, 
    compound SMILES/Fingerprints, a parameter of affinity (Kd, Ki, IC50, EC50)
    -for each file these columns have to be specified in the config file

add a config file that specifies:
The name of all input files
The wanted parameters: sequence similarity, smile similarity, affinity values...
The name of the output file
Whatever else comes to mind
'''

# PART 1

# Loading the data
raw_transformation, files, output, file_specifications = parse_config()
# files is writen as plural, however for now it only works with a single file
# params will very likely be moved to another file

cleaned_file = output['cleaned_frame']

if raw_transformation:
    print("Part 1: Performing raw transformation.")
else:
    print("Part 1: Skipping raw transformation.")

if raw_transformation:
    path = files['bindingDB_file']
    sep = file_specifications['separator_one']
    cols = ast.literal_eval(file_specifications['colnames_one'])

    '''
    Currently selecting:
    Those lacking UniProt/PubChem IDs were removed
    Those lacking ChEMBL IDs were removed
    Those having multiple/unknown Kd values were removed
    '''

    cleaned_frame = pd.DataFrame(columns=cols)

    chunksize = math.ceil(len(list(open(path)))/5)  # allows handeling large sets of input data
    for chunk in pd.read_csv(filepath_or_buffer=path, sep=sep, chunksize=chunksize, usecols=cols, on_bad_lines='skip',
                             engine='python'):
        # chunk = chunk[chunk['Target Source Organism According to Curator or DataSource'] == "Homo sapiens"] # removed from the config for now
        chunk = chunk.dropna(how='any', subset=['PubChem CID'])  # you do have Uniprot IDs already
        chunk = chunk.dropna(how='any', subset=[file_specifications['protein_IDs']])
        chunk = chunk.dropna(how='any', subset=[file_specifications['ligand_IDs']])
        chunk[file_specifications['interaction_value']] = pd.to_numeric(chunk[file_specifications['interaction_value']], errors='coerce')
        chunk = chunk.dropna(how='any', subset=[file_specifications['interaction_value']])
        cleaned_frame = pd.concat([cleaned_frame, chunk])
    # this needs to be adjusted to the config file
    cleaned_frame.to_csv(cleaned_file, sep='\t')

# PART 2
# move the other file here

print("Part 2: Creating necessary files.")

path = cleaned_file
sep = '\t'  # should come from config

protein_IDs = file_specifications['protein_IDs']
protein_sequence = file_specifications['protein_sequence']
ligand_IDs = file_specifications['ligand_IDs']
ligand_SMILE = file_specifications['ligand_SMILE']
interaction_value = file_specifications['interaction_value']

f = open(output['fasta_file'], 'w')
d = open('temp_drugs.txt', 'w')

file = pd.read_csv(filepath_or_buffer=path, sep=sep, engine='python')

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

    d.write(row[ligand_SMILE]+" "+row[ligand_IDs]+"\n")

# ensures the output is a csv file
# all duplicate lines are removed here from the drugs as well
temp_drugs = pd.read_csv('temp_drugs.txt', sep=' ').drop_duplicates(keep='first').reset_index()
temp_drugs = temp_drugs.iloc[:, 1:]
temp_drugs.to_csv(output['drug_file'], sep=' ', index=False)
os.remove('temp_drugs.txt')

interactions = file.pivot_table(index=ligand_IDs, columns=protein_IDs, values=interaction_value, aggfunc='sum')
interactions.to_csv(output['interaction_file'], sep='\t')

f.close()
d.close()
