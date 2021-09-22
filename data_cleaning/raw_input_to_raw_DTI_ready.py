from process_inputs import parse_config
import pandas as pd
import ast   # to go from string to list while parsing calnames
import math  # to calculate better chunk sizes
import os  # to specify if the user wants to keep the cleaned TSV

'''
One idea would be to expect the user to do a very basic level of preprocessing
This document will be much more powerful if we can expect the input to be the following:
    -a csv or tsv file including columns with portein IDs, protein sequences, compound IDs, 
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
raw_transformation, files, output, specifications, params = parse_config()
# files is writen as plural, however for now it only works with a single file

print(raw_transformation)

if raw_transformation:
    path = files['bindingDB_file']
    save_output = output['save_file']
    output_file = output['file_name']
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

    chunksize = math.ceil(len(list(open(path)))/5)  # allows handeling large sets of input data
    for chunk in pd.read_csv(filepath_or_buffer=path, sep=sep, chunksize=chunksize, usecols=cols, on_bad_lines='skip',
                             engine='python'):
        chunk = chunk[chunk['Target Source Organism According to Curator or DataSource'] == "Homo sapiens"]
        chunk = chunk.dropna(how='any', subset=['PubChem CID'])
        chunk = chunk.dropna(how='any', subset=['UniProt (SwissProt) Primary ID of Target Chain'])  # might need expanding
        chunk = chunk.dropna(how='any', subset=['ChEMBL ID of Ligand'])
        chunk['EC50 (nM)'] = pd.to_numeric(chunk['EC50 (nM)'], errors='coerce')
        chunk = chunk.dropna(how='any', subset=['EC50 (nM)'])
        cleaned_frame = pd.concat([cleaned_frame, chunk])


    cleaned_frame.to_csv(output_file, sep='\t')

    if not save_output:
        os.remove(output_file)

# PART 2
# move the other file here