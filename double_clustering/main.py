import numpy as np
from process_inputs import parse_config
from rdkit import Chem, DataStructs
import pandas as pd
from tqdm import tqdm

files, output, params = parse_config()

drug_file = pd.read_csv(files['drug_file'], sep='\t', header=None)

list_of_SMILES = list(dict.fromkeys(list(drug_file[2])))  # removes duplicates
length_of_list = len(list_of_SMILES)

'''
# There is so much redundancy here... 
print(len(list_of_SMILES))
print(len(list(drug_file[2])))
'''

# this is a matrix that will hold the fingerprint similarities of the SMILES
similarity_df = pd.DataFrame(columns=list_of_SMILES, index=list_of_SMILES)

ms = [Chem.MolFromSmiles(x) for x in list_of_SMILES]
fps = [Chem.RDKFingerprint(x) for x in ms]

# filling the matrix
for i in tqdm(range(length_of_list)):
    for j in range(length_of_list):
        if i < j:
            similarity_df.iloc[j, i] = np.nan
        elif i == j:
            similarity_df.iloc[j, i] = 0
        else:
            similarity_df.iloc[j, i] = DataStructs.FingerprintSimilarity(fps[i], fps[j])

similarity_df.to_csv(path_or_buf='SMILE_KD_similarity_matrix.csv', sep='\t')

