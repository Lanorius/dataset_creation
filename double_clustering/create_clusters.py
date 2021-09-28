import numpy as np
from process_inputs import parse_config
from rdkit import Chem, DataStructs
import pandas as pd
from tqdm import tqdm

tasks_to_perform, files, output, params = parse_config()


# Part 1 create drug matrix

if tasks_to_perform[0]:
    print('Creating Drug Matrix')
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

    similarity_df.to_csv(path_or_buf=output['drug_matrix'], sep='\t')
else:
    print('Skipping Drug Matrix')


# Part 2 create drug cluster

if tasks_to_perform[1]:
    print('Creating Drug Cluster')

else:
    print('Skipping Drug Cluster')

# Part 3 create target cluster

if tasks_to_perform[2]:
    print('Creating Target Cluster')

else:
    print('Skipping Target Cluster')


# Part 4 update drug target interactions

if tasks_to_perform[3]:
    print('Updating Drug Target Interactions')

else:
    print('Skipping Drug Target Interaction Update')
