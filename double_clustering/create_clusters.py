from process_inputs import parse_config
from rdkit import Chem, DataStructs  # everything fingerprint related
from tqdm import tqdm  # shows progress of for loops
import pandas as pd
import numpy as np
import subprocess  # to run CD-Hit and mayachemtools

tasks_to_perform, files, output, params = parse_config()


# Part 0 create drug matrix | looks like this step isn't needed

if tasks_to_perform[0]:  # this might not be needed since RDKit will do the clustering for you
    print('Creating Drug Matrix')

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
    similarity_df = pd.read_csv(output['drug_matrix'], sep='\t')
    print(similarity_df)
    print(similarity_df.shape)
'''

# Part 1 create drug cluster

if tasks_to_perform[1]:
    print('Creating Drug Cluster')

    # drug_file = pd.read_csv(files['drug_file'], sep=' ', header=None)
    # for this part you need the mayachemtools which you can find here:
    # http://www.mayachemtools.org/docs/scripts/html/index.html

    clustering_process = '../../mayachemtools/bin/RDKitClusterMolecules.py --butinaReordering=yes -i ' \
                         + files['drug_file'] + ' -o ' + output['drug_cluster']

    subprocess.call(clustering_process, shell=True)

else:
    print('Skipping Drug Cluster')

# Part 2 create target cluster

if tasks_to_perform[2]:
    print('Creating Target Cluster')

    # running CD-hit
    outside_python = "cd-hit -i "+files['target_file']+" -o "+output['target_cluster']
    subprocess.run(outside_python, shell=True)
    outside_python = "clstr2txt.pl "+output['target_cluster']+".clstr > " + output['target_representatives']
    subprocess.run(outside_python, shell=True)
else:
    print('Skipping Target Cluster')


# Part 3 update drug target interactions

if tasks_to_perform[2]:
    print('Updating Drug Target Interactions')

else:
    print('Skipping Drug Target Interaction Update')
