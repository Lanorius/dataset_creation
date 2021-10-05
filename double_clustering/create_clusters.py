from process_inputs import parse_config
from rdkit import Chem, DataStructs  # everything fingerprint related
from tqdm import tqdm  # shows progress of for loops
import pandas as pd
import numpy as np
import subprocess  # to run CD-Hit and mayachemtools

tasks_to_perform, files, output, params = parse_config()


# Part 1 create drug cluster

if tasks_to_perform[0]:
    print('Creating Drug Cluster')

    # drug_file = pd.read_csv(files['drug_file'], sep=' ', header=None)
    # for this part you need the mayachemtools which you can find here:
    # http://www.mayachemtools.org/docs/scripts/html/index.html

    clustering_process = '../../mayachemtools/bin/RDKitClusterMolecules.py' + \
                         ' --butinaSimilarityCutoff ' + params['smile_similarity'] + \
                         ' --butinaReordering=yes ' + \
                         '-i ' + files['drug_file'] + ' -o ' + output['drug_cluster']

    subprocess.call(clustering_process, shell=True)

else:
    print('Skipping Drug Cluster')

# Part 2 create target cluster

if tasks_to_perform[1]:
    print('Creating Target Cluster')

    # running CD-hit
    seq_sim = float(params['sequence_similarity'])
    if seq_sim < 0.4:
        print('Threshold for sequence similarity needs to be at least 0.4.')

    sim_dict = {0.5: 2, 0.6: 3, 0.7: 4}  # CD-Hit suggests to use these word sizes
    word_size = 5
    for i in sim_dict:
        if i > seq_sim:
            word_size = sim_dict[i]
            break

    outside_python = "cd-hit -i " + files['target_file'] + " -o "+output['target_cluster'] + \
                     " -c " + params['sequence_similarity'] + " -n " + str(word_size)
    subprocess.run(outside_python, shell=True)
    outside_python = "clstr2txt.pl "+output['target_cluster']+".clstr > " + output['target_representatives']
    subprocess.run(outside_python, shell=True)

    # removing duplicate rows from the target_representative file
    temp_target_reps = pd.read_csv(output['target_representatives'], sep='\t')
    temp_target_reps = temp_target_reps.drop_duplicates(subset='id', keep="first")
    temp_target_reps.to_csv(output['target_representatives'], sep='\t', index=False)
else:
    print('Skipping Target Cluster')


# Part 3 update drug target interactions

if tasks_to_perform[2]:
    print('Updating Drug Target Interactions')


else:
    print('Skipping Drug Target Interaction Update')
