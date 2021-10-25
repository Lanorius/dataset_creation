from process_inputs import parse_config
from tqdm import tqdm  # shows progress of for loops
import pandas as pd
import numpy as np
import subprocess  # to run CD-Hit and mayachemtools

import re #bug fixing

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
                         '-i ' + files['drug_file'] + ' -o ' + output['drug_representatives']

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

    # Make_dict creates dictionaries to know which drug/target is the representative.
    # Also creates lists of these representatives which will be used as row and col-names.

    # this works for both drugs and targets because the outputs
    # of Mayachemtools and CD-Hit accidentally use similar columns.
    # If either one of these tools is replaced this function might not work for the output anymore.
    def make_dict(data):
        out_dict = {}
        rows_or_cols = []
        clusternumber = ""
        clusterrep = ""
        for item in tqdm(range(data.shape[0])):
            if clusternumber != data.iat[item, 2]:
                clusternumber = data.iat[item, 2]
                clusterrep = data.iat[item, 1]
                rows_or_cols += [clusterrep]
                out_dict.update({clusterrep: clusterrep})
            else:
                out_dict.update({data.iat[item, 1]: clusterrep})
        return rows_or_cols, out_dict

    row_names, drug_dict = make_dict(pd.read_csv(output['drug_representatives'], sep=',', on_bad_lines='skip'))
    col_names, target_dict = make_dict(pd.read_csv(output['target_representatives'], sep='\t', on_bad_lines='skip'))

    # These two empty frames will be used to create the new averaged clean interaction frame.
    df_a = pd.DataFrame(0.0, columns=col_names, index=row_names, dtype=float)
    df_b = pd.DataFrame(0.0, columns=col_names, index=row_names, dtype=float)
    # you can change that to take min or max instead of the average, or experiment even further

    def update_interactions(data, frame_a, frame_b, dict_of_drugs, dict_of_targets):
        for name, _ in tqdm(data.iteritems()):
            for index, _ in data.iterrows():
                # if not np.isnan(data.at[index, name]):
                if type(data.at[index, name]) == float:
                    frame_a.at[dict_of_drugs[index], dict_of_targets[name]] += data.at[index, name]
                    frame_b.at[dict_of_drugs[index], dict_of_targets[name]] += 1
        frame_a.to_csv('../intermediate_files/frame_a.csv', sep='\t')
        frame_b.to_csv('../intermediate_files/frame_b.csv', sep='\t')
        for name, _ in frame_a.iteritems():
            for index, _ in frame_a.iterrows():
                # print(frame_a.at[index, name])
                if frame_a.at[index, name] != 0:
                    frame_a.at[index, name] = frame_a.at[index, name] / frame_b.at[index, name]
                else:
                    frame_a.at[index, name] = np.nan
        return frame_a

    # frame_a = pd.read_csv('../intermediate_files/frame_a.csv', sep='\t')
    # frame_b = pd.read_csv('../intermediate_files/frame_b.csv', sep='\t')
    i = 0

    # this version was only for finding a bug. The bug has been found but wasn't solved yet
    def update_interactions_(data, frame_a, frame_b, dict_of_drugs, dict_of_targets, i):  # remove the underscore when putting it back in.
        for name, _ in frame_a.iteritems():
            for index, _ in frame_a.iterrows():
                i += 1
                print(index, name, i)
                print(frame_a.at[index, name])
                if frame_a.at[index, name] != 0:
                    frame_a.at[index, name] = frame_a.at[index, name] / frame_b.at[index, name]
                else:
                    frame_a.at[index, name] = np.nan

        return frame_a

    interaction_file = pd.read_csv(files['interaction_file'], sep='\t')
    cleaned_interactions = update_interactions(interaction_file, df_a, df_b, drug_dict, target_dict)
    # cleaned_interactions = update_interactions(interaction_file, df_a, df_b, drug_dict, target_dict, i)
    cleaned_interactions.to_csv(output['cleaned_interaction_file'], sep='\t')

else:
    print('Skipping Drug Target Interaction Update')
