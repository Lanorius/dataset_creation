import traceback

from process_inputs import parse_config
from tqdm import tqdm  # shows progress of for loops
import pandas as pd
import numpy as np
import subprocess  # to run CD-Hit and mayachemtools

# import os
# import matplotlib.pyplot as plt  # for the boxplots

tasks_to_perform, file, output, params = parse_config()

'''
raw input_to_raw_DTI_ready and create_cluster are separate scripts, because I am currently using separate
environments for both for now
this is changing now with the main script
'''

# Part 1 create drug cluster
if tasks_to_perform[0]:
    print('Creating Drug Cluster')

    # RDKit offers the tools necessary to cluster SMILES.
    # For this part you need mayachemtools which uses RDKit and you can find it here:
    # http://www.mayachemtools.org/docs/scripts/html/index.html

    clustering_process = params['mayachemtools_path'] + ' --butinaSimilarityCutoff ' + params['smile_similarity'] + \
        ' --butinaReordering=yes ' + '-i ' + file['path'] + file['drug_file'] + \
        ' -o ' + file['path'] + output['intermediate_drug_representatives']
    subprocess.call(clustering_process, shell=True)

else:
    print('Skipping Drug Cluster')

# Part 2 create target cluster
if tasks_to_perform[1]:
    print('Creating Target Cluster')

    # running CD-hit
    seq_sim = float(params['sequence_similarity'])
    if seq_sim < 0.4:
        raise ValueError('Threshold for sequence similarity needs to be at least 0.4.')

    sim_dict = {0.5: 2, 0.6: 3, 0.7: 4}  # CD-Hit suggests to use these word sizes
    word_size = 5
    for i in sim_dict:
        if i > seq_sim:
            word_size = sim_dict[i]
            break

    outside_python = "cd-hit -i " + file['path'] + file['target_file'] + " -o " + \
                     file['path'] + output['target_cluster'] + \
                     " -c " + params['sequence_similarity'] + " -n " + str(word_size)
    subprocess.run(outside_python, shell=True)
    outside_python = "clstr2txt.pl " + file['path'] + output['target_cluster'] + ".clstr > " + \
                     file['path'] + output['target_representatives']
    subprocess.run(outside_python, shell=True)

    # removing duplicate rows from the target_representative file
    temp_target_reps = pd.read_csv(file['path'] + output['target_representatives'], sep='\t')
    temp_target_reps = temp_target_reps.drop_duplicates(subset='id', keep="first")
    temp_target_reps.to_csv(file['path'] + output['target_representatives'], sep='\t', index=False)
else:
    print('Skipping Target Cluster')

# Part 3 update drug target interactions

if tasks_to_perform[2]:
    print('Updating Drug Target Interactions')

    # Make_dict creates dictionaries to know which drug/target is the representative.
    # Also creates lists of these representatives which will be used as row and col-names.
    # this works for both drugs and targets because the outputs
    # of Mayachemtools and CD-Hit accidentally have a similar column structure.
    # If either one of these tools is replaced this function might not work for the output anymore.
    def make_dict_cd_hit(data):
        cols = []  # targets are always the columns
        out_dict = {}
        clusterrep = "No Target"
        for item in tqdm(range(data.shape[0])):
            if data.iat[item, 4] == 1:
                clusterrep = data.iat[item, 0]
                cols += [clusterrep]
                out_dict.update({clusterrep: clusterrep})
            else:
                out_dict.update({data.iat[item, 0]: clusterrep})

        return cols, out_dict


    def make_dict_mayachemtools(data):
        rows = []  # drugs are always the rows
        out_dict = {}
        last_cluster = data.iat[0, 2]  # first cluster id
        clusterrep = data.iat[0, 1]  # by mayechemtools logic the cluster center comes first
        for item in tqdm(range(data.shape[0])):
            current_cluster = data.iat[item, 2]
            if last_cluster == current_cluster:
                out_dict.update({data.iat[item, 1]: clusterrep})
            else:
                clusterrep = data.iat[item, 1]
                rows += [clusterrep]
                last_cluster = current_cluster
                out_dict.update({clusterrep: clusterrep})

        return rows, out_dict


    # update interactions takes the clusteres and the dictionaries of the cluster ids and
    # averages the interaction values for each cluster by using its members
    def update_interactions(data, frame_a, frame_b, dict_of_drugs, dict_of_targets):
        key_errors = []
        # print('Done by: ' + str(len(dict_of_targets)))
        for name, _ in tqdm(data.iteritems()):
            for index, _ in data.iterrows():
                # if not np.isnan(data.at[index, name]):
                if data.at[index, name] > 0:
                    try:  # there is at least one key error in here. Not sure where it comes from
                        frame_a.at[dict_of_drugs[index], dict_of_targets[name]] += data.at[index, name]
                        frame_b.at[dict_of_drugs[index], dict_of_targets[name]] += 1
                    except Exception:
                        error_msg = traceback.format_exc()
                        key_errors += [error_msg.split('\n')[-2][10:]]  # saves faulty keys
        # frame_a.to_csv('../intermediate_files/frame_a.csv', sep='\t')
        # frame_b.to_csv('../intermediate_files/frame_b.csv', sep='\t')
        for name, _ in tqdm(frame_a.iteritems()):
            for index, _ in frame_a.iterrows():
                if frame_a.at[index, name] != 0:
                    frame_a.at[index, name] = frame_a.at[index, name] / frame_b.at[index, name]
                else:
                    frame_a.at[index, name] = np.nan
        return frame_a, key_errors


    row_names, drug_dict = make_dict_mayachemtools(pd.read_csv(file['path'] +
                                                               output['intermediate_drug_representatives'], sep=','))
    col_names, target_dict = make_dict_cd_hit(pd.read_csv(file['path'] + output['target_representatives'], sep='\t'))

    # These two empty frames will be used to create the new averaged clean interaction frame.
    df_a = pd.DataFrame(0.0, columns=col_names, index=row_names, dtype=float)
    df_b = pd.DataFrame(0.0, columns=col_names, index=row_names, dtype=float)

    # RDKit has an issue with tautomeres. Even at 100% similarity you would expect all instances of the same compound,
    # to cluster, but in some tautomere cases they don't. In these cases pandas changes the datatype of rows with
    # repeating row_names to pd.core.series.Series
    # The worrisome part is, that we don't know if RDKit only fails at clustering tautomeres.
    compounds_appearing_more_than_once = []
    for i, _ in df_a.iterrows():
        if type(df_a.at[i, df_a.columns[0]]) == pd.core.series.Series:
            compounds_appearing_more_than_once += [i]
    compounds_appearing_more_than_once = list(set(compounds_appearing_more_than_once))

    intermediate_drugs = pd.read_csv(file['path'] + output['intermediate_drug_representatives'], sep=',', header=0,
                                     index_col=1)
    interaction_file = pd.read_csv(file['path'] + file['interaction_file'], sep='\t', header=0, index_col=0)

    df_a = df_a.drop(compounds_appearing_more_than_once)
    df_b = df_b.drop(compounds_appearing_more_than_once)
    intermediate_drugs = intermediate_drugs.drop(compounds_appearing_more_than_once)
    interaction_file = interaction_file.drop(compounds_appearing_more_than_once)

    # intermediate_interactions, key_Errors = update_interactions(file['path'], interaction_file, df_a, df_b, drug_dict,
    #                                                            target_dict)
    intermediate_drugs.to_csv(file['path'] + output['intermediate_drug_representatives'], sep=',')
    intermediate_interactions, key_Errors = update_interactions(interaction_file, df_a, df_b, drug_dict, target_dict)
    intermediate_interactions.to_csv(file['path'] + output['intermediate_interaction_file'], sep='\t')

    # Saving faulty indices to a separate file
    if tasks_to_perform[2] and tasks_to_perform[4]:
        lines_to_write = ["Drug ids of tautomeres, that RDKit doesn't put in the same cluster:\n"]
        for tautomere in compounds_appearing_more_than_once:
            lines_to_write += [tautomere + '\n']
        lines_to_write += ["\nDrug and Target ids that are not in any cluster:\n"]
        key_Errors = list(set(key_Errors))
        for key_Error in key_Errors:
            lines_to_write += [key_Error + '\n']
        key_error_file = open(file['path'] + output['key_errors'], 'w')
        key_error_file.writelines(lines_to_write)
        key_error_file.close()

else:
    print('Skipping Drug Target Interaction Update')

if tasks_to_perform[3]:
    print('THIS STEP REMOVES A LOT OF DATA, SINCE ONLY DRUGS WITH A STR LEN OF 120 OR LESS ARE ALLOWED. '
          'Removal of Drugs with characters that ChemVAE can\'t encode.')
    compound_file = pd.read_csv(file['path'] + output['intermediate_drug_representatives'], sep=',')
    interaction_file = pd.read_csv(file['path'] + output['intermediate_interaction_file'], sep='\t', index_col=0)

    print(compound_file.head())

    error_compounds = []
    count_long_ones = 0
    for i in range(compound_file.shape[0] - 1, 0, -1):
        if ((compound_file.iat[i, 1].find('.') != -1) or (compound_file.iat[i, 1].find('e') != -1) or
            (compound_file.iat[i, 1].find('i') != -1)) or len(compound_file.iat[i, 1]) > 120:
            if len(compound_file.iat[i, 1]) > 120:
                count_long_ones += 1
                print(count_long_ones)
            compound_file = compound_file.drop(index=compound_file.index[i])
            try:
                interaction_file = interaction_file.drop(index=compound_file.iat[i, 0])
            except:
                if compound_file.iat[i, 0] not in error_compounds:
                    error_compounds += [compound_file.iat[i, 0]]

    compound_error_file = open(file['path'] + "/chemVAE_bad_compounds.txt", 'w')
    for bad_compound in error_compounds:
        compound_error_file.write(bad_compound + '\n')
    compound_error_file.close()

    compound_file.to_csv(file['path'] + output['drug_representatives'], sep=',', index=False)
    interaction_file.to_csv(file['path'] + output['cleaned_interaction_file'], sep='\t')


# once all the preprocessing is done, the file with the drug cluster representatives can be created
# this step only works with the mayachemtools output
if tasks_to_perform[0]:  # TODO: This has to be moved
    intermediate_drugs = pd.read_csv(file['path'] + output['intermediate_drug_representatives'], sep=',', header=0,
                                     index_col=1)
    clustered_drugs = pd.DataFrame(columns=intermediate_drugs.columns)
    cluster_number = 1
    for i in range(intermediate_drugs.shape[0]):
        if intermediate_drugs.iat[i, 1] == cluster_number:
            clustered_drugs = clustered_drugs.append(intermediate_drugs.iloc[[i]])
            cluster_number += 1
    clustered_drugs.to_csv(file['path'] + output['intermediate_drug_representatives'], sep=',', index=False)

else:
    print('Skipping the removal of Drugs with characters that ChemVAE can\'t encode.')
    if tasks_to_perform[2]:
        cleaned_interactions = pd.read_csv(file['path'] + output['intermediate_interaction_file'], sep='\t')
        cleaned_interactions.to_csv(file['path'] + output['cleaned_interaction_file'], sep='\t')