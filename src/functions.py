import pandas as pd
import math  # to calculate better chunk sizes
import os  # to remove intermediate files
import subprocess  # to run CD-Hit and mayachemtools
from tqdm import tqdm  # shows progress of a few loops
import traceback  # needed in update_interactions
import numpy as np
from silx.io.dictdump import dicttoh5  # to save h5 files
from ast import literal_eval


def raw_transformer(files, file_specifications, output):
    # takes a data_set, like the one from BindingDB and does a basic cleanup to it
    cols = [file_specifications['protein_IDs'], file_specifications['ligand_IDs'],
            file_specifications['protein_sequence'], file_specifications['ligand_SMILE'],
            file_specifications['interaction_value']]

    cleaned_frame = pd.DataFrame(columns=cols)

    chunksize = math.ceil(len(list(open(files['raw_file']))) / 5)  # allows handeling large sets of input data
    for chunk in pd.read_csv(filepath_or_buffer=files['raw_file'], sep=file_specifications['separator'],
                             chunksize=chunksize, usecols=cols, on_bad_lines='skip', engine='python'):
        # chunk = chunk[chunk['Target Source Organism According to Curator or DataSource'] == "Homo sapiens"]
        # removed from the config for now, since we want a general understanding of binding
        chunk = chunk.dropna(how='any', subset=[file_specifications['protein_IDs']])
        chunk = chunk.dropna(how='any', subset=[file_specifications['ligand_IDs']])
        chunk[file_specifications['interaction_value']] = pd.to_numeric(chunk[file_specifications['interaction_value']],
                                                                        errors='coerce')
        chunk = chunk.dropna(how='any', subset=[file_specifications['interaction_value']])
        cleaned_frame = pd.concat([cleaned_frame, chunk])
    cleaned_frame.to_csv(files['path']+output['cleaned_frame'], sep='\t')

    return 0


def create_raw_files(files, file_specifications, output):
    f = open(files['path'] + output['target_file'], 'w')
    d = open('temp_drugs.txt', 'w')
    file = pd.read_csv(filepath_or_buffer=(files['path'] + output['cleaned_frame']), sep='\t', engine='python')

    for index, row in file.iterrows():

        f.write(">" + row[file_specifications['protein_IDs']] + "\n")
        i = 0
        while i < len(row[file_specifications['protein_sequence']]):
            if i % 40 != 39:
                f.write(row[file_specifications['protein_sequence']][i])
                i += 1
            else:
                f.write(row[file_specifications['protein_sequence']][i])
                f.write("\n")
                i += 1
        if i % 40 != 0:
            f.write("\n")

        d.write(row[file_specifications['ligand_SMILE']] + " " + row[file_specifications['ligand_IDs']] + "\n")

    # ensures the output is a csv file
    # all duplicate lines are removed here from the drugs as well
    temp_drugs = pd.read_csv('temp_drugs.txt', sep=' ').drop_duplicates(keep='first').reset_index()
    temp_drugs = temp_drugs.iloc[:, 1:]
    temp_drugs.to_csv(files['path'] + output['drug_file'], sep=' ', index=False)
    os.remove('temp_drugs.txt')

    interactions = file.pivot_table(index=file_specifications['ligand_IDs'], columns=file_specifications['protein_IDs'],
                                    values=file_specifications['interaction_value'], aggfunc='sum')
    interactions.to_csv(files['path'] + output['interaction_file'], sep='\t')

    f.close()
    d.close()

    return 0


def cluster_drugs(files, output, params):
    # RDKit offers the tools necessary to cluster SMILES.
    # For this part you need mayachemtools which uses RDKit and you can find it here:
    # http://www.mayachemtools.org/docs/scripts/html/index.html

    clustering_process = params['mayachemtools_path'] + ' --butinaSimilarityCutoff ' + \
                         params['smile_similarity'] + ' --butinaReordering=yes ' + \
                         '-i ' + files['path'] + output['drug_file'] + \
                         ' -o ' + files['path'] + output['clustered_drugs']
    subprocess.call(clustering_process, shell=True)

    return 0


def cluster_targets(files, output, params):
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

    outside_python = "cd-hit -i " + files['path'] + output['target_file'] + " -o " + \
                     files['path'] + output['target_cluster'] + \
                     " -c " + params['sequence_similarity'] + " -n " + str(word_size)
    subprocess.run(outside_python, shell=True)
    outside_python = "clstr2txt.pl " + files['path'] + output['target_cluster'] + ".clstr > " + \
                     files['path'] + output['target_representatives']
    subprocess.run(outside_python, shell=True)

    # removing duplicate rows from the target_representative file
    temp_target_reps = pd.read_csv(files['path'] + output['target_representatives'], sep='\t')
    temp_target_reps = temp_target_reps.drop_duplicates(subset='id', keep="first")
    temp_target_reps.to_csv(files['path'] + output['target_representatives'], sep='\t', index=False)

    return 0


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


def kd_to_pkd(data):
    data = pd.read_csv(data, sep='\t', header=0, index_col=0).apply(lambda x: -np.log10(x / 1e9))
    return data


def drop_unwanted_troublemakers(col_names, row_names, files, output, params):
    frame_a = pd.DataFrame(0.0, columns=col_names, index=row_names, dtype=float)
    frame_b = pd.DataFrame(0.0, columns=col_names, index=row_names, dtype=float)

    # First, removing compounds that are tautomeres, but RDKit didn't cluster them properly.
    compounds_appearing_more_than_once = []
    for i, _ in frame_a.iterrows():
        if type(frame_a.at[i, frame_a.columns[0]]) == pd.core.series.Series:
            compounds_appearing_more_than_once += [i]
    compounds_appearing_more_than_once = list(set(compounds_appearing_more_than_once))
    intermediate_drugs = pd.read_csv(files['path'] + output['clustered_drugs'], sep=',', header=0,
                                     index_col=1)
    interaction_file = pd.read_csv(files['path'] + output['interaction_file'], sep='\t', header=0, index_col=0)

    frame_a = frame_a.drop(compounds_appearing_more_than_once)
    frame_b = frame_b.drop(compounds_appearing_more_than_once)
    intermediate_drugs = intermediate_drugs.drop(compounds_appearing_more_than_once)
    interaction_file = interaction_file.drop(compounds_appearing_more_than_once)

    # Second, removing compounds that are either too long, or have unwanted characters. This step is optional.
    # Removing compounds that are too long.
    if len(params['drug_length']) > 0:
        for i in tqdm(range(intermediate_drugs.shape[0] - 1, -1, -1)):
            if len(intermediate_drugs.iat[i, 0]) > int(params['drug_length']):
                interaction_file = interaction_file.drop(index=intermediate_drugs.index[i])
                intermediate_drugs = intermediate_drugs.drop(index=intermediate_drugs.index[i])

    # Removing compounds that have a character that can't be processed
    if len(literal_eval(params['bad_characters'])) > 0:
        for i in tqdm(range(intermediate_drugs.shape[0] - 1, -1, -1)):
            if len([char for char in literal_eval(params['bad_characters'])
                    if (char in intermediate_drugs.iat[i, 0])]) > 0:
                interaction_file = interaction_file.drop(index=intermediate_drugs.index[i])
                intermediate_drugs = intermediate_drugs.drop(index=intermediate_drugs.index[i])

    # TODO: The code works fine even though there are some compounds getting lost on the way, check if you have time.
    '''
    print(list(set(interaction_file.index)-set(intermediate_drugs.index)))
    print(list(set(intermediate_drugs.index)-set(interaction_file.index)))

    print(len(list(set(interaction_file.index))))
    print(len(list(set(intermediate_drugs.index))))

    # Compounds which should be in the drug_file but got lost
    lost_compounds = list(set(interaction_file.index)-set(intermediate_drugs.index))
    interaction_file = interaction_file.drop(lost_compounds)

    compounds_lost = list(set(intermediate_drugs.index)-set(interaction_file.index))
    intermediate_drugs = intermediate_drugs.drop(compounds_lost)
    '''

    intermediate_drugs.to_csv(files['path'] + output['intermediate_drug_representatives'], sep='\t')
    interaction_file.to_csv(files['path'] + output['intermediate_interaction_file'], sep='\t')

    return frame_a, frame_b, compounds_appearing_more_than_once


'''
    compounds_appearing_more_than_once = []
    for i, _ in df_a.iterrows():
        if type(df_a.at[i, df_a.columns[0]]) == pd.core.series.Series:
            compounds_appearing_more_than_once += [i]
    compounds_appearing_more_than_once = list(set(compounds_appearing_more_than_once))

    intermediate_drugs = pd.read_csv(files['path'] + output['intermediate_drug_representatives'], sep=',', header=0,
                                     index_col=1)
    interaction_file = pd.read_csv(files['path'] + file['interaction_file'], sep='\t', header=0, index_col=0)

    df_a = df_a.drop(compounds_appearing_more_than_once)
    df_b = df_b.drop(compounds_appearing_more_than_once)
    intermediate_drugs = intermediate_drugs.drop(compounds_appearing_more_than_once)
    interaction_file = interaction_file.drop(compounds_appearing_more_than_once)
'''


def update_interactions(data, frame_a, frame_b, dict_of_drugs, dict_of_targets, files, output):
    key_errors = []
    box_plot_dict = {}  # to create some visualizations of the data
    print('Done by: ' + str(data.shape[1]))
    for name, _ in tqdm(data.iteritems()):
        for index, _ in data.iterrows():
            if data.at[index, name] > 0:
                try:
                    frame_a.at[dict_of_drugs[index], dict_of_targets[name]] += data.at[index, name]
                    frame_b.at[dict_of_drugs[index], dict_of_targets[name]] += 1
                    box_key = str(index)+'_'+str(name)
                    if box_key in box_plot_dict:
                        box_plot_dict[box_key] += [data.at[index, name]]
                    else:
                        box_plot_dict[box_key] = [data.at[index, name]]
                except Exception:
                    error_msg = traceback.format_exc()
                    key_errors += [error_msg.split('\n')[-2][10:]]  # saves faulty keys
    print('Done by: ' + str(frame_a.shape[0]))
    for name, _ in tqdm(frame_a.iteritems()):
        for index, _ in frame_a.iterrows():
            if frame_a.at[index, name] != 0:
                frame_a.at[index, name] = frame_a.at[index, name] / frame_b.at[index, name]
            else:
                frame_a.at[index, name] = np.nan
    dicttoh5(box_plot_dict, h5file=output['box_plot_dict'], h5path=files['path'], mode='w', overwrite_data=None,
             create_dataset_args=None, update_mode=None)
    # TODO: Save key_errors to file
    return frame_a


'''
    # intermediate_interactions, key_Errors = update_interactions(file['path'], interaction_file, df_a, df_b, drug_dict,
    #                                                            target_dict)
    intermediate_drugs.to_csv(file['path'] + output['intermediate_drug_representatives'], sep=',')
    intermediate_interactions, key_Errors = update_interactions(interaction_file, df_a, df_b, drug_dict, target_dict)
    intermediate_interactions.to_csv(file['path'] + output['intermediate_interaction_file'], sep='\t')
'''
