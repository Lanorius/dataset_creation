import pandas as pd
import math  # to calculate better chunk sizes, and remove nans from frequency lists
import os  # to remove intermediate files
import subprocess  # to run CD-Hit and mayachemtools
from tqdm import tqdm  # shows progress of a few loops
import traceback  # needed in update_interactions
import numpy as np
from silx.io.dictdump import dicttoh5, h5todict  # to save and load h5 files
from ast import literal_eval
import matplotlib.pyplot as plt
import time  # not needed for any function, only used to make the output nicer
import seaborn as sns
import random
import statistics


def raw_transformer(files, file_specifications, output, params):
    """
    :param files: input files and path parsed from the config
    :param file_specifications: information about the columns of the input, parsed from the config
    :param output: intermediate output or final output file names, names specified in the config
    :param params: only target_length is needed here
    :return: no value is returned, but a clean_frame file is created which holds only the binding information
    that is needed later
    """

    # takes a data_set, like the one from BindingDB and does a basic cleanup to it
    cols = [file_specifications['protein_IDs'], file_specifications['ligand_IDs'],
            file_specifications['protein_sequence'], file_specifications['ligand_SMILE'],
            file_specifications['interaction_value']]

    cleaned_frame = pd.DataFrame(columns=cols)

    chunksize = math.ceil(len(list(open(files['raw_file']))) / 5)  # allows handeling large sets of input data
    for chunk in pd.read_csv(filepath_or_buffer=files['raw_file'], sep=file_specifications['separator'],
                             chunksize=chunksize, usecols=cols, on_bad_lines='skip', engine='python'):
        chunk = chunk.dropna(how='any', subset=[file_specifications['protein_IDs']])
        chunk = chunk.dropna(how='any', subset=[file_specifications['ligand_IDs']])
        # chunk = chunk.dropna(how='any', subset=[file_specifications['interaction_value']])  # TODO: remove?
        if params['bad_characters'] != "":
            chunk = chunk[~chunk[file_specifications['ligand_SMILE']].str.contains('|'.join(literal_eval(
                params['bad_characters'])))]
        if params['drug_length'] != "":
            chunk = chunk[chunk[file_specifications['ligand_SMILE']].apply(len) <= int(params['drug_length'])]
        if params['target_length'] != "":
            chunk = chunk[chunk[file_specifications['protein_sequence']].apply(len) <= int(params['target_length'])]
        chunk[file_specifications['interaction_value']] = pd.to_numeric(chunk[file_specifications['interaction_value']],
                                                                        errors='coerce')
        chunk = chunk.dropna(how='any', subset=[file_specifications['interaction_value']])
        cleaned_frame = pd.concat([cleaned_frame, chunk])

    cleaned_frame.to_csv(files['path'] + output['cleaned_frame'], sep='\t')

    return 0


def create_raw_files(files, file_specifications, output, kd_pkd=False):
    """
    :param files: input files and path parsed from the config
    :param file_specifications: information about the columns of the input, parsed from the config
    :param output: intermediate output or final output file names, names specified in the config
    :param kd_pkd:  parameter that decides if Kd values should be also transformed into pKd values
    :return: no value is returned, creates the drug, target and interaction files, creates affinity value plot
    """
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

        d.write(row[file_specifications['ligand_SMILE']] + " " + str(row[file_specifications['ligand_IDs']]) + "\n")

    # ensures the output is a csv file
    # all duplicate lines are removed here from the drugs as well
    temp_drugs = pd.read_csv('temp_drugs.txt', sep=' ').drop_duplicates(keep='first').reset_index()
    temp_drugs = temp_drugs.iloc[:, 1:]
    temp_drugs.to_csv(files['path'] + output['drug_file'], sep=' ', index=False)
    os.remove('temp_drugs.txt')

    interactions = file.pivot_table(index=file_specifications['ligand_IDs'], columns=file_specifications['protein_IDs'],
                                    values=file_specifications['interaction_value'], aggfunc='sum')
    if kd_pkd:
        interactions = interactions.apply(lambda x: -np.log10(x / 1e9))
    interactions.to_csv(files['path'] + output['interaction_file'], sep='\t')

    f.close()
    d.close()

    return 0


def save_affinity_values_plot(files, output, before_after, create_plots):
    # TODO: maybe add constant bins
    """
    :param files: input files and path parsed from the config
    :param output: intermediate output or final output file names, names specified in the config
    :param before_after: either "before" or "after" the clustering process
    :param create_plots: True if values have been transformed from Kd to pKd, else False
    :return: no value is returned, creates affinity value plot
    """
    if before_after == "before":
        interactions = pd.read_csv(files['path'] + output['interaction_file'], sep='\t', header=0, index_col=0)
    elif before_after == "after":
        interactions = pd.read_csv(files['path'] + output['cleaned_interaction_file'], sep=',', header=0, index_col=0)
    else:
        raise ValueError("Wrong input: before_after can only either be a string value before or after.")

    flat_interactions = [val for sublist in interactions.values.tolist() for val in sublist]
    flat_interactions = [x for x in flat_interactions if math.isfinite(x)]

    if create_plots:
        plt.hist(flat_interactions, bins=(math.ceil(max(flat_interactions)) + math.ceil(min(flat_interactions))))
        plt.xlabel("Values")
        plt.ylabel("Frequencies")

    if before_after == "before":
        with open(files['path'] + output['binding_affinity_values_before'], "w") as g:
            for s in flat_interactions:
                g.write(str(s) + " ")
                g.write("\n")

        if create_plots:
            plt.title("Interaction Values before clustering.")
            plt.savefig(files['path'] + output['affinity_plot_before_clustering'])
            plt.clf()

    elif before_after == "after":
        with open(files['path'] + output['binding_affinity_values_after'], "w") as g:
            for s in flat_interactions:
                g.write(str(s) + " ")
                g.write("\n")

        if create_plots:
            plt.title("Interaction Values after clustering.")
            plt.savefig(files['path'] + output['affinity_plot_after_clustering'])
            plt.clf()

    else:
        pass

    return 0


def cluster_drugs(files, output, params):
    """
    :param files: for the file path as specified in the config
    :param output: intermediate output or final output file names, names specified in the config
    :param params: location of the mayachemtools folder and drug similarity parameter
    :return: no value is returned, creates the clustered drugs file
    """
    # RDKit offers the tools necessary to cluster SMILES.
    # For this part you need mayachemtools which uses RDKit and you can find it here:
    # http://www.mayachemtools.org/docs/scripts/html/index.html

    clustering_process = params['mayachemtools_path'] + ' --butinaSimilarityCutoff ' + params['smile_similarity'] + \
                         ' --butinaReordering=yes ' + '-i ' + files['path'] + output['drug_file'] + ' -o ' + \
                         files['path'] + output['clustered_drugs']
    subprocess.call(clustering_process, shell=True)

    return 0


def cluster_targets(files, output, params):
    """
    :param files: for the file path as specified in the config
    :param output: intermediate output or final output file names, names specified in the config
    :param params: target similarity parameter
    :return: no value is returned, creates the clustered target file
    """
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


def make_dict_mayachemtools(data):
    """
    :param data: clustered drugs created by mayachemtools
    :return: a list of row names, a dictionary of all drugs as keys and the representatives as values
    """
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


def make_dict_cd_hit(data):
    """
    :param data: target representatives created by CD-Hit
    :return: a list of column names, a dictionary of all targets as keys and the representatives as values
    """
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


def drop_unwanted_troublemakers(col_names, row_names, files, output):
    """
    :param col_names: column names returned from make_dict_cd_hit
    :param row_names: row names returned from make_dict_mayachemtools
    :param files: for the file path as specified in the config
    :param output: intermediate output or final output file names, names specified in the config
    :type: int, list
    :return: two frames required for the creation of the affinity matrix, a list of drugs that cause errors, also
    creates a new interaction file needed in update interactions
    """
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

    intermediate_drugs.to_csv(files['path'] + output['intermediate_drug_representatives'], sep='\t')
    # TODO: You need to create a proper clustered drugs output file, currently the intermediate is the finished one
    interaction_file.to_csv(files['path'] + output['intermediate_interaction_file'], sep='\t')

    return frame_a, frame_b, compounds_appearing_more_than_once


def update_interactions(data, frame_a, frame_b, dict_of_drugs, dict_of_targets, files, output):
    """
    :param data: interaction file created by drop_unwanted_troublemakers
    :param frame_a: frame created by drop_unwanted_troublemakers, holds the sum of the interaction values for
    one cluster
    :param frame_b: frame created by drop_unwanted_troublemakers, holds the number of the interaction values for
    one cluster
    :param dict_of_drugs: dictionary of drugs created by make_dict_mayachemtools
    :param dict_of_targets: dictionary of targets created by make_dict_cd_hit
    :param files: for the file path as specified in the config
    :param output: intermediate output or final output file names, names specified in the config
    :return: drugs and compounds that cause errors for saving, creates the final interaction file
    """

    key_errors = []
    boxplot_dict = {}  # to create some visualizations of the data

    print('Updating Interactions Part 1/2. Done by: ' + str(data.shape[1]))
    time.sleep(1)
    for name, _ in tqdm(data.iteritems()):
        for index, _ in data.iterrows():
            if data.at[index, name] > 0:
                try:
                    frame_a.at[dict_of_drugs[index], dict_of_targets[name]] += data.at[index, name]
                    frame_b.at[dict_of_drugs[index], dict_of_targets[name]] += 1
                    box_key = str(dict_of_drugs[index]) + '_' + str(dict_of_targets[name])
                    if box_key in boxplot_dict:
                        boxplot_dict[box_key] += [data.at[index, name]]
                    else:
                        boxplot_dict[box_key] = [data.at[index, name]]
                except (Exception,):  # Probably not 100% elegant
                    error_msg = traceback.format_exc()
                    key_errors += [error_msg.split('\n')[-2][10:]]  # saves faulty keys

    print('Updating Interactions Part 2/2. Done by: ' + str(frame_a.shape[1]))
    time.sleep(1)
    for name, _ in tqdm(frame_a.iteritems()):
        for index, _ in frame_a.iterrows():
            if frame_a.at[index, name] != 0:
                frame_a.at[index, name] = frame_a.at[index, name] / frame_b.at[index, name]
            else:
                frame_a.at[index, name] = np.nan
    dicttoh5(boxplot_dict, h5file=files['path'] + output['boxplot_dict'], h5path=files['path'], mode='w',
             overwrite_data=None, create_dataset_args=None, update_mode=None)

    frame_a.to_csv(files['path'] + output['cleaned_interaction_file'], sep=',')

    return key_errors


def save_problematic_drugs_targets(compounds_appearing_more_than_once, key_errors, files, output):
    """
    :param compounds_appearing_more_than_once: output of drop_unwanted_troublemakers
    :param key_errors: output of update_interactions
    :param files: for the file path as specified in the config
    :param output: intermediate output or final output file names, names specified in the config
    :return: saves a file with drugs and compounds that caused errors
    """
    lines_to_write = ["Drug ids of tautomeres, that RDKit doesn't put in the same cluster:\n"]
    for tautomere in compounds_appearing_more_than_once:
        lines_to_write += [str(tautomere) + '\n']
    lines_to_write += ["\nDrug and Target ids that are not in any cluster:\n"]
    key_errors = list(set(key_errors))
    for key_error in key_errors:
        lines_to_write += [key_error + '\n']
    key_error_file = open(files['path'] + output['key_errors'], 'w')
    key_error_file.writelines(lines_to_write)
    key_error_file.close()

    return 0


def sample_from_dict(d, sample):  # TODO: obsolete
    keys = random.sample(list(d), sample)
    values = [d[k] for k in keys]
    # return dict(zip(keys, values))
    return keys, values


def boxplot_creator(file, out_file, file_specifications, min_bin_size, sample_size):
    """
    :param file: boxplot data file
    :param out_file: boxplot image file
    :param file_specifications: for compound and protein ids
    :param min_bin_size: min amount of drug target pairs that a cluster needs to have
    :param sample_size: number of randomly drawn drug target clusters/plots that will be created
    :return: saves a boxplot png
    """

    raw_data = h5todict(file)
    # TODO: Maybe find a way to fix this. This step is needed, since silx
    #  saves the dict in a separate dict for every folder down the path the file is saved.
    while len(raw_data) == 1:
        raw_data = raw_data[list(raw_data.keys())[0]]

    dict_with_many_values = {}

    for key in raw_data.keys():
        # Still deciding on which bins have enough data for the boxplot.
        # The integer here decides how many values have to be at least present to be considered for the boxplot.
        if len(raw_data[key]) > min_bin_size:
            dict_with_many_values[key] = raw_data[key]

    keys = random.sample(list(dict_with_many_values), sample_size)
    values = [dict_with_many_values[k] for k in keys]
    # keys, values = sample_from_dict(dict_with_many_values, sample_size)

    means = [statistics.mean(x) for x in values]
    keys = [x for _, x in sorted(zip(means, keys))]
    values = [x for _, x in sorted(zip(means, values))]
    frame = pd.DataFrame(columns=["Keys", "Values", "Interacting"])
    for i in range(len(keys)):
        for value in values[i]:
            if file_specifications['ligand_IDs'] == "PubChem CID":
                key = " ".join(keys[i].split(".0_"))
            else:
                key = keys[i]
            new_row = {'Keys': key, 'Values': value, "Interacting": "Yes" if value < 7 else "No"}
            frame = frame.append(new_row, ignore_index=True)

    sns.boxplot(x='Keys', y='Values', data=frame, color="cornflowerblue")

    sns.stripplot(x='Keys', y='Values', data=frame, linewidth=1, edgecolor="black", hue="Interacting")
    plt.xticks(rotation=90)
    # adding cutoff line
    plt.axhline(y=7, color='r', linestyle='-')
    plt.xlabel(file_specifications['ligand_IDs'] + " and " + file_specifications['protein_IDs'])
    plt.title("Boxplot of pKd values of "+str(sample_size)+" random Clusters")
    plt.tight_layout()
    plt.savefig(out_file)
    plt.clf()

    print(frame)

    return 0
