import pandas as pd

from ast import literal_eval
from tqdm import tqdm


def drop_unwanted_troublemakers(col_names, row_names, files, output, params):
    """
    :param col_names: column names returned from make_dict_cd_hit
    :param row_names: row names returned from make_dict_mayachemtools
    :param files: for the file path as specified in the config
    :param output: intermediate output or final output file names, names specified in the config
    :param params: integer for maximum drug smile length, and characters that shouldn't appear in the drugs
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

    # Second, removing compounds that are either too long, or have unwanted characters. This step is optional.
    # Removing compounds that are too long.
    if params['drug_length'] != "":
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

    intermediate_drugs.to_csv(files['path'] + output['intermediate_drug_representatives'], sep='\t')
    # TODO: You need to create a proper clustered drugs output file, currently the intermediate is the finished one
    interaction_file.to_csv(files['path'] + output['intermediate_interaction_file'], sep='\t')

    return frame_a, frame_b, compounds_appearing_more_than_once
