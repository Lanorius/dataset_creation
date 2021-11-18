# import pandas as pd
# import ast  # to go from string to list while parsing col-names
# import math  # to calculate better chunk sizes

# import traceback
from process_inputs import parse_config
from src.functions import *
from ast import literal_eval

# import numpy as np

# from tqdm import tqdm  # shows progress of for loops

# import matplotlib.pyplot as plt  # for the boxplots


# TODO: What is the input? Format it the way Kyra did
'''
One idea would be to expect the user to do a very basic level of preprocessing
This document will be much more powerful if we can expect the input to be the following:
    -a csv or tsv file including columns with protein IDs, protein sequences, compound IDs, 
    compound SMILES/Fingerprints, a parameter of affinity (Kd, Ki, IC50, EC50)
    -for each file these columns have to be specified in the config file

add a config file that specifies:
The name of all input files
The wanted parameters: sequence similarity, smile similarity, affinity values...
The name of the output file
Whatever else comes to mind
'''

tasks_to_perform, sub_tasks_to_perform, files, file_specifications, output, params = parse_config()

if not os.path.isdir(files['path']):
    os.mkdir(files['path'])

# PART 1 Raw Transformation
if tasks_to_perform[0]:
    print("Part 1: Performing raw transformation and creating drug, target, and interaction files.")
    raw_transformer(files, file_specifications, output)
    create_raw_files(files, file_specifications, output)
else:
    print("Part 1: Skipping raw transformation, and the creating of drug, target, and interaction files.")

# PART 2 Clustering
if tasks_to_perform[1]:
    print("Part 2: Creating Drug Clusters.")
    cluster_drugs(files, output, params)
else:
    print("Part 2: Skipping Drug Clusters.")

if tasks_to_perform[2]:
    print("Part 2.1: Creating Target Clusters.")
    cluster_targets(files, output, params)
else:
    print("Part 2.2: Skipping Target Clusters.")

if tasks_to_perform[3]:
    print("Part 3: Preparing the files for DTI.")

    row_names, drug_dict = make_dict_mayachemtools(pd.read_csv(files['path'] +
                                                               output['intermediate_drug_representatives'], sep=','))
    col_names, target_dict = make_dict_cd_hit(pd.read_csv(files['path'] + output['target_representatives'], sep='\t'))


else:
    print("Part 3: Skipping the preparation of the files for DTI.")

print(literal_eval(params['bad_characters']))

