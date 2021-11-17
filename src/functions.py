import pandas as pd
import math  # to calculate better chunk sizes
import os  # to remove intermediate files
import subprocess  # to run CD-Hit and mayachemtools


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
                         ' -o ' + files['path'] + output['intermediate_drug_representatives']
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
