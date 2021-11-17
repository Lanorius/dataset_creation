import pandas as pd
import math  # to calculate better chunk sizes
import os  # to remove intermediate files


def raw_transformer(raw_info, specifications, out):
    # takes a data_set, like the one from BindingDB and does a basic cleanup to it
    cols = [specifications['protein_IDs'], specifications['ligand_IDs'],
            specifications['protein_sequence'], specifications['ligand_SMILE'],
            specifications['interaction_value']]

    cleaned_frame = pd.DataFrame(columns=cols)

    chunksize = math.ceil(len(list(open(raw_info['raw_file']))) / 5)  # allows handeling large sets of input data
    for chunk in pd.read_csv(filepath_or_buffer=raw_info['raw_file'], sep=specifications['separator'],
                             chunksize=chunksize, usecols=cols, on_bad_lines='skip', engine='python'):
        # chunk = chunk[chunk['Target Source Organism According to Curator or DataSource'] == "Homo sapiens"]
        # removed from the config for now, since we want a general understanding of binding
        chunk = chunk.dropna(how='any', subset=[specifications['protein_IDs']])
        chunk = chunk.dropna(how='any', subset=[specifications['ligand_IDs']])
        chunk[specifications['interaction_value']] = pd.to_numeric(chunk[specifications['interaction_value']],
                                                                   errors='coerce')
        chunk = chunk.dropna(how='any', subset=[specifications['interaction_value']])
        cleaned_frame = pd.concat([cleaned_frame, chunk])
    cleaned_frame.to_csv(raw_info['path']+out['cleaned_frame'], sep='\t')

    return 0


def create_raw_files(raw_info, specifications, out):
    f = open(raw_info['path'] + out['fasta_file'], 'w')
    d = open('temp_drugs.txt', 'w')
    file = pd.read_csv(filepath_or_buffer=(raw_info['path'] + out['cleaned_frame']), sep='\t', engine='python')

    for index, row in file.iterrows():

        f.write(">" + row[specifications['protein_IDs']] + "\n")
        i = 0
        while i < len(row[specifications['protein_sequence']]):
            if i % 40 != 39:
                f.write(row[specifications['protein_sequence']][i])
                i += 1
            else:
                f.write(row[specifications['protein_sequence']][i])
                f.write("\n")
                i += 1
        if i % 40 != 0:
            f.write("\n")

        d.write(row[specifications['ligand_SMILE']] + " " + row[specifications['ligand_IDs']] + "\n")

    # ensures the output is a csv file
    # all duplicate lines are removed here from the drugs as well
    temp_drugs = pd.read_csv('temp_drugs.txt', sep=' ').drop_duplicates(keep='first').reset_index()
    temp_drugs = temp_drugs.iloc[:, 1:]
    temp_drugs.to_csv(raw_info['path'] + out['drug_file'], sep=' ', index=False)
    os.remove('temp_drugs.txt')

    interactions = file.pivot_table(index=specifications['ligand_IDs'], columns=specifications['protein_IDs'],
                                    values=specifications['interaction_value'], aggfunc='sum')
    interactions.to_csv(raw_info['path'] + out['interaction_file'], sep='\t')

    f.close()
    d.close()

    return 0
