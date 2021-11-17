import pandas as pd
import math  # to calculate better chunk sizes


def raw_transformer(raw_info, specifications, out):
    sep = specifications['separator']
    cols = [specifications['protein_IDs'], specifications['ligand_IDs'],
            specifications['protein_sequence'], specifications['ligand_SMILE'],
            specifications['interaction_value']]

    cleaned_frame = pd.DataFrame(columns=cols)

    chunksize = math.ceil(len(list(open(raw_info['raw_file']))) / 5)  # allows handeling large sets of input data
    for chunk in pd.read_csv(filepath_or_buffer=raw_info['raw_file'], sep=sep, chunksize=chunksize, usecols=cols,
                             on_bad_lines='skip', engine='python'):
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
