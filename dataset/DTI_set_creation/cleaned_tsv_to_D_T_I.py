"""
This file takes a cleaned DTI tsv created by input_to_clean_tsv.py and divides it into a .fasta file for the proteins,
a file with ligands, and an interaction file.

"""

# will be part of the config soon

file = pd.read_csv(filepath_or_buffer=path, sep=sep, chunksize=chunksize, usecols=cols, on_bad_lines='skip', engine='python'):
