from data_loading.process_inputs import parse_config


# add a config file that specifies:
# The name of all input files
# The wanted parameters: sequence similarity, smile similarity, affinity values...
# The name of the output file
# Whatever else comes to mind with time

# Loading the data
name_of_file = parse_config()
print(name_of_file['bindingDB_file'])


# Preprocessing
# If you use multiple sources make sure each Drug-Target pair is unique
# Remove Sequences based on sequence similarity

# Remove SMILES based on (figure this out)

# Clustering?

# You could write an object for each source of data that you will have in the end

# The goal is to have a single dataset in the end that fullfills the given parameters