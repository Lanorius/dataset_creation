import pandas as pd
import numpy as np


# only checking if I have enough scores to compute Kiba scores from bindingDB 

cols = ['Ki (nM)','IC50 (nM)','Kd (nM)']
path = '/media/daniel/Elements/data/BindingDB_All_2021m7.tsv/BindingDB_All.tsv'

bindingdb_all = pd.read_csv(filepath_or_buffer=path, sep='\t', usecols = cols, error_bad_lines=False)


bindingdb_all = bindingdb_all.dropna( how='any', subset=['IC50 (nM)'])

bindingdb_all = bindingdb_all.dropna(subset=['Ki (nM)', 'Kd (nM)'], thresh=1)

print(bindingdb_all)




