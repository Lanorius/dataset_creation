import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np

path = '/media/lanorius/Elements/data/BindingDB_All.tsv'
output_species = '/media/lanorius/Elements/data/statistics_on_bDB/unique_species.csv'
sep = '\t'

output_Ki = '/media/lanorius/Elements/data/statistics_on_bDB/bindingDB_Ki.pdf'
output_IC50 = '/media/lanorius/Elements/data/statistics_on_bDB/bindingDB_IC50.pdf'
output_pKd = '/media/lanorius/Elements/data/statistics_on_bDB/bindingDB_pKd.pdf'
output_EC50 = '/media/lanorius/Elements/data/statistics_on_bDB/bindingDB_EC50.pdf'
output_kon = '/media/lanorius/Elements/data/statistics_on_bDB/bindingDB_kon.pdf'
output_koff = '/media/lanorius/Elements/data/statistics_on_bDB/bindingDB_koff.pdf'

'''
# short stats on the species
species = ['Target Source Organism According to Curator or DataSource']
unique_species = pd.read_csv(filepath_or_buffer=path, sep=sep, usecols=species)

freq_table = unique_species[species].squeeze().value_counts().to_frame()

row_num = unique_species.shape[0]
freq_table['percetage'] = freq_table[species] / row_num

freq_table.to_csv(output_species, sep='\t')
'''


# stats on interaction values
interaction_values = ['Ki (nM)', 'IC50 (nM)', 'Kd (nM)', 'EC50 (nM)', 'kon (M-1-s-1)', 'koff (s-1)']
interaction_tsv = pd.read_csv(filepath_or_buffer=path, sep=sep, usecols=interaction_values)
interaction_tsv['Ki (nM)'] = pd.to_numeric(interaction_tsv['Ki (nM)'], errors='coerce')
interaction_tsv['IC50 (nM)'] = pd.to_numeric(interaction_tsv['IC50 (nM)'], errors='coerce')
interaction_tsv['Kd (nM)'] = pd.to_numeric(interaction_tsv['Kd (nM)'], errors='coerce')
interaction_tsv['EC50 (nM)'] = pd.to_numeric(interaction_tsv['EC50 (nM)'], errors='coerce')
interaction_tsv['kon (M-1-s-1)'] = pd.to_numeric(interaction_tsv['kon (M-1-s-1)'], errors='coerce')
interaction_tsv['koff (s-1)'] = pd.to_numeric(interaction_tsv['koff (s-1)'], errors='coerce')

interaction_tsv = interaction_tsv.apply(lambda x : np.log10(x) if x.name == 'IC50 (nM)' else x)
interaction_tsv = interaction_tsv.apply(lambda x : -1*np.log10(x/10**9) if x.name == 'Kd (nM)' else x)
interaction_tsv = interaction_tsv.apply(lambda x : np.log10(x) if x.name == 'EC50 (nM)' else x)

interaction_tsv.replace([np.inf, -np.inf], np.nan, inplace=True)


print('\n'+'Ki')

interaction_tsv_Ki = interaction_tsv.dropna(how='any', subset=['Ki (nM)'])
interaction_tsv_Ki['Ki (nM)'].values[interaction_tsv_Ki['Ki (nM)'].values > 10] = 10
print(interaction_tsv_Ki.shape)
print(interaction_tsv_Ki['Ki (nM)'].min())
print(interaction_tsv_Ki['Ki (nM)'].max())

ax_Ki = interaction_tsv_Ki['Ki (nM)'].plot(kind='hist', bins=20, title="Frequency vs Ki")
ax_Ki.set_xlabel("Ki")  # plot has an xlabel argument, but for some reason that throws an error
ax_Ki.get_figure().savefig(output_Ki)
ax_Ki.remove()  # otherwise the artist is kept and the next plot will still include this plot


print('\n'+'IC50')

interaction_tsv_IC50 = interaction_tsv.dropna(how='any', subset=['IC50 (nM)'])
#interaction_tsv_IC50['IC50 (nM)'].values[interaction_tsv_Ki['IC50 (nM)'].values > 10] = 10
print(interaction_tsv_IC50.shape)
print(interaction_tsv_IC50['IC50 (nM)'].min())
print(interaction_tsv_IC50['IC50 (nM)'].max())

ax_IC50 = interaction_tsv_IC50['IC50 (nM)'].plot(kind='hist', bins=20, title="Frequency vs log10(IC50)")
ax_IC50.set_xlabel("log10(IC50)")  # plot has an xlabel argument, but for some reason that throws an error
ax_IC50.get_figure().savefig(output_IC50)
ax_IC50.remove()  # otherwise the artist is kept and the next plot will still include this plot

print('\n'+'pKd')

interaction_tsv_pKd = interaction_tsv.dropna(how='any', subset=['Kd (nM)'])
print(interaction_tsv_pKd.shape)
print(interaction_tsv_pKd['Kd (nM)'].min())
print(interaction_tsv_pKd['Kd (nM)'].max())

ax_pKd = interaction_tsv_pKd['Kd (nM)'].plot(kind='hist', bins=20, title="Frequency vs pKd")
ax_pKd.set_xlabel("pKd")  # plot has an xlabel argument, but for some reason that throws an error
ax_pKd.get_figure().savefig(output_pKd)
ax_pKd.remove()


print('\n'+'EC50')

interaction_tsv_EC50 = interaction_tsv.dropna(how='any', subset=['EC50 (nM)'])
print(interaction_tsv_EC50.shape)
print(interaction_tsv_EC50['EC50 (nM)'].min())
print(interaction_tsv_EC50['EC50 (nM)'].max())

ax_EC50 = interaction_tsv_EC50['EC50 (nM)'].plot(kind='hist', bins=20, title="Frequency vs log10(EC50)")
ax_EC50.set_xlabel("log10(EC50)")  # plot has an xlabel argument, but for some reason that throws an error
ax_EC50.get_figure().savefig(output_EC50)
ax_EC50.remove()  # otherwise the artist is kept and the next plot will still include this plot