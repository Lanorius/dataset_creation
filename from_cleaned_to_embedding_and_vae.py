import pandas as pd


# This file represents the final step before we create embeddings and use the vae 
# It belongs after the Updating Drug Target Interactions step in create_clusters.py

# unclustered_targets.fasta files are the output of CD-Hit with 100% similarity
file = open('intemediate_files/unclustered_targets_Kd.fasta','r')
Lines = file.readlines()

count = 0
good_headers = []
for line in Lines:
	if line.find(">") == 0:
		count+=1
		good_headers+=[line[1:-1]]

#################################################################################

compounds = pd.read_csv('intemediate_files/drugs_Kd.csv', sep=' ')

good_rows = []
for i in range(compounds.shape[0]):
	good_rows+=[compounds.iat[i,1]]

interactions = pd.read_csv('intemediate_files/interaction_Kd.csv', sep='\t', header=0, index_col=0)
interactions = interactions[interactions.columns.intersection(good_headers)] 
# removes proteins that are not in the target file

#################################################################################
	
good_interaction_indices = []
for i in range(interactions.shape[0]-1, 0, -1):
	if interactions.index[i] not in good_rows:
		interactions = interactions.drop(index=interactions.index[i])
	else:
		good_interaction_indices += [interactions.index[i]]
		
good_compound_indices = []
for i in range(compounds.shape[0]-1, 0, -1):
	if compounds.iat[i,1] not in good_interaction_indices:
		compounds = compounds.drop(index=i)
	else:
		good_compound_indices += [compounds.iat[i,1]]

# print(len(good_compound_indices))
# good_compound_indices = list(dict.fromkeys(good_compound_indices))		
# print(len(good_compound_indices))

print(count)
print(compounds.shape)
print(interactions.shape)

interactions.to_csv('intemediate_files/unclustered_interactions_Kd.csv', sep='\t')
compounds.to_csv('intemediate_files/unclustered_compounds_Kd.csv', sep='\t', index=False)



#Chemvae can't encode certain characters this file removes them. First I thought it is just the . character, but some others also don't work. The list might need expending, although I would suggest replacing chemvae.

#This file is currently hard coded, but since it might be useful, I should perhaps generalize it.

compound_file = pd.read_csv('unclustered_compounds_Kd.csv')
interaction_file = pd.read_csv('unclustered_interactions_Kd.csv', sep='\t', index_col=['ChEMBL ID of Ligand'])

for i in range(compound_file.shape[0]-1, 0, -1):
	print(i)
	if (compound_file.iat[i,0].find('.') != -1) or (compound_file.iat[i,0].find('e') != -1) or (compound_file.iat[i,0].find('i') != -1):
		interaction_file = interaction_file.drop(index=compound_file.iat[i,1])
		compound_file = compound_file.drop(index=compound_file.index[i])

		
compound_file.to_csv('unclustered_compounds_Kd_dot_free.csv', index=False)
interaction_file.to_csv('unclustered_interactions_Kd_dot_free.csv', sep='\t')


