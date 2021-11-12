from rdkit import Chem
from mol2vec.features import mol2alt_sentence, MolSentence, DfVec, sentences2vec
from gensim.models import word2vec

import numpy as np
import pandas as pd

import h5py
from silx.io.dictdump import dicttoh5

# this needs to be added to the main script

clustered_drugs = pd.read_csv('clustered_drugs_Kd.csv', sep=',')

def create_rdkit_embeddings(data, smile_column, id_column, output_file='smiles_embeddings_as_hd5.h5'):
	w2v_model = word2vec.Word2Vec.load('model_300dim.pkl')
	encoded_smiles = {}

	data['mol'] = data[smile_column].apply(lambda x: Chem.MolFromSmiles(x))
	data['sentence'] = data.apply(lambda x: MolSentence(mol2alt_sentence(x['mol'], radius=1)), axis=1)
	data['embedding'] = [DfVec(x) for x in sentences2vec(data['sentence'], w2v_model)]

	# print(data['embedding'].head())

	data_mol2vec = np.array([x.vec for x in data['embedding']])
	data_mol2vec = pd.DataFrame(data_mol2vec)

	ids = data[id_column].tolist()
	embeddings = data_mol2vec.values.tolist()


	for i in range(len(ids)):
		encoded_smiles[ids[i]] = np.array(embeddings[i]).astype('float32')

	dicttoh5(encoded_smiles, h5file=output_file, h5path='/', mode='w', overwrite_data=None, 
		create_dataset_args=None, update_mode=None)

create_rdkit_embeddings(clustered_drugs,'SMILES', 'Name')


