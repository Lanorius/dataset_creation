import pandas as pd

frame = pd.read_csv("final_affinity_matrix.csv", header=0, index_col=0, sep=",")

print(frame.shape)
print(frame.shape[0]*frame.shape[1]-frame.isnull().sum().sum())

