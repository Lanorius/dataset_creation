import pandas as pd
from silx.io.dictdump import h5todict


raw_data = h5todict("p4_2_boxplot_dict.h5")

while len(raw_data) == 1:
	raw_data = raw_data[list(raw_data.keys())[0]]
sum_data = 0

for i in raw_data:
	sum_data += len(raw_data[i])
	
print(sum_data/len(raw_data))
