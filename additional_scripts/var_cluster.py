import pandas as pd
from silx.io.dictdump import h5todict
import statistics


raw_data = h5todict("p4_2_boxplot_dict.h5")

while len(raw_data) == 1:
	raw_data = raw_data[list(raw_data.keys())[0]]
sum_data = 0

var_of_boxes = []
counter = 0
for i in raw_data:
	if len(raw_data[i]) > 1:
		var_of_boxes += [statistics.variance(raw_data[i])]
	else:
		counter+=1

var_b = [i for i in var_of_boxes if i > 0.5]

#print(counter)
print(len(raw_data)-counter)
print(len(var_b))
print("-----------")
print(statistics.median(var_of_boxes))
print(statistics.mean(var_of_boxes))
print(max(var_of_boxes))

