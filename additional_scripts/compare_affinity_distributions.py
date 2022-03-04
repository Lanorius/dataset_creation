import sys
import numpy as np
import scipy.stats as stats

data_a = np.loadtxt("p1_2_list_of_affinities_before.txt")
data_b = np.loadtxt("p4_2_list_of_affinities_after.txt")

out = stats.ttest_ind(data_a, data_b, equal_var=False)

print(np.var(data_a))
print(np.var(data_b))
print(out)
