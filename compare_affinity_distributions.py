import sys
import numpy as np
import scipy.stats as stats

data_a = np.loadtxt("./intermediate_files/Kd_d08_t08/list_of_affinities_before.txt")
data_b = np.loadtxt("./intermediate_files/Kd_d08_t08/list_of_affinities_after.txt")

out = stats.ttest_ind(data_a, data_b, equal_var=True)

print(out)
