import sys
import numpy as np
import scipy.stats as stats

#data_a = np.loadtxt("./intermediate_files/Kd_CID_d1_t1/list_of_affinities_before.txt")
#data_b = np.loadtxt("./intermediate_files/Kd_CID_d1_t1/list_of_affinities_after.txt")
data_a = np.loadtxt(sys.argv[1])
data_b = np.loadtxt(sys.argv[2])

out = stats.ttest_rel(data_a, data_b)

print(out)
