import numpy as np
from mlxtend.evaluate import permutation_test
import sys

# http://rasbt.github.io/mlxtend/api_subpackages/mlxtend.evaluate/#permutation_test

def load_data(filepath, fam=""):
    data = []
    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if fam in line:
                data.append(float(line.split()[-1]))
    
    return data

def significance_test(data1, data2):
    p_value = permutation_test(data1, data2,
                           paired=True, # paired
                           method='approximate', 
                           num_rounds=10000, # repetition number
                           func=lambda x, y: np.abs(np.mean(x) - np.mean(y)), # two sided, default
                           seed=0)
    return p_value

# parameters
params = len(sys.argv)
file1 = sys.argv[1]
file2 = sys.argv[2]
fam = ""
if params > 3:
    fam = sys.argv[3]

data1 = load_data(file1, fam)
data2 = load_data(file2, fam)
assert len(data1) == len(data2), (len(data1), len(data2))

if fam:
    print(fam, 'P value: %f' % significance_test(data1, data2))
else:
    print('P value: %f' % significance_test(data1, data2))