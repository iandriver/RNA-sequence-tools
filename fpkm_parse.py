import fnmatch
import os
import cPickle as pickle
import numpy as np
import itertools


with open('spc_fpkm.p', 'rb') as fp:
  data = pickle.load(fp)

num_genes = len(data['genes'])
cell_num = len(data)-1

gene_names = data['genes']
del data['genes']

fpkm_array = np.array([map(float,v) for k, v in data.items()])
cell_array = np.array([k for k, v in data.items()])


by_gene = fpkm_array.transpose()


gene_present_mask = [(np.count_nonzero(x)>5) for x in by_gene]

filtered_genes = np.ma.masked_array(by_gene, mask= gene_present_mask*96)

print cell_array
