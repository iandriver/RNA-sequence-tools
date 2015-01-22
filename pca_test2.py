import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
from sklearn.decomposition import PCA as skPCA
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.mlab import PCA
from mpl_toolkits.mplot3d import Axes3D
import scipy.cluster.hierarchy as sch
from sklearn.preprocessing import StandardScaler


def make_pickle(file_to_dump, p_name='gene_annotation.p'):
  with open(p_name, 'wb') as fp:
    pickle.dump(file_to_dump, fp)

def load_pickle(p_name='gene_annotation.p'):
  with open(p_name, 'rb') as fp:
    data = pickle.load(fp)
  return data    
    
path_to_file ='/Users/idriver/RockLab-files'

fpbcell = open(os.path.join(path_to_file,'pdgfra_mm10_outlier_by_cell.p'), 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open(os.path.join(path_to_file,'pdgfra_mm10_outlier_cell_list.p'), 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open(os.path.join(path_to_file,'pdgfra_mm10_outlier_by_gene.p'), 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open(os.path.join(path_to_file,'pdgfra_mm10_outlier_gene_list.p'), 'rb')
gene_list = pickle.load(fpgenelist)

print cell_list, len(gene_list)


df_by_gene = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
df_by_cell = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)
cols = df_by_cell.shape[0]
rows = df_by_cell.shape[1]

X= df_by_gene.ix[:,0:11169].values
y = gene_list

X_std = StandardScaler().fit_transform(X)
cov_mat = np.cov(X_std.T)

run_new = 1

if run_new == 0:
    eig_vals, eig_vecs = np.linalg.eig(cov_mat)
    make_pickle(eig_vals, p_name="eig_vals.p")
    make_pickle(eig_vecs, p_name="eig_vecs.p")
else:
    eig_vals = load_pickle(p_name="eig_vals.p")
    eig_vecs = load_pickle(p_name="eig_vecs.p")



u,s,v = np.linalg.svd(X_std.T)

for ev in eig_vecs:
    #np.testing.assert_array_almost_equal(1.0, np.linalg.norm(ev))
    print np.linalg.norm(ev)
print('Everything ok!')

# Make a list of (eigenvalue, eigenvector) tuples
eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]

# Sort the (eigenvalue, eigenvector) tuples from high to low
eig_pairs.sort(key=lambda tup: tup[0])
eig_pairs.reverse()

# Visually confirm that the list is correctly sorted by decreasing eigenvalues
print('Eigenvalues in descending order:')
for i in eig_pairs[0:10]:
    print(i[0])
