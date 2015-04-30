import cPickle as pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn import datasets
from sklearn import metrics
from sklearn.tree import DecisionTreeClassifier




fpbcell = open('spc_mm10_outlier_by_cell.p', 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open('spc_mm10_outlier_cell_list.p', 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open('spc_mm10_outlier_by_gene.p', 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open('spc_mm10_outlier_gene_list.p', 'rb')
gene_list = pickle.load(fpgenelist)

X = by_gene

pca = PCA(n_components=2)
X_r = pca.fit_transform(X)

print X_r
# Percentage of variance explained for each components
print('explained variance ratio (first two components): %s'
      % str(pca.explained_variance_ratio_))
plt.figure()
for v in X_r:
  plt.scatter(v[0], v[1])
plt.show()
