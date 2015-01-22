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

cov_mat_gene = df_by_gene.cov()
cov_mat_cell = df_by_cell.cov()
corrs_gene_np = np.asarray(cov_mat_gene)
corrs_cell_np = np.asarray(cov_mat_cell)
#pca_gene = PCA(n_components =2).fit_transform(corrs_gene_np)
#pca_cell = PCA(n_components =2).fit_transform(corrs_cell_np)
'''for i in pca_cell:
  plt.scatter(i[0], i[1])'''

all_samples = np.asarray(df_by_gene)
mlab_pca = PCA(all_samples.T)

print('PC axes in terms of the measurement axes'\
        ' scaled by the standard deviations:\n',\
          mlab_pca.Wt)

plt.plot(mlab_pca.Y[0:50,0],mlab_pca.Y[0:50,1], 'o', markersize=7,\
        color='blue', alpha=0.5, label='class1')
plt.plot(mlab_pca.Y[50:100,0], mlab_pca.Y[50:100,1], '^', markersize=7,\
        color='red', alpha=0.5, label='class2')

plt.xlabel('x_values')
plt.ylabel('y_values')
plt.xlim([-4,4])
plt.ylim([-4,4])
plt.legend()
plt.title('Transformed samples with class labels from matplotlib.mlab.PCA()')

plt.show()
