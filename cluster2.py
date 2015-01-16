import cPickle as pickle
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hcluster
import  pylab
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist,squareform
from scipy.cluster import hierarchy

fpbcell = open('spc_mm10_outlier_by_cell.p', 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open('spc_mm10_outlier_cell_list.p', 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open('spc_mm10_outlier_by_gene.p', 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open('spc_mm10_outlier_gene_list.p', 'rb')
gene_list = pickle.load(fpgenelist)
print by_cell.shape
df_by_cell = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
print df_by_cell.shape
row_dist = pd.DataFrame(squareform(pdist(df_by_cell, metric='euclidean')), columns=cell_list, index=cell_list)
row_clusters = linkage(row_dist, method='complete')
row_dendr = dendrogram(row_clusters)

# reorder rows with respect to the clustering
df_rowclust = df_by_cell.ix[row_dendr['leaves']]
pd.DataFrame.plot(df_rowclust)
plt.show()
# makes dendrogram black (1)
hierarchy.set_link_color_palette(['black'])

# plot row dendrogram
fig = plt.figure(figsize=(8,8))
axd = fig.add_axes([0.09,0.1,0.2,0.6])
row_dendr = dendrogram(row_clusters, orientation='right',
                   color_threshold=np.inf) # makes dendrogram black (2))
axd.set_xticks([])
axd.set_yticks([])


# remove axes spines from dendrogram
for i in axd.spines.values():
        i.set_visible(False)


# plot heatmap
axm = fig.add_axes([0.26,0.1,0.6,0.6]) # x-pos, y-pos, width, height
cax = axm.matshow(df_rowclust['Sftpc','Gapdh'], interpolation='nearest', cmap='hot_r')
fig.colorbar(cax)
axm.set_xticklabels([''] + list(df_rowclust.columns))
axm.set_yticklabels([''] + list(df_rowclust.index))

plt.show()
