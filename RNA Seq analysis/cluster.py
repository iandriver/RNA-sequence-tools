import cPickle as pickle
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

path_to_file ='/Volumes/Seq_data/Pdgfra2_all_fpkm_analysis'

fpbcell = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_by_cell.p'), 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_cell_list.p'), 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_by_gene.p'), 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_gene_list.p'), 'rb')
gene_list = pickle.load(fpgenelist)

print by_cell.shape
df_by_cell = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
print df_by_cell.shape
row_dist = pd.DataFrame(squareform(pdist(df_by_cell, metric='euclidean')), columns=cell_list, index=cell_list)
row_clusters = linkage(row_dist, method='complete')
link_mat = pd.DataFrame(row_clusters,
             columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
             index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])
row_dendr = dendrogram(row_clusters, labels=cell_list, leaf_rotation=90, leaf_font_size=8)


plt.show()
