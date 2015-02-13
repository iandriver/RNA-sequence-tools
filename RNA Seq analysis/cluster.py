import cPickle as pickle
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

path_to_file ='/Volumes/Seq_data/results_pdgfra1_ctrl_pnx'

fpbcell = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_by_cell.p'), 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_cell_list.p'), 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_by_gene.p'), 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_gene_list.p'), 'rb')
gene_list = pickle.load(fpgenelist)

print by_cell.shape
df_by_cell = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
print df_by_cell.shape
row_dist = pdist(np.array(by_cell))
row_clusters = linkage(row_dist)
dendrogram(row_clusters, labels=cell_list)
plt.show()