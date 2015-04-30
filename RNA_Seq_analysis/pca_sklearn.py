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
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, set_link_color_palette, to_tree, inconsistent
import seaborn as sns
from matplotlib.colors import rgb2hex, colorConverter
from collections import defaultdict

path_to_file ='/Volumes/Seq_data/Pdgfra2_all_fpkm_analysis'
run_pca =True


fpbcell = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_by_cell.p'), 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_cell_list.p'), 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_by_gene.p'), 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_gene_list.p'), 'rb')
gene_list = pickle.load(fpgenelist)


df_by_gene = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
df_by_cell = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)
cols = df_by_cell.shape[0]
rows = df_by_cell.shape[1]

row_dist = pd.DataFrame(squareform(pdist(df_by_cell, metric='euclidean')), columns=gene_list, index=gene_list)
row_clusters = linkage(row_dist, metric='euclidean', method='ward')
link_mat = pd.DataFrame(row_clusters,
             columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
             index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])
row_dendr = dendrogram(row_clusters, labels=cell_list, leaf_rotation=90, leaf_font_size=8);
plt.clf()
cluster_idxs = defaultdict(list)
for c, pi in zip(row_dendr['color_list'], row_dendr['icoord']):
    for leg in pi[1:3]:
        i = (leg - 5.0) / 10.0
        if abs(i - int(i)) < 1e-5:
            cluster_idxs[c].append(int(i))

class Clusters(dict):
    def _repr_html_(self):
        html = '<table style="border: 0;">'
        for c in self:
            hx = rgb2hex(colorConverter.to_rgb(c))
            html += '<tr style="border: 0;">' \
            '<td style="background-color: {0}; ' \
                       'border: 0;">' \
            '<code style="background-color: {0};">'.format(hx)
            html += c + '</code></td>'
            html += '<td style="border: 0"><code>'
            html += repr(self[c]) + '</code>'
            html += '</td></tr>'

        html += '</table>'

        return html

cluster_classes = Clusters()
for c, l in cluster_idxs.items():
    i_l = [row_dendr['ivl'][i] for i in l]
    for l2 in i_l:
        cluster_classes[l2] = c

clf = skPCA(2)
np_by_gene = np.asarray(by_gene)

if run_pca:
    by_gene_trans = clf.fit_transform(np_by_gene)
    with open(os.path.join(path_to_file,'gene_fit_trans.p'), 'wb') as fp:
        pickle.dump(by_gene_trans, fp)
else:
    p_gen_trans =  open(os.path.join(path_to_file,'gene_fit_trans.p'), 'rb')
    by_gene_trans = pickle.load(p_gen_trans)
    p_gen_trans.close()
print by_gene_trans
plt.scatter(by_gene_trans[:, 0], by_gene_trans[:, 1], c=[cluster_classes[x] for x in gene_list], edgecolor='none', alpha=0.5)

plt.savefig('skpca_1.png', bbox_inches='tight')
