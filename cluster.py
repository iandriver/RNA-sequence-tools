import cPickle as pickle
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hcluster
import  pylab
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage


fpbcell = open('spc_mm10_outlier_by_cell.p', 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open('spc_mm10_outlier_cell_list.p', 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open('spc_mm10_outlier_by_gene.p', 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open('spc_mm10_outlier_gene_list.p', 'rb')
gene_list = pickle.load(fpgenelist)


def augmented_dendrogram(*args, **kwargs):

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            plt.plot(x, y, 'ro')
            plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')

    return ddata

print cell_list, len(gene_list)

clusters =hcluster.fclusterdata(by_cell, 3, criterion='maxclust', metric='euclidean', depth=1, method='weighted')
print clusters
links = linkage(by_cell, "single")
print links

plt.subplot(1, 2, 1)
show_leaf_counts = False
ddata = augmented_dendrogram(links,
               color_threshold=100000,
               p=6,
               truncate_mode='lastp',
               show_leaf_counts=show_leaf_counts,
               )
plt.title("show_leaf_counts = %s" % show_leaf_counts)

plt.subplot(1, 2, 2)
show_leaf_counts = True
ddata = augmented_dendrogram(links,
               color_threshold=10000,
               p=6,
               truncate_mode='lastp',
               show_leaf_counts=show_leaf_counts,
               )
plt.title("show_leaf_counts = %s" % show_leaf_counts)

plt.show()
