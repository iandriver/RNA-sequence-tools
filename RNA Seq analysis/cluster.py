import cPickle as pickle
import numpy as np
import pandas as pd
import os

import matplotlib.pyplot as plt
import scipy
import json
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, set_link_color_palette, to_tree, inconsistent
import seaborn as sns
from matplotlib.colors import rgb2hex, colorConverter
from pprint import pprint

path_to_file ='/Volumes/Seq_data/Pdgfra2_all_fpkm_analysis'

fpbcell = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_by_cell.p'), 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_cell_list.p'), 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_by_gene.p'), 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra2_outlier_gene_list.p'), 'rb')
gene_list = pickle.load(fpgenelist)

sns.set_palette('Set1', 10, 0.65)
palette = sns.color_palette()
set_link_color_palette(map(rgb2hex, palette))
sns.set_style('white')

df_by_cell = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
df_by_gene = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)
row_dist = pd.DataFrame(squareform(pdist(df_by_cell, metric='euclidean')), columns=cell_list, index=cell_list)
row_clusters = linkage(row_dist, metric='euclidean', method='ward')
link_mat = pd.DataFrame(row_clusters,
             columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
             index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])
row_dendr = dendrogram(row_clusters, labels=cell_list, leaf_rotation=90, leaf_font_size=8)



def augmented_dendrogram(*args, **kwargs):

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord'], ):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y >= 200000:
                plt.plot(x, y, 'ro')
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')

    plt.show()

def cluster_indices(cluster_assignments):
    n = cluster_assignments.max()
    indices = []
    for cluster_number in range(1, n + 1):
        indices.append(np.where(cluster_assignments == cluster_number)[0])
    return indices

def clust_members(r_link, cutoff):
    clust = fcluster(r_link,cutoff)
    num_clusters = clust.max()
    indices = cluster_indices(clust)
    return  num_clusters, indices

def print_clust_membs(indices, cell_list):
    for k, ind in enumerate(indices):
        print "cluster", k + 1, "is", [cell_list[x] for x in ind]

def plot_tree(dendr, pos=None, save=False):
    icoord = scipy.array(dendr['icoord'])
    dcoord = scipy.array(dendr['dcoord'])
    color_list = scipy.array(dendr['color_list'])
    xmin, xmax = icoord.min(), icoord.max()
    ymin, ymax = dcoord.min(), dcoord.max()
    if pos:
        icoord = icoord[pos]
        dcoord = dcoord[pos]
    for xs, ys, color in zip(icoord, dcoord, color_list):
        plt.plot(xs, ys, color)
    plt.xlim(xmin-10, xmax + 0.1*abs(xmax))
    plt.ylim(ymin, ymax + 0.1*abs(ymax))
    if save:
        plt.savefig('plot_dendrogram.png')
    plt.show()


# Create a nested dictionary from the ClusterNode's returned by SciPy
def add_node(node, parent):
	# First create the new node and append it to its parent's children
	newNode = dict( node_id=node.id, children=[] )
	parent["children"].append( newNode )

	# Recursively add the current node's children
	if node.left: add_node( node.left, newNode )
	if node.right: add_node( node.right, newNode )


cc = []
# Label each node with the names of each leaf in its subtree
def label_tree(n, id2name):
    # If the node is a leaf, then we have its name
    if len(n["children"]) == 0:
        leafNames = [ id2name[n["node_id"]] ]

    # If not, flatten all the leaves in the node's subtree
    else:
        leafNames = reduce(lambda ls, c: ls + label_tree(c,id2name), n["children"], [])

    cc.append((len(leafNames), [x.strip('\n') for x in leafNames]))
    cc.sort(key=lambda tup: tup[0], reverse = True)

    # Delete the node id since we don't need it anymore and
    # it makes for cleaner JSON
    del n["node_id"]

    # Labeling convention: "-"-separated leaf names
    n["name"] = name = "-".join(sorted(map(str, leafNames)))

    return leafNames

#Makes labeled json tree for visulaization in d3
def make_tree_json(row_clusters, df_by_cell):
    T= to_tree(row_clusters)

    # Create dictionary for labeling nodes by their IDs
    labels = list(df_by_cell.index)
    id2name = dict(zip(range(len(labels)), labels))

    # Initialize nested dictionary for d3, then recursively iterate through tree
    d3Dendro = dict(children=[], name="Root1")
    add_node( T, d3Dendro )
    label_tree( d3Dendro["children"][0], id2name )
    # Output to JSON
    json.dump(d3Dendro, open("d3-dendrogram.json", "w"), sort_keys=True, indent=4)

make_tree_json(row_clusters, df_by_cell)

#makes
def find_twobytwo(cc, threshold_num = 14):
    pair_dict = {}
    parent = cc[0][1]
    p_num = cc[0][0]
    l_nums = [x[0] for x in cc]
    c_lists = [c[1] for c in cc[1:]]
    for i, c in enumerate(c_lists):
        for i2, c2 in enumerate(c_lists):
            if i != i2 and len(c)>threshold_num and len(c2)>threshold_num:
                if c+c2 in c_lists[:max(i,i2)] or c+c2 == parent:
                    pair_dict[len(c)+len(c2)]= [c, c2]
    g_pvalue_dict = {}
    pvalue_by_level_dict = {}
    sig_gene_list = []
    for v, k in pair_dict.items():
        cell_list1 = [x.strip('\n') for x in k[0]]
        cell_list2 = [xx.strip('\n') for xx in k[1]]
        df_by_cell_1 = df_by_gene[cell_list1]
        df_by_cell_2 = df_by_gene[cell_list2]
        df_by_gene_1 = df_by_cell_1.transpose()
        df_by_gene_2 = df_by_cell_2.transpose()
        for g in gene_list:
            g_pvalue = scipy.stats.f_oneway(df_by_gene_1[g], df_by_gene_2[g])
            if g_pvalue[0] > 3 and g_pvalue[1] <= 0.05:
                g_pvalue_dict[g] = g_pvalue
                if g not in [s[0] for s in sig_gene_list]:
                    sig_gene_list.append((g, g_pvalue[1]))
        pvalue_by_level_dict[v] = g_pvalue_dict
    sig_gene_list.sort(key=lambda tup: tup[1])
    sig_just_genes = [sig[0] for sig in sig_gene_list]
    sig_by_cell = df_by_cell[sig_just_genes[:-75]]
    cg = sns.clustermap(sig_by_cell.transpose(), z_score=1)
    plt.show()

#plot_tree(row_dendr)
find_twobytwo(cc)
#augmented_dendrogram(row_clusters, labels=cell_list, leaf_rotation=90, leaf_font_size=8)
