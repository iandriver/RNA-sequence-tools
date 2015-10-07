import cPickle as pickle
import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import scipy
import json
from sklearn.decomposition import PCA as skPCA
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, set_link_color_palette, to_tree, inconsistent
import seaborn as sns
from matplotlib.colors import rgb2hex, colorConverter
from pprint import pprint


#base path to pickle files with fpkm or count matrix
path_to_file = '/Volumes/Seq_data/count-picard_combined_ips17_BU3'
#for labeling all output files
start_file_name = 'ips17_BU3_normERCC'

make_go_matrix = False

#load file gene
by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file,'ips_hu_normalized_Genes_by_ERCC.txt'), sep='\t')
by_gene = by_cell.transpose()
#create list of genes
gene_list = by_cell.index.tolist()
#create cell list
cell_list = [x for x in list(by_cell.columns.values)]


df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list, index=cell_list)
df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list)


def make_new_matrix(org_matrix_by_cell, gene_list_file):
    split_on='_'
    gene_df = pd.read_csv(os.path.join(path_to_file, gene_list_file), delimiter= '\t')
    gene_list = gene_df['GeneID'].tolist()
    group_list = gene_df['GroupID'].tolist()
    gmatrix_df = org_matrix_by_cell[gene_list]
    cmatrix_df = gmatrix_df.transpose()
    cell_list1 = []
    for cell in cmatrix_df.columns.values:
        if exclude:
            if cell.split(split_on)[1] == 'ctrl' or cell.split(split_on)[1] == 'pnx':
                if cell.split(split_on)[2][0] =='C':
                    print cell, 'cell'
                    cell_list1.append(cell)
        else:
            cell_list1.append(cell)
    new_cmatrix_df = cmatrix_df[cell_list1]
    new_gmatrix_df = new_cmatrix_df.transpose()
    return new_cmatrix_df, new_gmatrix_df
if make_go_matrix:
    df_by_cell2, df_by_gene2 = make_new_matrix(df_by_gene1, gene_file_source)
else:
    df_by_cell2, df_by_gene2 = df_by_cell1, df_by_gene1


def preprocess_df(np_by_cell, gen_list, number_expressed=4):
    g_todelete = []
    print gen_list
    for g1, gene in enumerate(np_by_cell):
        cells_exp = (gene >= 1.0).sum()
        print cells_exp
        if cells_exp < number_expressed:
            g_todelete.append(g1)
    g1_todelete = sorted(g_todelete, reverse = True)
    print np_by_cell.shape
    for pos in g1_todelete:
        if type(gen_list[pos]) != float:
            print 'Gene '+gen_list[pos]+' not expressed in '+str(number_expressed)+' cells.'
            pass
        del gen_list[pos]
    n_by_cell = np.delete(np_by_cell, g1_todelete, axis=0)
    print n_by_cell.shape
    return n_by_cell, gen_list

np_by_cell2 = np.array(df_by_cell2.values, dtype='f')
gen_list = df_by_cell2.index.tolist()
np_by_cell, n_gene_list = preprocess_df(np_by_cell2, gen_list)
df_by_gene = pd.DataFrame(np_by_cell.transpose(), index = df_by_cell2.columns.values, columns= n_gene_list)
df_by_cell = df_by_gene.transpose()
log2_df = np.log2(df_by_cell+1)

log2_df2= pd.DataFrame(log2_df.convert_objects(convert_numeric=True))
print log2_df2.columns.values
log_mean = log2_df.mean(axis=0).order(ascending=False)
log2_sorted = log2_df.reindex_axis(log2_df.mean(axis=0).order(ascending=False).index, axis=1)
ax = sns.boxplot(data=log2_sorted, whis= .75, notch=True)
ax = sns.stripplot(x=log2_sorted.columns.values, y=log2_sorted.mean(axis=0), size=4, jitter=True, edgecolor="gray")
xtickNames = plt.setp(ax, xticklabels=log2_sorted.columns.values)
plt.setp(xtickNames, rotation=90, fontsize=9)
plt.show()

def run_cluster(by_gene_matrix):
    cell_dist = pdist(np.array(by_gene_matrix), metric='euclidean')
    row_dist = pd.DataFrame(squareform(cell_dist), columns=cell_list, index=cell_list)
    row_clusters = linkage(row_dist, metric='euclidean', method='ward')
    link_mat = pd.DataFrame(row_clusters,
                 columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
                 index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])
    row_dendr = dendrogram(row_clusters, labels=cell_list, leaf_rotation=90, leaf_font_size=8)

    plt.savefig(os.path.join(path_to_file,'dendrogram_gene.png'))
    plt.clf()
    return cell_dist, row_dist, row_clusters, link_mat, row_dendr

def augmented_dendrogram(*args, **kwargs):
    plt.clf()
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
    plt.savefig(os.path.join(path_to_file,'augmented_dendrogram.png'))

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
        plt.savefig(os.path.join(path_to_file,'plot_dendrogram.png'))
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
def make_tree_json(row_clusters, df_by_gene):
    T= to_tree(row_clusters)

    # Create dictionary for labeling nodes by their IDs
    labels = list(df_by_gene.index)
    id2name = dict(zip(range(len(labels)), labels))

    # Initialize nested dictionary for d3, then recursively iterate through tree
    d3Dendro = dict(children=[], name="Root1")
    add_node( T, d3Dendro )
    label_tree( d3Dendro["children"][0], id2name )
    # Output to JSON
    json.dump(d3Dendro, open(os.path.join(path_to_file,"d3-dendrogram.json"), "w"), sort_keys=True, indent=4)

    return cc


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
        df_by_cell_1 = df_by_cell[cell_list1]
        df_by_cell_2 = df_by_cell[cell_list2]
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
    return sig_just_genes

def plot_PCA(gene_list_filter, df_by_gene):
    sig_by_cell = df_by_gene[gene_list_filter]
    clf = skPCA(2)
    np_by_gene = np.asarray(sig_by_cell)
    by_gene_trans = clf.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(clf.components_.T, columns=['PC-1', 'PC-2'], index=sig_by_cell.columns)
    Pc_1_sort_df = pd.DataFrame.sort(pd.DataFrame.abs(Pc_df), columns = 'PC-1', ascending=False)
    Pc_2_sort_df = pd.DataFrame.sort(pd.DataFrame.abs(Pc_df), columns = 'PC-2', ascending=False)
    pc_1_gene_list = Pc_1_sort_df.index.values
    pc_2_gene_list = Pc_2_sort_df.index.values
    top_pca_list = []
    for i, x in enumerate(pc_1_gene_list):
        if x not in top_pca_list and x[0:2] != 'Rp':
            top_pca_list.append(x)
        elif pc_2_gene_list[i] not in top_pca_list and pc_2_gene_list[i][0:2] != 'Rp':
            top_pca_list.append(pc_2_gene_list[i])
    print top_pca_list[0:150]
    top_by_cell = df_by_gene[top_pca_list[0:150]]
    clf_top = skPCA(n_components=2)
    np_top_gene = np.asarray(top_by_cell.transpose())
    print top_by_cell
    print np_top_gene
    top_gene_trans = clf_top.fit_transform(np_top_gene)
    plt.clf()
    fig, ax = plt.subplots()
    ax.scatter(top_gene_trans[:, 0], top_gene_trans[:, 1], alpha=0.75)
    ax.set_xlim([min(top_gene_trans[:, 0])-1, max(top_gene_trans[:, 0]+1)])
    ax.set_ylim([min(top_gene_trans[:, 1])-1, max(top_gene_trans[:, 1]+1)])
    print np_top_gene.shape
    print top_gene_trans.shape
    print len(top_by_cell.columns), len(top_gene_trans[:, 0]), len(top_gene_trans[:, 1])
    for label, x, y in zip(top_by_cell.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
        ax.annotate(label, (x, y))
    plt.show()
    plt.savefig(os.path.join(path_to_file,'skpca_2.png'), bbox_inches='tight')
    plt.clf()
    sns.set_palette("RdBu_r", 10, 1)
    return top_pca_list

def clust_heatmap(top_pca_list, df_by_gene, num_to_plot=75):
    cg = sns.clustermap(df_by_gene[top_pca_list[0:num_to_plot]].transpose(), z_score=0)
    plt.show()
    plt.savefig(os.path.join(path_to_file,'clustermap_1.png'), bbox_inches='tight')

#plot_tree(row_dendr)
'''cell_dist, row_dist, row_clusters, link_mat, row_dendr = run_cluster(df_by_gene)
cc = make_tree_json(row_clusters, df_by_gene)
sig_gene_list = find_twobytwo(cc)
top_pca = plot_PCA(sig_gene_list, df_by_gene)
clust_heatmap(top_pca, df_by_gene)
augmented_dendrogram(row_clusters, labels=cell_list, leaf_rotation=90, leaf_font_size=8)'''
