import cPickle as pickle
import numpy as np
import pandas as pd
import os
from subprocess import call
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
import difflib
from operator import itemgetter
import itertools


#base path to pickle files with fpkm or count matrix
path_to_file = '/Volumes/Seq_data/cuffnorm_pdgfra_1_and_2'
#for labeling all output files
base_name = 'pdgfra2_all_n2'

filename = os.path.join(path_to_file, base_name+'subgroups_200_E15unsorted')
call('mkdir -p '+filename, shell=True)

#if you want to restrict the genes inlcuded to a specific genelist, requires 'GeneID' and 'GroupID' header
make_gene_matrix = True
if make_gene_matrix:
    gene_list_file = 'short_pdgfra_list.txt'
#if you want to restrict the cell matrix file to a subset of cells, expects 'SampleID' header
make_cell_matrix = False
if make_cell_matrix:
    cell_file = 'E15.5_unsorted.txt'
    cell_file_source = os.path.join(path_to_file, cell_file)
#choose metric and method for scipy clustering (also used in seaborn clustermap)
metric='euclidean'
method='average'

#load file gene
if path_to_file.split('/')[-1][0:8] == 'cuffnorm':
    by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file,base_name+'_outlier_filtered.txt'), sep='\t')
elif path_to_file.split('/')[-1][0:5] == 'count':
    by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file,base_name+'_normalized_cpm_all.txt'), sep='\t')
by_gene = by_cell.transpose()
#create list of genes
gene_list = by_cell.index.tolist()
#create cell list
cell_list = [x for x in list(by_cell.columns.values)]


df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list, index=cell_list)
df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list)
hu_cc_gene_df = pd.DataFrame.from_csv('/Volumes/Seq_data/cell_cycle_genes.txt', sep='\t')

def cell_cycle(cell_cycle_gene_df, df_by_gene):
    gene_list = df_by_gene.columns.tolist()
    for g_sym, alt_g_name in zip(cell_cycle_gene_df['Symbol'], cell_cycle_gene_df['Gene Name']):
        if g_sym not in gene_list:
            print g_sym
            for g in alt_g_name.split(','):
                if g.strip() in gene_list:
                    cell_cycle_gene_df['Symbol'][g_sym] = g.strip()
                else:
                    cell_cycle_gene_df = cell_cycle_gene_df[cell_cycle_gene_df.Symbol != g_sym]
    cc_gene_df = df_by_gene[cell_cycle_gene_df['Symbol']]
    return cc_gene_df


def make_new_matrix_gene(org_matrix_by_gene, gene_list_file):
    split_on='_'
    gene_df = pd.read_csv(os.path.join(path_to_file, gene_list_file), delimiter= '\t')
    gene_list = gene_df['GeneID'].tolist()
    group_list = gene_df['GroupID'].tolist()
    gmatrix_df = org_matrix_by_gene[gene_list]
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

def make_new_matrix_cell(org_matrix_by_cell, cell_list_file):
    cell_df = pd.read_csv(os.path.join(path_to_file, cell_list_file), delimiter= '\t')
    cell_list_new = [cell.strip('\n') for cell in cell_df['Sample ID'].tolist()]
    cell_list_old = org_matrix_by_cell.columns.values
    cell_list = [c for c in cell_list_new if c in cell_list_old]
    timepoint = cell_df['Timepoint'].tolist()
    cell_type = cell_df['Type'].tolist()
    new_cmatrix_df = org_matrix_by_cell[cell_list]
    new_name_list = ['_'.join([x,y,z]) for x,y,z in zip(cell_list,timepoint,cell_type)]
    new_name_dict = {k:v for k,v in zip(cell_list,new_name_list)}
    print new_name_dict
    new_cmatrix_df = new_cmatrix_df.rename(columns = new_name_dict)
    new_gmatrix_df = new_cmatrix_df.transpose()
    return new_cmatrix_df, new_gmatrix_df



if make_cell_matrix:
    df_by_cell2, df_by_gene2 = make_new_matrix_cell(df_by_cell1, cell_file_source)
else:
    df_by_cell2, df_by_gene2 = df_by_cell1, df_by_gene1

def preprocess_df(np_by_cell, gen_list, number_expressed=3):
    g_todelete = []
    for g1, gene in enumerate(np_by_cell):
        cells_exp = (gene >= 1.0).sum()
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

def log2_oulierfilter(df_by_cell, log2_cutoff = 0.3, plot=False):
    log2_df = np.log2(df_by_cell+1)
    log2_df2= pd.DataFrame(log2_df.convert_objects(convert_numeric=True))
    log_mean = log2_df.mean(axis=0).order(ascending=False)
    log2_sorted = log2_df.reindex_axis(log2_df.mean(axis=0).order(ascending=False).index, axis=1)
    xticks = []
    keep_col= []
    for col, m in zip(log2_sorted.columns.tolist(),log2_sorted.mean()):
        print m
        if m > log2_cutoff:
            keep_col.append(col)
            xticks.append(col+' '+str("%.2f" % m))
    filtered_df_by_cell = df_by_cell[keep_col]
    filtered_df_by_gene = filtered_df_by_cell.transpose()
    filtered_log2 = np.log2(filtered_df_by_cell[filtered_df_by_cell>0])
    if plot:
        ax = sns.boxplot(data=filtered_log2, whis= .75, notch=True)
        ax = sns.stripplot(x=filtered_log2.columns.values, y=filtered_log2.mean(axis=0), size=4, jitter=True, edgecolor="gray")
        xtickNames = plt.setp(ax, xticklabels=xticks)
        plt.setp(xtickNames, rotation=90, fontsize=9)
        plt.show()
        plt.clf()
        sns.distplot(filtered_log2.mean())
        plt.show()
    log2_expdf_cell = np.log2(filtered_df_by_cell+1)
    log2_expdf_gene = log2_expdf_cell.transpose()
    return log2_expdf_cell, log2_expdf_gene

def run_cluster(by_gene_matrix):
    cell_list = [x for x in list(by_gene_matrix.index.values)]
    cell_dist = pdist(np.array(by_gene_matrix), metric='euclidean')
    row_dist = pd.DataFrame(squareform(cell_dist), columns=cell_list, index=cell_list)
    row_clusters = linkage(cell_dist, metric=metric, method='average')
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
def find_twobytwo(cc, df_by_cell, threshold_num = 14):
    gene_list = df_by_cell.index.tolist()
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

def plot_PCA(df_by_gene, num_genes=100, gene_list_filter=False, title='', plot=False):
    gene_list = df_by_gene.columns.tolist()
    print len(gene_list)
    sns.set_palette("RdBu_r", 10, 1)
    if gene_list_filter:
        sig_by_gene = df_by_gene[gene_list_filter]
    else:
        sig_by_gene = df_by_gene
    clf = skPCA(3)
    np_by_gene = np.asarray(sig_by_gene)
    try:
        by_gene_trans = clf.fit_transform(np_by_gene)
        Pc_df = pd.DataFrame(clf.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=sig_by_gene.columns.tolist())
        pca_rank_df = Pc_df.abs().sum(axis=1)
        Pc_sort_df = pca_rank_df.nlargest(len(sig_by_gene.columns.tolist()))
        top_pca_list = Pc_sort_df.index.tolist()
        print top_pca_list[0:num_genes], 'top_pca_list'
        top_by_gene = df_by_gene[top_pca_list[0:num_genes]]
        clf_top = skPCA(n_components=2)
        np_top_cell = np.asarray(top_by_gene.transpose())

        top_gene_trans = clf_top.fit_transform(np_top_cell)
        fig, ax = plt.subplots(figsize=(15, 16))
        ax.scatter(top_gene_trans[:, 0], top_gene_trans[:, 1], alpha=0.75)
        ax.set_xlim([min(top_gene_trans[:, 0])-1, max(top_gene_trans[:, 0]+1)])
        ax.set_ylim([min(top_gene_trans[:, 1])-1, max(top_gene_trans[:, 1]+1)])
        ax.set_title(title)
        print len(top_by_gene.columns), len(top_gene_trans[:, 0]), len(top_gene_trans[:, 1])
        for label, x, y in zip(top_by_gene.columns, top_gene_trans[:, 0], top_gene_trans[:, 1]):
            ax.annotate(label, (x, y))
        if plot:
            plt.show()
        if title != '':
            save_name = '_'.join(title.split(' ')[0:2])
            plt.savefig(os.path.join(filename,save_name+'_skpca.pdf'), bbox_inches='tight')
        else:
            plt.savefig(os.path.join(filename,'non_group_skpca.pdf'), bbox_inches='tight')
        plt.close()
        return top_pca_list
    except:
        print 'SVD did not converge: PCA'
        return []

def clust_heatmap(gene_list, df_by_gene, num_to_plot=len(gene_list), title='', plot=False):
    if num_to_plot >175:
        sns.set(context= 'poster', font_scale = 0.65/(num_to_plot/100))
    else:
        sns.set(context= 'poster', font_scale = 0.37)
    cg = sns.clustermap(df_by_gene[gene_list[0:num_to_plot]].transpose(), metric=metric, method=method, z_score=0, figsize=(25, 18))
    cg.ax_heatmap.set_title(title)
    if plot:
        plt.show()
    cell_linkage = cg.dendrogram_col.linkage
    link_mat = pd.DataFrame(cell_linkage,
                columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
                index=['cluster %d' %(i+1) for i in range(cell_linkage.shape[0])])
    col_order = cg.dendrogram_col.reordered_ind
    if title != '':
        save_name = '_'.join(title.split(' ')[0:2])
        cg.savefig(os.path.join(filename, save_name+'_heatmap.pdf'), bbox_inches='tight')
    else:
        cg.savefig(os.path.join(filename,'Non_group_heatmap.pdf'), bbox_inches='tight')
    plt.close()
    return cell_linkage, df_by_gene[gene_list[0:num_to_plot]], col_order

def make_subclusters(cc, log2_expdf_cell, gene_corr_list=False, fraction_to_plot=6, filename=filename, base_name=base_name):
    parent = cc[0][1]
    p_num = cc[0][0]
    l_nums = [x[0] for x in cc]
    c_lists = [c[1] for c in cc]
    group_ID = 0

    for num_members, cell_list in zip(l_nums, c_lists):
        if num_members < p_num and num_members >= p_num/fraction_to_plot:
            group_ID+=1
            title = 'Group_'+str(group_ID)+'_with_'+str(num_members)+'_cells'
            cell_subset = log2_expdf_cell[cell_list]
            gene_subset = cell_subset.transpose()
            norm_df_cell1 = np.exp2(cell_subset)
            norm_df_cell = norm_df_cell1 -1
            norm_df_cell.to_csv(os.path.join(filename, base_name+'_'+title+'_matrix.txt'), sep = '\t', index_col=0)
            top_pca = plot_PCA(gene_subset, num_genes=gene_number, title=title, plot=False)
            if top_pca != []:
                top_pca_by_gene = gene_subset[top_pca]
                top_pca_by_cell = top_pca_by_gene.transpose()
                if gene_corr_list:
                    top_genes_search = [x for x in top_pca if x not in cc_gene_df.columns.tolist()]
                    corr_plot(gene_corr_list+top_genes_search[0:3], gene_subset, title = title)
                cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, num_to_plot=gene_number, title=title, plot=False)
                plt.close()
            else:
                pass

def clust_stability(log2_expdf_gene, iterations=16):
    sns.set(context='poster', font_scale = 1)
    sns.set_palette("RdBu_r")
    stability_ratio = []
    total_genes = len(log2_expdf_gene.columns.tolist())
    end_num = 1000
    iter_list = range(100,int(round(end_num)),int(round(end_num/iterations)))
    for gene_number in iter_list:
        title= str(gene_number)+' genes plot.'
        top_pca = plot_PCA(log2_expdf_gene, num_genes=gene_number, title=title)
        top_pca_by_gene = log2_expdf_gene[top_pca]
        top_pca_by_cell = top_pca_by_gene.transpose()
        cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, num_to_plot=gene_number, title=title)
        if gene_number == 100:
            s1 = col_order
            s0 = col_order
        else:
            s2= col_order
            sm_running = difflib.SequenceMatcher(None,s1,s2)
            sm_first = difflib.SequenceMatcher(None,s0,s2)
            stability_ratio.append((sm_running.ratio(), sm_first.ratio()))
            s1=col_order
        plt.close()
    x= iter_list[1:]
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    y1= [m[0] for m in stability_ratio]
    y2= [m[1] for m in stability_ratio]
    sns.barplot(x, y1, palette="RdBu_r", ax=ax1)
    ax1.set_ylabel('Running ratio (new/last)')
    sns.barplot(x, y2, palette="RdBu_r", ax=ax2)
    ax2.set_ylabel('Ratio to 100')
    plt.savefig(os.path.join(filename,'clustering_stability.pdf'), bbox_inches='tight')
    plt.show()
    plt.close()
    return stability_ratio

#

gene_number= 200
log2_expdf_cell, log2_expdf_gene = log2_oulierfilter(df_by_cell, plot=False)
gene_df = pd.read_csv(os.path.join(path_to_file, gene_list_file), delimiter= '\t')
gene_list = gene_df['GeneID'].tolist()
#stability_ratio = clust_stability(log2_expdf_gene)
#print stability_ratio
cc_gene_df = cell_cycle(hu_cc_gene_df, log2_expdf_gene)


new_by_gene = log2_expdf_gene[gene_list]
new_by_cell = new_by_gene.transpose()
cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(gene_list, new_by_gene, num_to_plot=len(gene_list))
#cell_dist, row_dist, row_clusters, link_mat, row_dendr = run_cluster(top_pca_by_gene)
#cc = make_tree_json(cell_linkage, plotted_df_by_gene)
#make_subclusters(cc, log2_expdf_cell, gene_corr_list=['Krt7', 'Dpt'])
#sig_gene_list = find_twobytwo(cc, plotted_df_by_gene.transpose())
#augmented_dendrogram(row_clusters, labels=top_pca_by_cell.columns.tolist(), leaf_rotation=90, leaf_font_size=8)
