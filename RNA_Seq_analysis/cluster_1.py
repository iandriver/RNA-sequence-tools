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
#if you want to use the single cell file (created by make_monocle)
singlecell_file = True

filename = os.path.join(path_to_file, base_name+'_short1_z1')
call('mkdir -p '+filename, shell=True)

#if you want to restrict the genes inlcuded to a specific genelist, requires 'GeneID' and 'GroupID' header
make_gene_matrix = True
if make_gene_matrix:
    gene_list_file = 'short_pdgfra_list1.txt'
#if you want to restrict the cell matrix file to a subset of cells, expects 'SampleID' header
make_cell_matrix = False
if make_cell_matrix:
    cell_file = 'E15.5_unsorted.txt'
    cell_file_source = os.path.join(path_to_file, cell_file)
#choose metric and method for scipy clustering (also used in seaborn clustermap)
metric='euclidean'
method='average'

#load file gene
if singlecell_file:
    by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file,'single_cell_matrix.txt'), sep='\t')
elif path_to_file.split('/')[-1][0:8] == 'cuffnorm':
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
#run correlation matrix and save only those above threshold
def run_corr(df_by_gene, title, method_name='pearson', sig_threshold= 0.35, run_new=True, min_period=3):
    if run_new:
        if method_name != 'kendall':
            corr_by_gene = df_by_gene.corr(method=method_name, min_periods=min_period)
        else:
            corr_by_gene = df_by_gene.corr(method=method_name)
        corr_by_cell = df_by_cell.corr()

        cor = corr_by_gene
        cor.loc[:,:] =  np.tril(cor.values, k=-1)
        cor = cor.stack()
        corr_by_gene_pos = cor[cor >=sig_threshold]
        corr_by_gene_neg = cor[cor <=(sig_threshold*-1)]

        with open(os.path.join(path_to_file,'gene_correlations_sig_neg_'+method_name+'.p'), 'wb') as fp:
            pickle.dump(corr_by_gene_neg, fp)
        with open(os.path.join(path_to_file,'gene_correlations_sig_pos_'+method_name+'.p'), 'wb') as fp0:
            pickle.dump(corr_by_gene_pos, fp0)
        with open(os.path.join(path_to_file,'by_gene_corr.p'), 'wb') as fp1:
            pickle.dump(corr_by_gene, fp1)
        with open(os.path.join(path_to_file,'by_cell_corr.p'), 'wb') as fp2:
            pickle.dump(corr_by_cell, fp2)
    else:
        corr_by_g_pos =  open(os.path.join(path_to_file,'gene_correlations_sig_pos_'+method_name+'.p'), 'rb')
        corr_by_g_neg =  open(os.path.join(path_to_file,'gene_correlations_sig_neg_'+method_name+'.p'), 'rb')
        corr_by_gene_pos = pickle.load(corr_by_g_pos)
        corr_by_gene_neg = pickle.load(corr_by_g_neg)

    cor_pos_df = pd.DataFrame(corr_by_gene_pos)
    cor_neg_df = pd.DataFrame(corr_by_gene_neg)
    sig_corr = cor_pos_df.append(cor_neg_df)
    sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])

    if run_new:
        sig_corrs.to_csv(os.path.join(path_to_file, title+'_counts_corr_sig_'+method_name+'.txt'), sep = '\t')
    return sig_corrs

#corr_plot finds and plots all correlated genes, log turns on log scale, sort plots the genes in the rank order of the gene searched
def corr_plot(terms_to_search, df_by_gene, title, log=False, sort=True, sig_threshold=0.5):
    sig_corrs = run_corr(df_by_gene, title, sig_threshold=sig_threshold)
    for term_to_search in terms_to_search:
        corr_tup = [(term_to_search, 1)]
        neg = True
        fig, ax = plt.subplots()
        marker = itertools.cycle(('+', 'o', '*'))
        linestyles = itertools.cycle(('--', '-.', '-', ':'))
        for index, row in sig_corrs.iterrows():
            if term_to_search in index:
                neg = False
                if index[0]==term_to_search:
                    corr_tup.append((index[1],row['corr']))
                else:
                    corr_tup.append((index[0],row['corr']))

        if neg:
            print term_to_search+' not correlated.'
        corr_tup.sort(key=itemgetter(1), reverse=True)
        corr_df = pd.DataFrame(corr_tup, columns=['GeneID', 'Correlation'])
        corr_df.to_csv(os.path.join(filename, title+'_Corr_w_'+term_to_search+'_list.txt'), sep = '\t', index=False)
        for c in corr_tup:
            print c
        to_plot = [x[0] for x in corr_tup]
        sns.set_palette(sns.cubehelix_palette(len(to_plot), start=1, rot=-.9, reverse=True))
        try:
            sorted_df = df_by_gene.sort([term_to_search])
            log2_df = np.log2(df_by_gene[to_plot])
            sorted_log2_df=np.log2(sorted_df[to_plot])
            ylabel='CPM (log2)'
            if sort and log:
                ax = sorted_log2_df.plot()
                xlabels = sorted_log2_df[to_plot].index.values
            elif sort:
                ax =sorted_df[to_plot].plot()
                xlabels = sorted_df[to_plot].index.values
            elif log:
                ax = log2_df.plot()
                ylabel= 'log2 FPKM'
                xlabels = log2_df.index.values
            else:
                ax = df_by_gene[to_plot].plot()
                xlabels = df_by_gene[to_plot].index.values
            ax.set_xlabel('Cell #')
            ax.set_ylabel(ylabel)
            ax.set_title('Correlates with '+term_to_search, loc='right')
            ax.xaxis.set_minor_locator(LinearLocator(numticks=len(xlabels)))
            ax.set_xticklabels(xlabels, minor=True, rotation='vertical', fontsize=6)
            ax.set_ylim([0, df_by_gene[to_plot].values.max()])
            ax.tick_params(axis='x', labelsize=1)
            if len(corr_tup) > 15:
                l_labels = [str(x[0])+' '+"%.2f" % x[1] for x in corr_tup]
                ax.legend(l_labels, loc='upper left', bbox_to_anchor=(0.01, 1.05), ncol=6, prop={'size':6})
            else:
                l_labels = [str(x[0])+' '+"%.2f" % x[1] for x in corr_tup]
                ax.legend(l_labels, loc='upper left', bbox_to_anchor=(0.01, 1.05), ncol=4, prop={'size':8})
            fig = plt.gcf()
            fig.subplots_adjust(bottom=0.08, top=0.95, right=0.98, left=0.03)
            plt.savefig(os.path.join(filename, title+'_corr_with_'+term_to_search+'.pdf'), bbox_inches='tight')
            plt.close()
        except KeyError:
            print term_to_search+' not in this matrix'
            pass
def single_gene_sig(gene, by_cell_df):
    gene_list= by_cell_df.index.tolist()
    by_gene_df = by_cell_df.transpose()
    sorted_df = by_gene_df.sort([gene])
    cutoffs = sorted_df[gene].quantile([.2,.8])

    low_cut = cutoffs[0.2]
    high_cut = cutoffs[0.8]
    high_g_df = sorted_df[sorted_df[gene]>high_cut]
    mid_g_df = sorted_df[(sorted_df[gene]<high_cut) | (sorted_df[gene]>low_cut)]
    low_g_df = sorted_df[sorted_df[gene]<low_cut]
    high_cells = high_g_df.index.values
    mid_cells = mid_g_df.index.values
    low_cells = low_g_df.index.values
    high_df = by_cell_df[high_cells]
    mid_df = by_cell_df[mid_cells]
    low_df = by_cell_df[low_cells]
    high_gene_df = high_df.transpose()
    mid_gene_df = mid_df.transpose()
    low_gene_df = low_df.transpose()
    print high_cells, mid_cells, low_cells
    g_pvalue_dict = {}
    sig_gene_list = []
    for g in gene_list:
        g_pvalue = scipy.stats.f_oneway(high_gene_df[g], mid_gene_df[g], low_gene_df[g])
        if g_pvalue[0] > 3 and g_pvalue[1] <= 0.05:
            g_pvalue_dict[g] = g_pvalue
            if g not in [s[0] for s in sig_gene_list]:
                sig_gene_list.append([g, g_pvalue[1]])

    sig_gene_list.sort(key=lambda tup: tup[1])
    pvalues = [p[1] for p in sig_gene_list]
    gene_index = [ge[0] for ge in sig_gene_list]
    tcf_mean_1 = high_gene_df['Tcf21'].mean()
    tcf_mean_2 = mid_gene_df['Tcf21'].mean()
    tcf_mean_3 = low_gene_df['Tcf21'].mean()
    print tcf_mean_1, tcf_mean_2, tcf_mean_3
    mean_log2_exp_list = []
    sig_high_low_list = []
    for sig_gene in gene_index:
        sig_gene_df = by_gene_df[sig_gene]
        mean_log2_exp_list.append(sig_gene_df.mean())
        sig_cell_df = sig_gene_df.transpose()
        mean_cell_high = sig_cell_df[high_cells].mean()
        mean_celll_low = sig_cell_df[low_cells].mean()
        ratio_high_low = (mean_cell_high+1)/(mean_celll_low+1)
        sig_high_low_list.append(ratio_high_low)
    sig_df = pd.DataFrame({'pvalues':pvalues,'mean':mean_log2_exp_list,'ratio_hi_low':sig_high_low_list}, index=gene_index)
    sig_df.to_csv(os.path.join(filename,'sig_'+gene+'_pvalues.txt'), sep = '\t')
    sig_gene_good = sig_df[(sig_df['mean']>=1.1) & ((sig_df['ratio_hi_low']>=1.45) | (sig_df['ratio_hi_low']<=.7))]
    print sig_gene_good
    good_genes = sig_gene_good.index.tolist()
    df_to_plot = sorted_df[good_genes]
    ax = df_to_plot.plot()
    xlabels = df_to_plot.index.values
    ax.xaxis.set_minor_locator(LinearLocator(numticks=len(xlabels)))
    ax.set_xticklabels(xlabels, minor=True, rotation='vertical')
    ax.tick_params(axis='x', labelsize=1)
    l_labels = good_genes
    ax.legend(l_labels, loc='upper left', bbox_to_anchor=(0.01, 1.05), ncol=6, prop={'size':6})
    plt.show()

#finds significant genes between subclusters
def find_twobytwo(cc, df_by_cell, full_by_cell_df, fraction_to_plot=8):
    gene_list = full_by_cell_df.index.tolist()
    pair_dict = {}
    parent = cc[0][1]
    p_num = cc[0][0]
    l_nums = [x[0] for x in cc]
    c_lists = [c[1] for c in cc[1:]]
    for i, c in enumerate(c_lists):
        for i2, c2 in enumerate(c_lists):
            if i != i2 and len(c)>=p_num/fraction_to_plot and len(c2)>=p_num/fraction_to_plot:
                if c+c2 in c_lists[:max(i,i2)] or c+c2 == parent:
                    pair_dict[str(len(c))+'cells_vs_'+str(len(c2))+'cells']= [c, c2]

    for v, k in pair_dict.items():
        g_pvalue_dict = {}
        index_list = []
        sig_gene_list = []
        cell_list1 = [x.strip('\n') for x in k[0]]
        cell_list2 = [xx.strip('\n') for xx in k[1]]
        group1 = str(len(cell_list1))
        group2 = str(len(cell_list2))
        df_by_cell_1 = full_by_cell_df[cell_list1]
        df_by_cell_2 = full_by_cell_df[cell_list2]
        df_by_gene_1 = df_by_cell_1.transpose()
        df_by_gene_2 = df_by_cell_2.transpose()
        for g in gene_list:
            g_pvalue = scipy.stats.f_oneway(df_by_gene_1[g], df_by_gene_2[g])
            if g_pvalue[0] > 3 and g_pvalue[1] <= 0.05:
                g_pvalue_dict[g] = g_pvalue
                if g not in [s[0] for s in sig_gene_list]:
                    sig_gene_list.append([g, g_pvalue[1]])

        sig_gene_list.sort(key=lambda tup: tup[1])
        pvalues = [p[1] for p in sig_gene_list]
        gene_index = [ge[0] for ge in sig_gene_list]
        tcf_mean_1 = df_by_gene_1['Tcf21'].mean()
        tcf_mean_2 = df_by_gene_2['Tcf21'].mean()
        print v, gene_index, pvalues, tcf_mean_1, tcf_mean_2
        sig_df = pd.DataFrame(pvalues, index=gene_index)
        sig_df.columns = ['pvalues']
        cell_names_df = pd.DataFrame({'cells1':pd.Series(cell_list1, index=range(len(cell_list1))), 'cells2':pd.Series(cell_list2, index=range(len(cell_list2)))})
        print sig_df
        sig_df.to_csv(os.path.join(filename,'sig_'+v+'_pvalues.txt'), sep = '\t')
        cell_names_df.to_csv(os.path.join(filename,'sig_'+v+'_cells.txt'), sep = '\t')

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
        sns.set(context= 'poster', font_scale = 0.65)
    sns.set_palette('RdBu',4,0.1)
    cg = sns.clustermap(df_by_gene[gene_list[0:num_to_plot]].transpose(), metric=metric, method=method, z_score=1, figsize=(25, 18))
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
        cg.savefig(os.path.join(filename,'Non_group_heatmap_z1_deleted.pdf'), bbox_inches='tight')
    plt.close()
    return cell_linkage, df_by_gene[gene_list[0:num_to_plot]], col_order

def make_subclusters(cc, log2_expdf_cell, gene_corr_list=False, fraction_to_plot=4, filename=filename, base_name=base_name):
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
            if len(top_pca)<gene_number:
                plot_num = len(top_pca)
            else:
                plot_num = gene_number
            if top_pca != []:
                top_pca_by_gene = gene_subset[top_pca]
                top_pca_by_cell = top_pca_by_gene.transpose()
                if gene_corr_list:
                    top_genes_search = [x for x in top_pca if x not in cc_gene_df.columns.tolist()]
                    corr_plot(gene_corr_list+top_genes_search[0:3], gene_subset, title = title)
                cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(top_pca, top_pca_by_gene, num_to_plot=plot_num, title=title, plot=False)
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
single_gene_sig('Tcf21', log2_expdf_cell)
#stability_ratio = clust_stability(log2_expdf_gene)
#print stability_ratio
cc_gene_df = cell_cycle(hu_cc_gene_df, log2_expdf_gene)


new_by_gene = log2_expdf_gene[gene_list]
new_by_cell = new_by_gene.transpose()
cell_linkage, plotted_df_by_gene, col_order = clust_heatmap(gene_list, new_by_gene, num_to_plot=len(gene_list))
#cell_dist, row_dist, row_clusters, link_mat, row_dendr = run_cluster(top_pca_by_gene)
cc = make_tree_json(cell_linkage, plotted_df_by_gene)
make_subclusters(cc, new_by_cell, gene_corr_list=['Tcf21', 'Pdgfra'])
find_twobytwo(cc, new_by_cell, log2_expdf_cell)
#augmented_dendrogram(row_clusters, labels=top_pca_by_cell.columns.tolist(), leaf_rotation=90, leaf_font_size=8)
