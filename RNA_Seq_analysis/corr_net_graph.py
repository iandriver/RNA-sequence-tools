import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from operator import itemgetter
from matplotlib.ticker import LinearLocator
import networkx as nx
import plotly.plotly as py
from plotly.graph_objs import *
import plotly.tools as tls
import itertools

#plotly sign in information
py.sign_in('driver.ian', '0oql1n8y2r')
#base path to pickle files with fpkm or count matrix
path_to_file ='/Volumes/Seq_data/counts_sheppard_all/cpm_norm_kidney'
start_file_name = 'kidney_norm_cpm'

#if new cell pickle files need to be loaded
load_new_cells = True
#if you want to load a new signifcant correlation file
load_new_sig = True
#if you want save a new significant correlation file (pickle)
save_new_sig = False
#if you want to run a new correlation (can take a while)
run_corr = False

cmaps = {'Sequential':['Blues', 'BuGn', 'BuPu',
'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'],
'Sequential (2)': ['afmhot', 'autumn', 'bone', 'cool', 'copper',
'gist_heat', 'gray', 'hot', 'pink',
'spring', 'summer', 'winter'],
'Diverging':['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
'seismic'],
'Qualitative':['Accent', 'Dark2', 'Paired', 'Pastel1',
'Pastel2', 'Set1', 'Set2', 'Set3'],
'Miscellaneous':['gist_earth', 'terrain', 'ocean', 'gist_stern',
'brg', 'CMRmap', 'cubehelix',
'gnuplot', 'gnuplot2', 'gist_ncar',
'nipy_spectral', 'jet', 'rainbow',
'gist_rainbow', 'hsv', 'flag', 'prism']}

if load_new_cells:
    fpbcell = open(os.path.join(path_to_file, start_file_name+'_outlier_by_cell.p'), 'rb')
    by_cell = pickle.load(fpbcell)
    fpbcell.close()
    fpcelllist = open(os.path.join(path_to_file, start_file_name+'_outlier_cell_list.p'), 'rb')
    cell_list = pickle.load(fpcelllist)
    fpcelllist.close()
    fpbgene = open(os.path.join(path_to_file, start_file_name+'_outlier_by_gene.p'), 'rb')
    by_gene = pickle.load(fpbgene)
    fpbgene.close()
    fpgenelist = open(os.path.join(path_to_file, start_file_name+'_outlier_gene_list.p'), 'rb')
    gene_list = pickle.load(fpgenelist)
    fpgenelist.close()


    df_by_gene = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
    df_by_cell = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)
    cols = df_by_cell.shape[0]
    rows = df_by_cell.shape[1]

#run correlation matrix and save only those above threshold
if run_corr:
    corr_by_gene = df_by_gene.corr(method='spearman', min_periods=3)
    corr_by_cell = df_by_cell.corr()

    cor = corr_by_gene
    cor.loc[:,:] =  np.tril(cor.values, k=-1)
    cor = cor.stack()
    sig_corr_pos = cor[cor >=0.5]
    sig_corr_neg = cor[cor <=-0.5]

    with open(os.path.join(path_to_file,'gene_correlations_sig_neg_spearman.p'), 'wb') as fp:
        pickle.dump(sig_corr_neg, fp)
    with open(os.path.join(path_to_file,'gene_correlations_sig_pos_spearman.p'), 'wb') as fp0:
        pickle.dump(sig_corr_pos, fp0)
    with open(os.path.join(path_to_file,'by_gene_corr.p'), 'wb') as fp1:
        pickle.dump(corr_by_gene, fp1)
    with open(os.path.join(path_to_file,'by_cell_corr.p'), 'wb') as fp2:
        pickle.dump(corr_by_cell, fp2)

#load or save significant correlations as needed
if load_new_sig:
    corr_by_gene_pos =  open(os.path.join(path_to_file,'gene_correlations_sig_pos_spearman.p'), 'rb')
    corr_by_gene_neg =  open(os.path.join(path_to_file,'gene_correlations_sig_neg_spearman.p'), 'rb')
    cor_pos = pickle.load(corr_by_gene_pos)
    cor_neg = pickle.load(corr_by_gene_neg)
    cor_pos_df = pd.DataFrame(cor_pos)
    cor_neg_df = pd.DataFrame(cor_neg)
    sig_corr = cor_pos_df.append(cor_neg_df)
    sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])
if save_new_sig:
    sig_corrs.to_csv(os.path.join(path_to_file,'kidney_counts_corr_sig_spearman.txt'), sep = '\t')

def find_gen_rank(g, split_on='_', pos=1, cat_name=['d4pnx', 'ctrl']):

    sorted_df = by_cell.sort([g])
    score_on = 'd4pnx'
    g_df = sorted_df[g]
    ranked_cells = sorted_df.index.values
    ranked_cat = [x.split(split_on)[pos] for x in ranked_cells]
    div_by = int(len(ranked_cat)/len(cat_name))
    start = div_by *(len(cat_name)-1)
    score1 = len([x for x in ranked_cat[start:len(ranked_cat)] if x == score_on])
    tot = len([x for x in ranked_cat if x == score_on])
    res_score = float(score1)/float(tot)
    return "%.2f" % res_score

#gene to search
term_to_search ='Sparc'

#corr_plot finds and plots all correlated genes, log turns on log scale, sort plots the genes in the rank order of the gene searched
def corr_plot(term_to_search, log=False, sort=False):
    corr_tup = [(term_to_search, 1)]
    neg = True
    rank = False
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
    corr_df.to_csv(os.path.join(path_to_file, 'Corr_w_'+term_to_search+'_list.txt'), sep = '\t', index=False)
    for c in corr_tup:
        print c
    to_plot = [x[0] for x in corr_tup]
    sorted_df = df_by_gene.sort([term_to_search])
    log2_df = np.log2(df_by_gene[to_plot])
    sorted_log2_df=np.log2(sorted_df[to_plot])
    ylabel='Counts'
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
    if rank:
        ax.set_title('Correlates with '+term_to_search+'.  Percent seperate PNX: '+find_gen_rank(term_to_search))
    else:
        ax.set_title('Correlates with '+term_to_search)
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
    plt.show()

#network plot plots correlated genes as network
def network_plot(term_to_search):
    MG_genes=nx.DiGraph()
    MG_genes.add_node(term_to_search, size=10, color='red')
    network_tup = []
    neg = True
    for index, row in sig_corrs.iterrows():
        if term_to_search in index:
            neg = False
            if index[0]==term_to_search:
                network_tup.append((term_to_search, index[1], row['corr']))
            else:
                network_tup.append((term_to_search, index[0], row['corr']))
    if neg:
        print term_to_search+' not correlated.'
    else:
        network_tup.sort(key=itemgetter(2), reverse=True)
        labels= {term_to_search:term_to_search}
        for n1, n2, w in network_tup:
            MG_genes.add_node(n2, size=10, label=n2, color='blue')
            labels[n2]=n2
        MG_genes.add_weighted_edges_from(network_tup)
        pos = nx.spring_layout(MG_genes)
        nx.draw(MG_genes, pos, with_labels=True, node_color='#A0CBE2')
        plt.show()

corr_plot(term_to_search, sort=True)

corr_by_gene_pos.close()
corr_by_gene_neg.close()
