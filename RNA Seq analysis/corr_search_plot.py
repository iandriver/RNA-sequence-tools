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

py.sign_in('driver.ian', '0oql1n8y2r')
path_to_file ='/Volumes/Seq_data/results_macs_pnx'

load_new_cells = True
load_new_sig = True
save_new_sig = False
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
    fpbcell = open(os.path.join(path_to_file,'fpkm_cuff_macs1_outlier_by_cell.p'), 'rb')
    by_cell = pickle.load(fpbcell)
    fpcelllist = open(os.path.join(path_to_file,'fpkm_cuff_macs1_outlier_cell_list.p'), 'rb')
    cell_list = pickle.load(fpcelllist)
    fpbgene = open(os.path.join(path_to_file,'fpkm_cuff_macs1_outlier_by_gene.p'), 'rb')
    by_gene = pickle.load(fpbgene)
    fpgenelist = open(os.path.join(path_to_file,'fpkm_cuff_macs1_outlier_gene_list.p'), 'rb')
    gene_list = pickle.load(fpgenelist)


    df_by_gene = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
    df_by_cell = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)
    cols = df_by_cell.shape[0]
    rows = df_by_cell.shape[1]

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
    sig_corrs.to_csv(os.path.join(path_to_file,'macs_corr_sig_spearman.txt'), sep = '\t')

term_to_search ='Vegfb'
def corr_plot(term_to_search, log=False, sort=False):
    corr_tup = [(term_to_search, 1)]
    neg = True
    fig, ax = plt.subplots()
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
    for c in corr_tup:
        print c
    to_plot = [x[0] for x in corr_tup]
    sorted_df = df_by_gene.sort([term_to_search])
    log2_df = np.log2(df_by_gene[to_plot])
    sorted_log2_df=np.log2(sorted_df[to_plot])
    ylabel='FPKM'
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
    ax.xaxis.set_minor_locator(LinearLocator(numticks=len(xlabels)))
    ax.set_xticklabels(xlabels, minor=True, rotation='vertical')
    ax.tick_params(axis='x', labelsize=5)
    ax.legend(loc=1,prop={'size':6})
    plt.show()

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
