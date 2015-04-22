import os
import cPickle as pickle
import pandas as pd
import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.ticker import LinearLocator
import seaborn as sns
import numpy as np

shared_df = pd.DataFrame.from_csv('/Volumes/Seq_data/counts_sheppard_all/shared_genes_all.txt', sep='\t')
shared_list = []
for x in shared_df.index:
    shared_list.append(x.rstrip())

path_to_file ='/Volumes/Seq_data/counts_sheppard_all'
start_file_name = 'norm_cpm'
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

def sort_df_bylist(num_to_show=50):
    find_exp_list =[]
    for g in shared_list:
        if g in list(by_gene.columns.values):
            find_exp_list.append(g)

    shared_df = by_gene[find_exp_list]
    shared_g_df = shared_df.transpose()
    sorted_df = shared_df.sort(find_exp_list,ascending=False)
    plot_df = shared_g_df[sorted_df.index.values[0:num_to_show]]
    print sorted_df.index.values[0:num_to_show]

gene_list_toplot = ['Pdgfra', 'Pdgfrb', 'Col1a1', 'Acta2', 'Col3a1', 'Ctss', 'Cd48', 'Ptprc', 'Tgfbr1', 'Tgfbr2', 'Tgfbr3']

term_to_sort = 'Col1a1'
sorted_df = by_cell.sort([term_to_sort])

selected_df = sorted_df[gene_list_toplot]

def plot_genes(df_to_plot, log=False):
    fig, ax = plt.subplots()
    if log:
        log2_df = np.log2(df_to_plot)
        ax = log2_df.plot()
    else:
        ax = df_to_plot.plot()
    xlabels = df_to_plot.index.values
    ax.xaxis.set_minor_locator(LinearLocator(numticks=len(xlabels)))
    ax.set_xticklabels(xlabels, minor=True, rotation='vertical')
    ax.tick_params(axis='x', labelsize=5)
    ax.legend(loc=1,prop={'size':6})
    plt.show()

    plt.clf()
    cg = sns.clustermap(df_to_plot.transpose(), z_score=0)
    plt.show()

plot_genes(selected_df, log=True)
