import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from operator import itemgetter
from matplotlib.ticker import LinearLocator
import itertools

#***fill in these variables***



#base path to pickle files with fpkm or count matrix
path_to_file = '/Volumes/Seq_data/results_pdgfra_all_n2'
#for labeling all output files
start_file_name = 'fpkm_pdgfra_n2'

#gene to search
term_to_search =raw_input('Enter gene name to search correlation:')

#if you need to run a new correlation (can take a while)
run_corr = False

#difine the threshold for significant correlation (0-1 one being perfect correlation)
sig_threshold =0.5
#define correlation method options are: 'pearson', 'kendall', 'spearman'
method_name = 'spearman'
#Minimum number of observations required per pair of columns to have a valid result. Currently only available for pearson and spearman correlation.
min_period = 3

#if you want the correlation plot to be sorted by expression of the searched gene
plot_sort = True
#if you want the plot to be plotted on log2 scale
plot_log = False
#if you want save a new significant correlation file (pickle)
save_new_sig = True

#can rank genes by category separation (single gene clustering) if true define categories
#in find_gen_rank function
rank = False
#***only edit find_gen_rank categories***

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



#name of file containing gene
gene_file_source = 'go_search_genes_lung_all.txt'
df_by_gene1 = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
df_by_cell1 = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)

def make_new_matrix(org_matrix_by_cell, gene_list_file):
    split_on='_'
    gene_df = pd.read_csv(os.path.join(path_to_file, gene_list_file), delimiter= '\t')
    gene_list = gene_df['GeneID'].tolist()
    group_list = gene_df['GroupID'].tolist()
    gmatrix_df = org_matrix_by_cell[gene_list]
    cmatrix_df = gmatrix_df.transpose()
    cell_list1 = []
    for cell in cmatrix_df.columns.values:
        if cell.split(split_on)[1] == 'ctrl' or cell.split(split_on)[1] == 'pnx':
            if cell.split(split_on)[2][0] =='C':
                print cell, 'cell'
                cell_list1.append(cell)
    new_cmatrix_df = cmatrix_df[cell_list1]
    new_gmatrix_df = new_cmatrix_df.transpose()
    return new_cmatrix_df, new_gmatrix_df

df_by_cell, df_by_gene = make_new_matrix(df_by_gene1, gene_file_source)

#run correlation matrix and save only those above threshold
if run_corr:
    if method_name != 'kendall':
        corr_by_gene = df_by_gene.corr(method=method_name, min_periods=min_period)
    else:
        corr_by_gene = df_by_gene.corr(method=method_name)
    corr_by_cell = df_by_cell.corr()

    cor = corr_by_gene
    cor.loc[:,:] =  np.tril(cor.values, k=-1)
    cor = cor.stack()
    sig_corr_pos = cor[cor >=sig_threshold]
    sig_corr_neg = cor[cor <=(sig_threshold*-1)]

    with open(os.path.join(path_to_file,'gene_correlations_sig_neg_'+method_name+'.p'), 'wb') as fp:
        pickle.dump(sig_corr_neg, fp)
    with open(os.path.join(path_to_file,'gene_correlations_sig_pos_'+method_name+'.p'), 'wb') as fp0:
        pickle.dump(sig_corr_pos, fp0)
    with open(os.path.join(path_to_file,'by_gene_corr.p'), 'wb') as fp1:
        pickle.dump(corr_by_gene, fp1)
    with open(os.path.join(path_to_file,'by_cell_corr.p'), 'wb') as fp2:
        pickle.dump(corr_by_cell, fp2)


corr_by_gene_pos =  open(os.path.join(path_to_file,'gene_correlations_sig_pos_'+method_name+'.p'), 'rb')
corr_by_gene_neg =  open(os.path.join(path_to_file,'gene_correlations_sig_neg_'+method_name+'.p'), 'rb')
cor_pos = pickle.load(corr_by_gene_pos)
cor_neg = pickle.load(corr_by_gene_neg)
cor_pos_df = pd.DataFrame(cor_pos)
cor_neg_df = pd.DataFrame(cor_neg)
sig_corr = cor_pos_df.append(cor_neg_df)
sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])

if save_new_sig:
    sig_corrs.to_csv(os.path.join(path_to_file, start_file_name.split('_')[0]+'_counts_corr_sig_'+method_name+'.txt'), sep = '\t')

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


#corr_plot finds and plots all correlated genes, log turns on log scale, sort plots the genes in the rank order of the gene searched
def corr_plot(term_to_search, log=plot_log, sort=plot_sort):
    plt.clf()
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
    corr_df.to_csv(os.path.join(path_to_file, 'Corr_w_'+term_to_search+'_list.txt'), sep = '\t', index=False)
    for c in corr_tup:
        print c
    to_plot = [x[0] for x in corr_tup]
    sorted_df = df_by_gene.sort([term_to_search])
    print to_plot
    log2_df = np.log2(df_by_cell[to_plot])
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
    plt.savefig(os.path.join(path_to_file, start_file_name+'_corr_with_'+term_to_search), bbox_inches='tight')

    plt.show()



corr_plot(term_to_search)

corr_by_gene_pos.close()
corr_by_gene_neg.close()
