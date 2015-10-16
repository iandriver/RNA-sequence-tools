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


#run correlation matrix and save only those above threshold
def run_corr(df_by_gene, method_name='pearson', run_new=True):
    if run_new:
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
    else:
        corr_by_gene_pos =  open(os.path.join(path_to_file,'gene_correlations_sig_pos_'+method_name+'.p'), 'rb')
        corr_by_gene_neg =  open(os.path.join(path_to_file,'gene_correlations_sig_neg_'+method_name+'.p'), 'rb')
        cor_pos = pickle.load(corr_by_gene_pos)
        cor_neg = pickle.load(corr_by_gene_neg)

    cor_pos_df = pd.DataFrame(cor_pos)
    cor_neg_df = pd.DataFrame(cor_neg)
    sig_corr = cor_pos_df.append(cor_neg_df)
    sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])

    if run_new:
        sig_corrs.to_csv(os.path.join(path_to_file, base_name+'_counts_corr_sig_'+method_name+'.txt'), sep = '\t')
    return sig_corrs

#corr_plot finds and plots all correlated genes, log turns on log scale, sort plots the genes in the rank order of the gene searched
def corr_plot(terms_to_search, df_by_gene, log=False, sort=True):
    sig_corrs = run_corr(df_by_gene)
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
        plt.savefig(os.path.join(path_to_file, base_name+'_corr_with_'+term_to_search), bbox_inches='tight')
        plt.show()


if __name__ == '__main__':
    corr_plot()
