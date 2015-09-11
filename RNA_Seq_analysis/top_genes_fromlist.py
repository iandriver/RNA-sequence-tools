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
from operator import itemgetter

shared_df = pd.DataFrame.from_csv('/Volumes/Seq_data/cuffnorm_spc_d0_4_7', sep='\t')
shared_list = []
for x in shared_df.index:
    shared_list.append(x.rstrip())

path_to_file ='/Volumes/Seq_data/Spc2_counts'
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

def singular_gene_list(path=path_to_file, file = '/go_search_genes.txt'):
    gene_df = pd.read_csv(path+file, delimiter= '\t')
    gene_list = gene_df['GeneID'].tolist()
    final_list = []
    for g in gene_list:
        final_list.append(g.split('_')[0])
    return final_list

def sort_df_bylist(shared_list, num_to_show=50):
    find_exp_list =[]
    for g in shared_list:
        if g in list(by_gene.columns.values):
            find_exp_list.append(g)

    shared_df = by_gene[find_exp_list]
    shared_g_df = shared_df.transpose()
    sorted_df = shared_df.sort(find_exp_list,ascending=False)
    plot_df = shared_g_df[sorted_df.index.values[0:num_to_show]]
    print sorted_df.index.values[0:num_to_show]

#a list of genes to search for ranking categories (how well does expression level define categories)
#split_on='sybol to split cell name', pos=index number of category name after split, cat_name=[list of category names]
def find_gen_rank_one(gen_list, split_on='_', pos=1, cat_name=['pnx', 'ctrl']):

    score_tup =[]
    for g in gen_list:
        sorted_df = by_cell.sort([g])
        g_df = sorted_df[g]
        ranked_cells = sorted_df.index.values
        ranked_cat = [x.split(split_on)[pos] for x in ranked_cells]
        div_by = int(len(ranked_cat)/len(cat_name))
        start = 0
        last_div = div_by
        score_list = np.zeros(len(cat_name)*len(cat_name))
        store_ind = 0
        for i,n in enumerate(cat_name):
            if last_div > len(ranked_cat):
                last_div = len(ranked_cat)
            for z in range(len(cat_name)):
                score1 = len([x for x in ranked_cat[start:last_div] if x == cat_name[z]])
                tot = len([x for x in ranked_cat if x == cat_name[z]])
                res_score = float(score1)/float(tot)
                score_list[store_ind] = res_score
                store_ind+=1

            start= last_div
            last_div +=div_by
        if max(score_list) >= 0.7 and np.count_nonzero(g_df.values) >= 15:
            score_tup.append((g,max(score_list)))
    score_tup.sort(key=itemgetter(1), reverse=True)
    score_df = pd.DataFrame(score_tup, columns=['GeneID', 'Maxscore'])
    score_df.to_csv(os.path.join(path_to_file, 'sep_score_genes_all.txt'), sep = '\t', index=False)

def find_gen_rank(g_list, split_on='_', pos=1, cat_name=['pnx', 'ctrl']):
    score_tup =[]
    for g in g_list:
        sorted_df = by_cell.sort([g])
        score_on = 'pnx'
        g_df = sorted_df[g]
        ranked_cells = sorted_df.index.values
        ranked_cat = [x.split(split_on)[pos] for x in ranked_cells]
        div_by = int(len(ranked_cat)/len(cat_name))
        start = div_by *(len(cat_name)-1)
        score1 = len([x for x in ranked_cat[start:len(ranked_cat)] if x == score_on])
        tot = len([x for x in ranked_cat if x == score_on])
        res_score = float(score1)/float(tot)
        score = "%.3f" % res_score
        score_tup.append((g,float(score)))
    score_tup.sort(key=itemgetter(1), reverse=True)
    score_df = pd.DataFrame(score_tup, columns=['GeneID', 'Maxscore'])
    score_df.to_csv(os.path.join(path_to_file, 'sep_score_genes_all.txt'), sep = '\t', index=False)

gene_list_toplot = singular_gene_list()[0:100]

term_to_sort = 'Hbegf'
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
    ax.tick_params(axis='x', labelsize=1)
    ax.legend(loc=1,prop={'size':6})
    plt.show()

    plt.clf()
    cg = sns.clustermap(df_to_plot.transpose(), z_score=1)
    plt.show()

plot_genes(selected_df, log=False)
#find_gen_rank(gene_list)
