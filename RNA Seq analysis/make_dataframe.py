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

shared_df = pd.DataFrame.from_csv('/Volumes/Seq_data/Pdgfra2_all_fpkm_analysis/fpkm_cuff_pdgfra2_outlier_filtered.txt', sep='\t')
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

def make_new_matrix(org_matrix_by_cell, gene_list):
    new_gmatrix_df = org_matrix_by_cell[gene_list]
    new_cmatrix_df = new_gmatrix_df.transpose()
    new_cmatrix_df.to_csv(os.path.join(path_to_file, 'goterms_top5000_count_matrix.txt'), sep = '\t', index=gene_list)


make_new_matrix(by_cell, singular_gene_list())
