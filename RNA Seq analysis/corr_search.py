import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
import matplotlib.pyplot as plt
import matplotlib.pylab as pl


path_to_file ='/Volumes/Seq_data/results_pdgfra1_ctrl_pnx'

load_new_cells = False
load_new_sig = True
save_new_sig = False

if load_new_cells:
    fpbcell = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_by_cell.p'), 'rb')
    by_cell = pickle.load(fpbcell)
    fpcelllist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_cell_list.p'), 'rb')
    cell_list = pickle.load(fpcelllist)
    fpbgene = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_by_gene.p'), 'rb')
    by_gene = pickle.load(fpbgene)
    fpgenelist = open(os.path.join(path_to_file,'fpkm_cuff_pdgfra1_outlier_gene_list.p'), 'rb')
    gene_list = pickle.load(fpgenelist)


    df_by_gene = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
    df_by_cell = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)
    cols = df_by_cell.shape[0]
    rows = df_by_cell.shape[1]
    corr_by_gene = df_by_gene.corr()
    corr_by_cell = df_by_cell.corr()

    cor = corr_by_gene
    cor.loc[:,:] =  np.tril(cor.values, k=-1)
    cor = cor.stack()
    sig_corr_pos = cor[cor >=0.5]
    sig_corr_neg = cor[cor <=-0.5]

    with open(os.path.join(path_to_file,'gene_correlations_sig_neg.p'), 'wb') as fp:
        pickle.dump(sig_corr_neg, fp)
    with open(os.path.join(path_to_file,'gene_correlations_sig_pos.p'), 'wb') as fp0:
        pickle.dump(sig_corr_pos, fp0)
    with open(os.path.join(path_to_file,'by_gene_corr.p'), 'wb') as fp1:
        pickle.dump(corr_by_gene, fp1)
    with open(os.path.join(path_to_file,'by_cell_corr.p'), 'wb') as fp2:
        pickle.dump(corr_by_cell, fp2)

if load_new_sig:
    corr_by_gene_pos =  open(os.path.join(path_to_file,'gene_correlations_sig_pos.p'), 'rb')
    corr_by_gene_neg =  open(os.path.join(path_to_file,'gene_correlations_sig_neg.p'), 'rb')
    cor_pos = pickle.load(corr_by_gene_pos)
    cor_neg = pickle.load(corr_by_gene_neg)
    cor_pos_df = pd.DataFrame(cor_pos)
    cor_neg_df = pd.DataFrame(cor_neg)
    sig_corr = cor_pos_df.append(cor_neg_df)
    sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])
if save_new_sig:
    sig_corrs.to_csv(os.path.join(path_to_file,'pdgfra_corr_sig.txt'), sep = '\t')

term_to_search ='Dcn'
neg = True
for index, row in sig_corrs.iterrows():
    if term_to_search in index:
        neg = False
        if index[0]==term_to_search:
            print index[1], row['corr']
        else:
            print index[0], row['corr']
if neg:
    print term_to_search+' not correlated.'


corr_by_gene_pos.close()
corr_by_gene_neg.close()
