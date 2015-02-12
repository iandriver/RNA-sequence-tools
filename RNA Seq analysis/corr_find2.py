import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
import matplotlib.pyplot as plt
import matplotlib.pylab as pl



path_to_file = '/Volumes/Seq_data/results_pdgfra1_ctrl_pnx'
corr_by_gene_pos =  open(os.path.join(path_to_file,'gene_correlations_sig_pos.p'), 'rb')
corr_by_gene_neg =  open(os.path.join(path_to_file,'gene_correlations_sig_neg.p'), 'rb')
cor_pos = pickle.load(corr_by_gene_pos)
cor_neg = pickle.load(corr_by_gene_neg)
cor_pos_df = pd.DataFrame(cor_pos)
cor_neg_df = pd.DataFrame(cor_neg)
cor_pos_df.to_csv(os.path.join(path_to_file,'pdgfra_cor_pos_sig.txt'), sep = '\t')
cor_neg_df.to_csv(os.path.join(path_to_file,'pdgfra_cor_neg_sig.txt'), sep = '\t')

sig_corr = cor_pos_df.append(cor_neg_df)
sig_corrs = pd.DataFrame(sig_corr[0], columns=["corr"])
term_to_search ='Dcn'
for index, row in sig_corrs.iterrows():
    if term_to_search in index:
        if index[0]=='Dcn':
            print index[1], row['corr']
        else:
            print index[0], row['corr']


corr_by_gene_pos.close()
corr_by_gene_neg.close()
