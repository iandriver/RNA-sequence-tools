import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
import matplotlib.pyplot as plt
import matplotlib.pylab as pl



path_to_file = '/Volumes/Seq_data/results_pdgfra1_ctrl_pnx'
corr_by_gene =  open(os.path.join(path_to_file,'by_gene_corr.p'), 'rb')
corr_gene = pickle.load(corr_by_gene)
corr_by_cell = open(os.path.join(path_to_file,'by_cell_corr.p'), 'rb')
corr_cell = pickle.load(corr_by_cell)

cor = corr_by_gene
cor.loc[:,:] =  np.tril(cor.values, k=-1)
cor = cor.stack()
sig_corr = cor[cor >=0.5]

print sig_corr
with open(os.path.join(path_to_file,'gene_correlations_sig.p'), 'wb') as fp1:
    pickle.dump(sig_corr, fp1)


corr_by_gene.close()
corr_by_cell.close()
