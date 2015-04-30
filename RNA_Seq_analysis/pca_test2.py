import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
import plotly.plotly as py
from plotly.graph_objs import *
import plotly.tools as tls
from sklearn.decomposition import PCA as skPCA
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.mlab import PCA
from mpl_toolkits.mplot3d import Axes3D
import scipy.cluster.hierarchy as sch
from sklearn.preprocessing import StandardScaler

path_to_file ='/Volumes/Seq_data/results_spc1_ctrl_pnx'
def make_pickle(file_to_dump, p_name='gene_annotation.p'):
    path_to_dump= os.path.join(path_to_file, p_name)
    with open(path_to_dump, 'wb') as fp:
        pickle.dump(file_to_dump, fp)

def load_pickle(p_name='gene_annotation.p'):
    path_to_load= os.path.join(path_to_file, p_name)
    with open(path_to_load, 'rb') as fp:
        data = pickle.load(fp)
    return data



fpbcell = open(os.path.join(path_to_file,'fpkm_cuff_spc1_outlier_by_cell.p'), 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open(os.path.join(path_to_file,'fpkm_cuff_spc1_outlier_cell_list.p'), 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open(os.path.join(path_to_file,'fpkm_cuff_spc1_outlier_by_gene.p'), 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open(os.path.join(path_to_file,'fpkm_cuff_spc1_outlier_gene_list.p'), 'rb')
gene_list = pickle.load(fpgenelist)


df_by_gene = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
df_by_cell = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)
cols = df_by_cell.shape[0]
rows = df_by_cell.shape[1]

X= df_by_gene.ix[:,0:].values
y = gene_list

X_std = StandardScaler().fit_transform(X)
cov_mat = np.cov(X_std.T)

run_new = False
eig_run = False
load_evs = False
eig_pair_load =False
if run_new:
    eig_vals, eig_vecs = np.linalg.eig(cov_mat)
    make_pickle(eig_vals, p_name="eig_vals.p")
    make_pickle(eig_vecs, p_name="eig_vecs.p")
elif load_evs:
    eig_vecs = load_pickle(p_name="eig_vecs.p")
else:
    eig_vals = load_pickle(p_name="eig_vals.p")

if eig_run:
    # Make a list of (eigenvalue, eigenvector) tuples
    eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]

    # Sort the (eigenvalue, eigenvector) tuples from high to low
    eig_pairs.sort(key=lambda tup: tup[0])
    eig_pairs.reverse()
    make_pickle(eig_pairs, p_name="eig_pairs.p")
elif eig_pair_load:
    eig_pairs = load_pickle(p_name="eig_pairs.p")

    # Visually confirm that the list is correctly sorted by decreasing eigenvalues
    print('Eigenvalues in descending order:')
    for i in eig_pairs[0:10]:
        print i[0]
eig_val = np.abs(eig_vals)
print eig_val[0:10]
tot = sum(eig_val)
var_exp = [(i / tot)*100 for i in sorted(eig_val, reverse=True)]
print var_exp
cum_var_exp = np.cumsum(var_exp)

trace1 = Bar(x=['PC %s' %i for i in range(1,20)], y=var_exp, showlegend=False)

trace2 = Scatter(x=['PC %s' %i for i in range(1,20)],y=cum_var_exp,name='cumulative explained variance')

data = Data([trace1, trace2])

layout=Layout(yaxis=YAxis(title='Explained variance in percent'),title='Explained variance by different principal components')

fig = Figure(data=data, layout=layout)
py.iplot(fig)
