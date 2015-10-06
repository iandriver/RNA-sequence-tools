import sys
import scipy as SP
import pylab as PL
import numpy as np
from matplotlib import cm
import h5py
import os
import GPy
import pandas as pd
import cPickle as pickle

#adjust path
scLVM_BASE = '/Volumes/Seq_data/count-picard_combined_ips17_BU3'
from scLVM import scLVM

#sys.path.append(scLVM_BASE)
#sys.path.append( scLVM_BASE +'..')
#sys.path.append(scLVM_BASE + 'scLVM/utils') #this is not included in the github repo
#sys.path.append(scLVM_BASE +'CFG')
#from misc import *
#from barplot import *
#from default import *
from scLVM.utils.barplot import *
from scLVM.utils.misc import *

from IPython.display import Latex

data = os.path.join(scLVM_BASE,'data_ips__normCounts.h5f')
f = h5py.File(data,'r')
# gene expression matrix
Y = f['LogNcountsGene'][:]
# technical noise
tech_noise = f['LogVar_techGene'][:]
# index of heterogeneous genes
genes_het_bool=f['genes_heterogen'][:]
# gene names
geneID = f['gene_names'][:]
cell_names = f['cell_names'][:]
# idx of cell cycle genes from GO
cellcyclegenes_filter = SP.unique(f['cellcyclegenes_filter'][:].ravel() -1)
# idx of cell cycle genes from cycle base
cellcyclegenes_filterCB = f['ccCBall_gene_indices'][:].ravel() -1

# filter cell cycle genes
idx_cell_cycle = SP.union1d(cellcyclegenes_filter,cellcyclegenes_filterCB)
# determine non-zero counts
idx_nonzero = SP.nonzero((Y.mean(0)**2)>0)[0]
idx_cell_cycle_noise_filtered = SP.intersect1d(idx_cell_cycle,idx_nonzero)
# subset gene expression matrix
Ycc = Y[:,idx_cell_cycle_noise_filtered]


plt = PL.subplot(1,1,1)
PL.imshow(Ycc,cmap=cm.RdBu,vmin=-3,vmax=+3,interpolation='None')
#PL.colorbar()
plt.set_xticks([])
plt.set_yticks([])
PL.xlabel('genes')
PL.ylabel('cells')
PL.ylabel('Variance explained')
PL.show()

k = 80                    # number of latent factors
out_dir = scLVM_BASE + 'cache'       # folder where results are cached
file_name = 'Kcc.hdf5'    # name of the cache file
recalc = True             # recalculate X and Kconf
use_ard = True            # use automatic relevance detection
sclvm = scLVM(Y)
#Fit model with 80 factors
X_ARD,Kcc_ARD,varGPLVM_ARD = sclvm.fitGPLVM(idx=idx_cell_cycle_noise_filtered,k=k,out_dir=out_dir,file_name=file_name,recalc=recalc, use_ard=use_ard)

#Plot variance contributions from ARD
plt = PL.subplot(1,1,1)
PL.title('Variance explained by latent factors')
PL.scatter(SP.arange(k)+1,varGPLVM_ARD['X_ARD'])
PL.xlim([0,k+1])
PL.xlabel('# Factor')
PL.ylabel('Variance explained')
PL.show()

#Fit model with a single factor (rank 1 covariance matrix)
X,Kcc,varGPLVM = sclvm.fitGPLVM(idx=idx_cell_cycle_noise_filtered,k=1,out_dir='./cache',file_name=file_name,recalc=True, use_ard=False)

#Plot inferred similarity matrix
plt = PL.subplot(1,1,1)
PL.title('Similarity matrix based on cell cycle')
PL.imshow(Kcc,cmap=cm.RdBu,vmin=-3,vmax=+3,interpolation='None')
PL.colorbar()
plt.set_xticks([])
plt.set_yticks([])
PL.xlabel('cells')
PL.ylabel('cells')
PL.ylabel('Variance explained')
PL.show()

# considers only heterogeneous genes
Ihet = genes_het_bool==1
Y    = Y[:,Ihet]
tech_noise = tech_noise[Ihet]
geneID = geneID[Ihet]
print geneID
print len(geneID)
#optionally: restrict range for the analysis
i0 = 0    # gene from which the analysis starts
i1 = 6025 # gene at which the analysis ends

# construct sclvm object
sclvm = scLVM(Y,geneID=geneID,tech_noise=tech_noise)

# fit the model from i0 to i1
sclvm.varianceDecomposition(K=Kcc,i0=i0,i1=i1)

normalize=True    # variance components are normalizaed to sum up to one

# get variance components
var, var_info = sclvm.getVarianceComponents(normalize=normalize)
var_filtered = var[var_info['conv']] # filter out genes for which vd has not converged

# get corrected expression levels
Ycorr = sclvm.getCorrectedExpression()
with open(os.path.join(scLVM_BASE,'scLVM_Ycorr.p'), 'wb') as fp1:
    pickle.dump(Ycorr, fp1)
print Ycorr.shape
Ycorr_df = pd.DataFrame(Ycorr, index=cell_names, columns=geneID)
Ycorr_df.to_csv(os.path.join(scLVM_BASE,'scLVM_Ycorr.txt'), sep = '\t')
#calculate average variance components across all genes and visualize
var_mean = var_filtered.mean(0)
colors = ['Green','MediumBlue','Gray']
pp=PL.pie(var_mean,labels=var_info['col_header'],autopct='%1.1f%%',colors=colors,shadow=True, startangle=0)
PL.show()

H2=1-var_filtered[:,2]
var_comp_fileds = SP.array([[0, 'cell cycle', 'Peru'],
                            [1, 'biol. var', 'DarkMagenta'],
                            [2, 'tech. var', '#92c5de']], dtype=object)
var_plot(var_filtered,H2,var_comp_fileds,normalize=True, figsize=[5,4])
PL.show()

i0 = 0     # gene from which the analysis starts
i1 = 20    # gene to which the analysis ends

# fit lmm without correction
pv0,beta0,info0 = sclvm.fitLMM(K=None,i0=i0,i1=i1,verbose=False)
# fit lmm with correction
pv1,beta1,info1 = sclvm.fitLMM(K=Kcc,i0=i0,i1=i1,verbose=False)

plt=PL.subplot(2,2,1)
PL.title('Without Correction')
p=PL.imshow(beta0[:,i0:i1],cmap=cm.RdBu,vmin=-0.6,vmax=+1,interpolation='None')
PL.colorbar()
plt.set_xticks([])
plt.set_yticks([])
PL.xlabel('gene'),PL.ylabel('gene')
plt=PL.subplot(2,2,2)
PL.title('With Correction')
p=PL.imshow(beta1[:,i0:i1],cmap=cm.RdBu,vmin=-0.6,vmax=+1,interpolation='None')
PL.colorbar()
plt.set_xticks([])
plt.set_yticks([])
PL.xlabel('gene'),PL.ylabel('gene')
PL.show()
np.savetxt(os.path.join(scLVM_BASE,'Ycorr.txt'),Ycorr)

# Model optimization
Ystd = Ycorr-Ycorr.mean(0)
Ystd/=Ystd.std(0)
input_dim = 2 # How many latent dimensions to use
kern = GPy.kern.RBF(input_dim,ARD=True) # ARD kernel
m = GPy.models.BayesianGPLVM(Ystd, input_dim=input_dim, kernel=kern, num_inducing=40)
m.optimize('scg', messages=0, max_iters=2000)
m.kern.plot_ARD()
PL.show()
gene_to_plot = 'NKX2-1'
i_nkx = SP.where(geneID==gene_to_plot)
color = Ycorr[:,i_nkx]
#color = Ycorr[:,0]
PL.scatter(m.X[:,0]['mean'], m.X[:,1]['mean'], 40, color)
PL.xlabel('PC1')
PL.ylabel('PC2')
PL.title(gene_to_plot)
PL.colorbar()
PL.show()

[S,W] = PCA(Ystd,2)
PL.scatter(S[:,0],S[:,1], 40, color)
PL.xlabel('PC1')
PL.ylabel('PC2')
PL.colorbar()
PL.title(gene_to_plot)
PL.show()
