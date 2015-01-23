import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
from sklearn.decomposition import PCA as skPCA
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.mlab import PCA
from mpl_toolkits.mplot3d import Axes3D
import scipy.cluster.hierarchy as sch



fpbcell = open('pdgfra_outlier_by_cell.p', 'rb')
by_cell = pickle.load(fpbcell)
fpcelllist = open('pdgfra_outlier_cell_list.p', 'rb')
cell_list = pickle.load(fpcelllist)
fpbgene = open('pdgfra_outlier_by_gene.p', 'rb')
by_gene = pickle.load(fpbgene)
fpgenelist = open('pdgfra_outlier_gene_list.p', 'rb')
gene_list = pickle.load(fpgenelist)

print cell_list, len(gene_list)


def fcluster( pts, ncluster, method="average", criterion="maxclust" ):
    """ -> (pts, Y pdist, Z linkage, T fcluster, clusterlists)
        ncluster = n1 + n2 + ... (including n1 singletons)
        av cluster size = len(pts) / ncluster
    """
    pts = np.asarray(pts)
    Y = scipy.spatial.distance.pdist( pts )  # ~ N^2 / 2
    Z = hier.linkage( Y, method )  # N-1
    T = hier.fcluster( Z, ncluster, criterion=criterion )
        # clusters = clusterlists(T)
    return (pts, Y, Z, T)

# Compute and plot dendrogram.
fig = pl.figure()
axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
Y = sch.linkage(by_gene, method='centroid')
Z = sch.dendrogram(Y, orientation='right')
axdendro.set_xticks([])
axdendro.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
index = Z['leaves']
by_cell = by_gene[index,:]
by_cell = by_gene[:,index]
im = axmatrix.matshow(by_gene, aspect='auto', origin='lower')
axmatrix.set_xticks([])
axmatrix.set_yticks([])

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
pl.colorbar(im, cax=axcolor)

# Display and save figure.
fig.show()
fig.savefig('dendrogram.png')

def twod_pca(by_gene, by_cell):
  pca = skPCA()
  print by_cell.shape, 'by_cell'
  print by_gene.shape, 'by gene'

  X_pca = pca.fit_transform(by_gene)
  pc1 = []
  pc2 = []
  for pcs in X_pca:
    pc1.append(pcs[0])
    pc2.append(pcs[1])
  np_pc1 = np.array(pc1)
  np_pc2 = np.array(pc2)
  plt.scatter(np_pc1, np_pc2)
  plt.show()


def threeD_pca_plot(to_pca):
  pca = PCA(to_pca)

  pca.Y

  x = []
  y = []
  z = []
  for item in pca.Y:
    x.append(item[0])
    y.append(item[1])
    z.append(item[2])

  pl.close('all') # close all latent plotting windows
  fig1 = pl.figure() # Make a plotting figure
  ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
  pltData = [x,y,z]
  ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data

  # make simple, bare axis lines through space:
  xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis
  ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
  yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
  ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
  zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
  ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.

  # label the axes
  ax.set_xlabel("x-axis label")
  ax.set_ylabel("y-axis label")
  ax.set_zlabel("y-axis label")
  ax.set_title("Cell PCA")
  pl.show() # show the plot

  fpbcell.close()
  fpcelllist.close()
  fpbgene.close()
  fpgenelist.close()
