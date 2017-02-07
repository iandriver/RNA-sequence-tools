import datetime
import numpy as np
from sklearn.covariance import GraphLasso
from sklearn.decomposition import PCA as skPCA
import pandas as pd
import matplotlib

from sklearn.preprocessing import scale
import community
import matplotlib.pyplot as plt
import os
import networkx as nx
from sklearn import cluster, covariance, manifold
from matplotlib import collections


def return_top_pca_gene(by_cell_matrix, user_num_genes = None):
    gene_number = 100
    gene_pca = skPCA(n_components=3)
    np_by_gene = np.asarray(by_cell_matrix.transpose())
    gene_index = by_cell_matrix.index.tolist()

    if user_num_genes is not None:
        num_genes = user_num_genes
    else:
        num_genes = min(gene_number, len(gene_index))
    by_gene_trans = gene_pca.fit_transform(np_by_gene)
    Pc_df = pd.DataFrame(gene_pca.components_.T, columns=['PC-1', 'PC-2', 'PC-3'], index=gene_index)
    pca_rank_df = Pc_df.abs().sum(axis=1)
    Pc_sort_df = pca_rank_df.nlargest(len(gene_index))
    top_pca_list = Pc_sort_df.index.tolist()
    new_cell_matrix = by_cell_matrix.ix[top_pca_list[0:num_genes],:]
    return new_cell_matrix.transpose(), top_pca_list[0:num_genes]

def save_network_graph( matrix, labels, filename, title, scale=8, layout = "circular", weight = lambda x: abs(4*x)**(2.5) ):
    labels = dict( zip( range( len(labels) ), labels) )
    d = matrix.shape[0]
    D = nx.Graph(matrix)
    #D.add_nodes_from( range(d) )
    #for i in range(d):
    #	for j in range(i):
    #			if matrix[i,j] != 0:
    #				D.add_edge( i, j, weight = matrix[i,j])
    weights = [ D[x][y]['weight'] for x,y in D.edges() ]

    elarge=[(u,v) for (u,v,d) in D.edges(data=True) if d['weight'] >0.2]
    esmall=[(u,v) for (u,v,d) in D.edges(data=True) if d['weight'] <=0.2]
    #weights = weights/np.max( np.abs( weights ) )
    cmap = plt.get_cmap( "Reds" ) #or some other one

    fig = plt.figure(figsize=(50,50))
    ax = fig.add_subplot(111)
    if layout == "circular":
    	pos = nx.circular_layout( D, scale =scale )
    elif layout == "spring":
    	pos = nx.spring_layout( D ,scale = scale, iterations = 35 )
    #bweights = [ 1+100*(x-min(weights))/( max(weights)- min(weights) ) for x in weights ]
    bweights = [ 'k'*(z<0) + 'r'*(z>0) for z in weights ]
    width_small = [ weight(w) for w in weights if w < 0.2]
    width_large = [ weight(w) for w in weights if w >= 0.2]
    nx.draw_networkx_edges(D,pos,edgelist=elarge,
                    edge_color= "red", width=width_large, )
    nx.draw_networkx_edges(D,pos,edgelist=esmall,width=width_small,alpha=0.5,edgedge_color="blue",style='dashed')

    nx.draw_networkx_nodes( D, pos, ax=ax, node_size = 0, node_color="red")
    nx.draw_networkx_labels( D, pos,font_size=15, labels = labels, ax = ax)
    plt.axis("off")
    plt.title(title)
    plt.savefig( filename, bbox_inches="tight")

path_to_file = '/Users/iandriver/Downloads/monocle2_5state_groups_all_genes_nocc_scicast_analysis/count_matrix_after_filtering.txt'


#load file gene
by_cell = pd.DataFrame.from_csv(path_to_file, sep='\t')
by_gene = by_cell.transpose()
#create list of genes
gene_list = by_cell.index.tolist()
#create cell list
cell_list = [x for x in list(by_cell.columns.values)]


df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list, index=cell_list)
df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list)
log2_df_cell = np.log2(df_by_cell1+1)
alpha=0.4

top_pca_matrix, top_pca_genes = return_top_pca_gene(log2_df_cell,user_num_genes=100)
gl = covariance.GraphLassoCV()
gene_data = scale(top_pca_matrix.as_matrix())

gl.fit(gene_data)
_, labels = cluster.affinity_propagation(gl.covariance_)
n_labels = labels.max()
names = np.array(top_pca_genes)
prec_sp = gl.precision_


for i in range(n_labels + 1):
    print('Cluster %i: %s' % ((i + 1), ', '.join(names[labels == i])))



node_position_model = manifold.LocallyLinearEmbedding(
    n_components=2, eigen_solver='dense', n_neighbors=7)

embedding = node_position_model.fit_transform(gene_data.T).T


plt.figure(1, facecolor='w', figsize=(10, 8))
plt.clf()
ax = plt.axes([0., 0., 1., 1.])
plt.axis('off')

# Display a graph of the partial correlations
partial_correlations = gl.precision_.copy()
d = 1 / np.sqrt(np.diag(partial_correlations))
partial_correlations *= d
partial_correlations *= d[:, np.newaxis]
non_zero = (np.abs(np.triu(partial_correlations, k=1)) > 0.02)

# Plot the nodes using the coordinates of our embedding
cm = plt.cm.get_cmap('spectral')
#cmap = [cm(1.*i/len(embedding[0])) for i in range(len(embedding[0]))]


plt.scatter(embedding[0], embedding[1], s=100 * d ** 2, c=labels,
            cmap=cm)

# Plot the edges
start_idx, end_idx = np.where(non_zero)
#a sequence of (*line0*, *line1*, *line2*), where::
#            linen = (x0, y0), (x1, y1), ... (xm, ym)
segments = [[embedding[:, start], embedding[:, stop]]
            for start, stop in zip(start_idx, end_idx)]
values = np.abs(partial_correlations[non_zero])
lc = collections.LineCollection(segments,
                    zorder=0, cmap=plt.cm.hot_r,
                    norm=plt.Normalize(0, .7 * values.max()))
lc.set_array(values)
lc.set_linewidths(15 * values)
ax.add_collection(lc)

# Add a label to each node. The challenge here is that we want to
# position the labels to avoid overlap with other labels
for index, (name, label, (x, y)) in enumerate(
        zip(names, labels, embedding.T)):

    dx = x - embedding[0]
    dx[index] = 1
    dy = y - embedding[1]
    dy[index] = 1
    this_dx = dx[np.argmin(np.abs(dy))]
    this_dy = dy[np.argmin(np.abs(dx))]
    if this_dx > 0:
        horizontalalignment = 'left'
        x = x + .002
    else:
        horizontalalignment = 'right'
        x = x - .002
    if this_dy > 0:
        verticalalignment = 'bottom'
        y = y + .002
    else:
        verticalalignment = 'top'
        y = y - .002
    plt.text(x, y, name, size=10,
             horizontalalignment=horizontalalignment,
             verticalalignment=verticalalignment,
             bbox=dict(facecolor='w',
                       edgecolor=plt.cm.spectral(label / float(n_labels)),
                       alpha=.6))

plt.xlim(embedding[0].min() - .15 * embedding[0].ptp(),
         embedding[0].max() + .10 * embedding[0].ptp(),)
plt.ylim(embedding[1].min() - .03 * embedding[1].ptp(),
         embedding[1].max() + .03 * embedding[1].ptp())

plt.savefig(os.path.join(os.path.dirname(path_to_file),'linedraw_graph.pdf'),bbox_inches="tight")

def community_cluster(cov_sp, symbols):
	G = nx.Graph( cov_sp )
	partition = community.best_partition( G )
	for i in set(partition.values() ):
		print("Community: ",i)
		members = [ symbols[node] for node  in partition.keys() if partition[node] == i]
		print(members)


def affinity_cluster( cov_sp, symbols):
	print("Affinity cluster")
	ap = cluster.AffinityPropagation()
	ap.fit( cov_sp )
	labels = np.array(  ap.labels_ )
	for i in set(ap.labels_):
		print("Community: ",i)
		members = [ symbols[node] for node in np.nonzero( labels == i)[0]]
		print(members)
#community_cluster(gl.covariance_, top_pca_genes)
affinity_cluster(gl.covariance_, top_pca_genes)
save_network_graph( -prec_sp + np.diag( np.diagonal( prec_sp) ), names, os.path.join(os.path.dirname(path_to_file),"LargeNetworkNo_SP.pdf"), title="spring_sp_prec_test",layout="spring", scale= 10, weight = lambda x: abs(5*x)**(2.5) )
save_network_graph( gl.covariance_, names, os.path.join(os.path.dirname(path_to_file),"cov_diagram.pdf"), title="cov_test", scale = 8, layout = "spring" )
save_network_graph( gl.precision_, names, os.path.join(os.path.dirname(path_to_file),"precision.pdf") , title="Precision Matrix Network", layout="spring")
