import cPickle as pickle
import numpy as np
import pandas as pd
import scipy
import os
import matplotlib.pyplot as plt
import networkx as nx

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


def draw_graph(graph, labels=None, graph_layout='shell',
    node_size=500, node_color='blue', node_alpha=0.3,
    node_text_size=12,
    edge_color='blue', edge_alpha=0.3, edge_tickness=1,
    edge_text_pos=0.3,
    text_font='sans-serif'):

    # create networkx graph
    G=nx.Graph()

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

        # these are different layouts for the network you may try
        # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)

            # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,font_family=text_font)

    if labels is None:
        labels = range(len(graph))

    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, label_pos=edge_text_pos)

                # show graph
    plt.show()



term_to_search ='Plin2'
gen_network_dict ={}
graph = []
term_list = [term_to_search]
term_next = term_to_search
hops_to_go = 3
hops = 0
while term_next and hops_to_go>0:
    neg = True
    hops_to_go -= 1
    hops += 1
    print hops, hops_to_go, term_next
    if term_next not in term_list[0:hops-1]:
        gen_network_dict[term_next]=[]
        for index, row in sig_corrs.iterrows():
            if term_next in index:
                neg = False
                if index[0]==term_next:
                    gen_corr, corr = index[1], row['corr']
                else:
                    gen_corr, corr = index[0], row['corr']
                gen_network_dict[term_next].append((gen_corr, corr))
                if gen_corr not in term_list:
                    term_list.append(gen_corr)
                if hops != term_list.index(gen_corr):
                    graph.append((hops, term_list.index(gen_corr)))
                print gen_corr
    if not neg:
        term_next = term_list[hops]
    if neg:
        print term_to_search+' not correlated.'
        print gen_network_dict
        term_list = False
print gen_network_dict
print graph
draw_graph(graph, labels=term_list)
corr_by_gene_pos.close()
corr_by_gene_neg.close()
