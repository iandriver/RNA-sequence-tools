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

#the file path where gene list will be and where new list will output
path_to_file = '/Volumes/Seq_data/cuffnorm_ips17_BU3_GFP'
#name of file containing gene
gene_file_source = 'go_search_genes_lung_all.txt'

base_name = 'ips17_BU3_GFP'
#load file gene
by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file, base_name+'_outlier_filtered.txt'), sep='\t')
by_gene = by_cell.transpose()
#create list of genes
gene_list = by_cell.index.tolist()
#create cell list
cell_list = [x for x in list(by_cell.columns.values)]


df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list, index=cell_list)
df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list)

def make_new_matrix(org_matrix_by_cell, gene_list_file):
    split_on='_'
    gene_df = pd.read_csv(os.path.join(path_to_file, gene_list_file), delimiter= '\t')
    gene_list = gene_df['GeneID'].tolist()
    group_list = gene_df['GroupID'].tolist()
    gmatrix_df = org_matrix_by_cell[gene_list]
    cmatrix_df = gmatrix_df.transpose()
    score_df = pd.DataFrame(zip(gene_list, group_list), columns=['GeneID', 'GroupID'])
    sample_data = pd.read_csv(os.path.join(path_to_file, 'samples.table'), delimiter= '\t', index_col=0)
    by_sample = sample_data.transpose()
    map_data = pd.read_csv(os.path.join(path_to_file, 'results_'+base_name+'_align.txt'), delimiter= '\t', index_col=0)
    by_cell_map = map_data.transpose()
    loading_data = pd.read_csv(os.path.join(path_to_file, 'ips17_BU3_GFP_cell_loading.txt'), delimiter= '\t', index_col=0)
    l_data = loading_data.transpose()
    cell_list = gmatrix_df.index.tolist()
    cell_data = []
    cell_label_dict ={'BU3':('BU3'), 'ips17':('ips17')}
    new_cell_list = []
    old_cell_list = []
    for cell in cell_list:
        match = False
        if cell[0:2] == 'BU3':
            k = 'BU3'
            match = True
            tracking_id = cell
            c_split = cell.split('_')
            num = int(c_split[2].strip('C_'))
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:5] == 'ips17':
            k = 'ips17'
            c_split = cell.split('_')
            tracking_id = cell
            num = int(c_split[2].strip('C_'))
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        if match:
            nkx_level = cmatrix_df[cell]['NKX2-1']
            if int(nkx_level) >100:
                nkx='high'
            elif int(nkx_level) >25:
                nkx='med'
            elif int(nkx_level) > 1:
                nkx='low'
            else:
                nkx='none'
            cell_type = k
            loading_df = loading_data[cell_label_dict[k][0]]
            loading = loading_df.iloc[num-1]
            print num, tracking_id
            if loading == '1':
                single_cell = 'yes'
            else:
                single_cell = 'no'
            print by_cell_map[align_name]
            total_mass = by_sample[cell+'_0'][1]
            input_mass = by_cell_map[align_name][0]
            per_mapped = by_cell_map[align_name][4]
            c_data_tup = (tracking_id,total_mass,input_mass,per_mapped,cell_type,nkx,single_cell)
            print c_data_tup
            cell_data.append(c_data_tup)
    score_df.to_csv(os.path.join(path_to_file, 'gene_feature_data.txt'), sep = '\t', index=False)
    new_cmatrix_df = cmatrix_df[old_cell_list]
    new_cmatrix_df.columns = new_cell_list
    new_cmatrix_df.to_csv(os.path.join(path_to_file, 'goterms_monocle_count_matrix.txt'), sep = '\t', index_col=0)
    cell_data_df = pd.DataFrame(cell_data, columns=['tracking_id','total_mass','input_mass','per_mapped','cell_type','nkx','single_cell'])
    cell_data_df.to_csv(os.path.join(path_to_file, 'cell_feature_data.txt'), sep = '\t', index=False)

make_new_matrix(df_by_gene1, gene_file_source)
