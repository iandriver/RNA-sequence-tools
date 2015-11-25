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
path_to_file = '/Users/idriver/RockLab-files/RNA-sequence-tools/Gene_Ontology/gene_lookup.py'
#name of file containing gene
gene_file_source = 'go_search_genes_lung_all.txt'
plate_map = 'plate_map_grid.txt'

base_name = 'hu_alpha6_ht280_combined'
#load file gene
by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file, base_name+'_outlier_filtered.txt'), sep='\t')
by_gene = by_cell.transpose()
#create list of genes
gene_list = by_cell.index.tolist()
#create cell list
cell_list = [x for x in list(by_cell.columns.values)]

plate_map_df = pd.DataFrame.from_csv(os.path.join(path_to_file, plate_map), sep='\t')
def ret_loading_pos(pos, plate_map_df):
    let_dict = {'A':0, 'B':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7}
    let = pos[0]
    num = pos[1:]
    num_col = plate_map_df[num]
    cnum = num_col.iloc[let_dict[let]]
    return int(cnum.strip('C'))

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
    loading_data = pd.read_csv(os.path.join(path_to_file, 'Cell_loading_hu_ht280_alpha6_all.txt'), delimiter= '\t', index_col=0)
    l_data = loading_data.transpose()
    cell_list = gmatrix_df.index.tolist()
    cell_data = []
    cell_label_dict ={'norm':('norm_ht280', 'ctrl', 'norm'),
                    'scler':('scler_ht280', 'diseased', 'scler'),
                    'IPF':('hu_IPF_HTII_280','diseased', 'IPF'),
                    'DK':('DK_ht280','diseased', 'DK'),
                    'alpha6_norm':('norm_alpha6', 'ctrl', 'norm'),
                    'alpha6_scler':('scler_alpha6', 'diseased', 'scler')}
    new_cell_list = []
    old_cell_list = []
    for cell in cell_list:
        match = False
        if cell[0:6] == 'norm_h':
            k = 'norm'
            tracking_id = cell
            match = True
            num = cell.split('_')[2]
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:7] == 'scler_h':
            k='scler'
            tracking_id = cell
            num = cell.split('_')[2]
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:2] == 'hu':
            k = 'IPF'
            tracking_id = cell
            num = cell.split('_')[4]
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:2] == 'DK':
            k = 'DK'
            tracking_id = cell
            num = cell.split('_')[2]
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:7] == 'scler_a':
            k = 'alpha6_scler'
            tracking_id = cell
            num = cell.split('_')[2]
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:6] == 'norm_a':
            k = 'alpha6_norm'
            tracking_id = cell
            num = cell.split('_')[2]
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        if match:
            condition = cell_label_dict[k][1]
            disease = cell_label_dict[k][2]
            loading_df = loading_data[cell_label_dict[k][0]]
            print ret_loading_pos(num, plate_map_df), num
            loading = loading_df.iloc[ret_loading_pos(num, plate_map_df)-1]
            print num, tracking_id
            if 'single' in loading:
                single_cell = 'yes'
            else:
                single_cell = 'no'
            print by_cell_map[cell]
            total_mass = by_sample[cell+'_0'][1]
            input_mass = by_cell_map[cell][0]
            per_mapped = by_cell_map[cell][4]
            c_data_tup = (tracking_id,total_mass,input_mass,per_mapped,condition,disease,single_cell)
            print c_data_tup
            cell_data.append(c_data_tup)
    score_df.to_csv(os.path.join(path_to_file, 'gene_feature_data.txt'), sep = '\t', index=False)
    new_cmatrix_df = cmatrix_df[old_cell_list]
    new_cmatrix_df.columns = new_cell_list
    new_cmatrix_df.to_csv(os.path.join(path_to_file, 'goterms_monocle_count_matrix.txt'), sep = '\t', index_col=0)
    cell_data_df = pd.DataFrame(cell_data, columns=['tracking_id','total_mass','input_mass','per_mapped','condition','disease','single_cell'])
    cell_data_df.to_csv(os.path.join(path_to_file, 'cell_feature_data.txt'), sep = '\t', index=False)

make_new_matrix(df_by_gene1, gene_file_source)
