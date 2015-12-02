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
path_to_file = '/Volumes/Seq_data/cuffnorm_pdgfra_1_and_2'
#name of file containing gene
gene_file_source = 'go_search_genes_lung_all.txt'

base_name = 'pdgfra2_all_n2'
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
    loading_data = pd.read_csv(os.path.join(path_to_file, 'Pdgfra_cell_loading_all.txt'), delimiter= '\t', index_col=0)
    l_data = loading_data.transpose()
    cell_list = gmatrix_df.index.tolist()
    cell_data = []
    cell_label_dict ={'pnx1':('pnx', 4, 'high', 'd4_Pdgfra_ctrl'), 'ctrl1':('ctrl',0, 'high', 'd0_Pdgfra_ctrl'), 'Low_pnx':('pnx',4, 'low', 'd4_PDGFRalpha_low'), 'Low_ctrl':('ctrl', 0, 'ctrl', 'd0_PDGFRalpha_low')}
    new_cell_list = []
    old_cell_list = []
    single_cell_list = []
    single_cell_rename = []
    for cell in cell_list:
        match = False
        if cell[0:4] == 'pnx1':
            k = 'pnx1'
            tracking_id = 'd4_'+cell.split('_')[1]
            match = True
            num = int(cell.split('_')[1].strip('C'))
            if cell.split('_')[1][1]== '0':
                no_z = cell.split('_')[1][0]+cell.split('_')[1][2]
                align_k = no_z+'_'+'Pdgfra-pnxd4'
            else:
                align_k = cell.split('_')[1]+'_'+'Pdgfra-pnxd4'
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:5] == 'ctrl1':
            k = 'ctrl1'
            tracking_id = 'd0_'+cell.split('_')[1]
            num = int(cell.split('_')[1].strip('C'))
            if cell.split('_')[1][1]== '0':
                no_z = cell.split('_')[1][0]+cell.split('_')[1][2]
                align_k = no_z+'_'+'Pdgfra-ctrl1'
            else:
                align_k = cell.split('_')[1]+'_'+'Pdgfra-ctrl1'
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:7] == 'Low_pnx':
            k = 'Low_pnx'
            tracking_id = 'd4_low_'+cell.split('_')[2]
            num = int(cell.split('_')[2].strip('C'))
            if cell.split('_')[2][1]== '0':
                no_z = cell.split('_')[2][0]+cell.split('_')[2][2]
                align_k = 'pdgfra_low_d4pnx_'+no_z
            else:
                align_k = 'pdgfra_low_d4pnx_'+cell.split('_')[2]
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        elif cell[0:8] == 'Low_ctrl':
            k = 'Low_ctrl'
            tracking_id = 'd0_low_'+cell.split('_')[2]
            num = int(cell.split('_')[2].strip('C'))
            if cell.split('_')[2][1]== '0':
                no_z = cell.split('_')[2][0]+cell.split('_')[2][2]
                align_k = 'pdgfra_low_ctrl_'+no_z
            else:
                align_k = 'pdgfra_low_ctrl_'+cell.split('_')[2]
            match = True
            old_cell_list.append(cell)
            new_cell_list.append(tracking_id)
        if match:
            tcf21_level = cmatrix_df[cell]['Tcf21']
            if int(tcf21_level) >= 50:
                Tcf21='high'
            elif int(tcf21_level) < 50 and int(tcf21_level) >= 5 :
                Tcf21='med'
            elif int(tcf21_level) < 5:
                Tcf21='low'
            day = cell_label_dict[k][1]
            condition = cell_label_dict[k][0]
            hi_low = cell_label_dict[k][2]
            loading_df = loading_data[cell_label_dict[k][3]]
            loading = loading_df.iloc[num-1]
            print num, tracking_id
            if loading == '1':
                single_cell = 'yes'
                single_cell_list.append(cell)
                single_cell_rename.append(tracking_id)
            else:
                single_cell = 'no'
            total_mass = by_sample[cell+'_0'][1]
            input_mass = by_cell_map[align_k][0]
            per_mapped = by_cell_map[align_k][4]
            c_data_tup = (tracking_id,total_mass,input_mass,per_mapped,condition,day,hi_low, Tcf21, single_cell)
            print c_data_tup
            cell_data.append(c_data_tup)
    singlecell_cmatrix_df = cmatrix_df[single_cell_list]
    singlecell_cmatrix_df.columns = single_cell_rename
    singlecell_cmatrix_df.to_csv(os.path.join(path_to_file, 'single_cell_matrix.txt'), sep = '\t', index_col=0)
    score_df.to_csv(os.path.join(path_to_file, 'gene_feature_data.txt'), sep = '\t', index=False)
    new_cmatrix_df = cmatrix_df[old_cell_list]
    new_cmatrix_df.columns = new_cell_list
    new_cmatrix_df.to_csv(os.path.join(path_to_file, 'goterms_monocle_count_matrix.txt'), sep = '\t', index_col=0)
    cell_data_df = pd.DataFrame(cell_data, columns=['tracking_id','total_mass','input_mass','per_mapped','condition','day', 'hi_low', 'Tcf21','single_cell'])
    cell_data_df.to_csv(os.path.join(path_to_file, 'cell_feature_data.txt'), sep = '\t', index=False)

make_new_matrix(df_by_gene1, gene_file_source)
