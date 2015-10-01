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
path_to_file = '/Volumes/Seq_data/cuffnorm_sca_spc_combined'
#name of file containing gene
gene_file_source = 'go_search_genes_lung_all.txt'

base_name = 'combined_sca_spc'
#load file gene
by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file, base_name+'_outlier_filtered.txt'), sep='\t')
by_gene = by_cell.transpose()
#create list of genes
gene_list = by_cell.index.tolist()
#create cell list
cell_list = [x for x in list(by_cell.columns.values)]


df_by_gene1 = pd.DataFrame(by_gene, columns=gene_list, index=cell_list)
df_by_cell1 = pd.DataFrame(by_cell, columns=cell_list, index=gene_list)

def make_new_matrix(org_matrix_by_cell, gene_list_file, rankscore=['pnx', 'ctrl']):
    split_on='_'
    gene_df = pd.read_csv(os.path.join(path_to_file, gene_list_file), delimiter= '\t')
    gene_list = gene_df['GeneID'].tolist()
    group_list = gene_df['GroupID'].tolist()
    gmatrix_df = org_matrix_by_cell[gene_list]
    cmatrix_df = gmatrix_df.transpose()
    cell_list1 = []
    for cell in cmatrix_df.columns.values:
        if cell.split(split_on)[1] == 'ctrl' or cell.split(split_on)[1] == 'pnx':
            if cell.split(split_on)[2][0] =='C':
                print cell, 'cell'
                cell_list1.append(cell)
    new_cmatrix_df = cmatrix_df[cell_list1]
    new_gmatrix_df = new_cmatrix_df.transpose()
    new_cmatrix_df.to_csv(os.path.join(path_to_file, 'goterms_monocle_count_matrix.txt'), sep = '\t', index_col=0)
    pos=1
    cat_name=['pnx', 'ctrl']
    score_tup =[]
    for gene, go_term in zip(gene_list, group_list):
        sorted_df = new_gmatrix_df.sort([gene])
        score_on = 'pnx'
        g_df = sorted_df[gene]
        ranked_cells = sorted_df.index.values
        ranked_cat = [x.split(split_on)[pos] for x in ranked_cells]
        div_by = int(len(ranked_cat)/len(cat_name))
        start = div_by *(len(cat_name)-1)
        score1 = len([x for x in ranked_cat[start:len(ranked_cat)] if x == score_on])
        tot = len([x for x in ranked_cat if x == score_on])
        res_score = float(score1)/float(tot)
        score = "%.3f" % res_score
        score_tup.append((gene, gene, float(score), go_term))
    score_df = pd.DataFrame(score_tup, columns=['', 'GeneID', 'Rankscore', 'GroupID'])
    score_df.to_csv(os.path.join(path_to_file, 'gene_feature_data.txt'), sep = '\t', index=False)
    sample_data = pd.read_csv(os.path.join(path_to_file, 'samples.table'), delimiter= '\t', index_col=0)
    by_sample = sample_data.transpose()
    map_data = pd.read_csv(os.path.join(path_to_file, 'results_spc_d0_4_7_align.txt'), delimiter= '\t', index_col=0)
    by_cell_map = map_data.transpose()
    loading_data = pd.read_csv(os.path.join(path_to_file, 'Cell_loading_SPC_all.txt'), delimiter= '\t', index_col=0)
    l_data = loading_data.transpose()
    new_cell_list = new_gmatrix_df.index.values
    cell_data = []
    for cell in new_cell_list:

        match_all = False
        day_text = cell.split(split_on)[0]
        tracking_id = cell
        if day_text == 'D0':
            day = 0
            condition = 'ctrl'
        elif day_text =='D4':
            day = 4
            condition = 'pnx'
        elif day_text == 'D7':
            day = 7
            condition = 'pnx'
        for s_cell in by_sample.columns.values:
            if cell.split(split_on)[2] == s_cell.split('_')[2]:
                total_mass = by_sample[s_cell][1]
                for m_cell in by_cell_map.columns.values:
                    if cell.split(split_on)[2] in m_cell:
                        input_mass = by_cell_map[m_cell][0]
                        per_mapped = by_cell_map[m_cell][4]
                        for c_num in l_data.columns.values:
                            if int(c_num) == int(cell.split(split_on)[2].strip('C')):
                                match_all = True
                                if day_text == 'D0':
                                    loading = l_data[c_num][1]
                                elif day_text == 'D4':
                                    loading = l_data[c_num][0]
                                elif day_text == 'D7':
                                    loading = l_data[c_num][2]
                                if loading == '1':
                                    single_cell = 'yes'
                                else:
                                    single_cell = 'no'
        if match_all:
            c_data_tup = (tracking_id,total_mass,input_mass,per_mapped,condition,day,single_cell)
            print c_data_tup
            cell_data.append(c_data_tup)
    cell_data_df = pd.DataFrame(cell_data, columns=['tracking_id','total_mass','input_mass','per_mapped','condition','day','single_cell'])
    cell_data_df.to_csv(os.path.join(path_to_file, 'cell_feature_data.txt'), sep = '\t', index=False)

make_new_matrix(df_by_gene1, gene_file_source)
