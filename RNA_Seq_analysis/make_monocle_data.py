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
path_to_file ='/Volumes/Seq_data/cuffnorm_Spc2_all_RS'
#name of file containing gene
gene_file_source = 'go_search_genes_lung_all.txt'


path_to_file ='/Volumes/Seq_data/cuffnorm_Spc2_all_RS'
start_file_name = 'fpkm_cuff_spc2'
fpbcell = open(os.path.join(path_to_file, start_file_name+'_outlier_by_cell.p'), 'rb')
by_cell = pickle.load(fpbcell)
fpbcell.close()
fpcelllist = open(os.path.join(path_to_file, start_file_name+'_outlier_cell_list.p'), 'rb')
cell_list = pickle.load(fpcelllist)
fpcelllist.close()
fpbgene = open(os.path.join(path_to_file, start_file_name+'_outlier_by_gene.p'), 'rb')
by_gene = pickle.load(fpbgene)
fpbgene.close()
fpgenelist = open(os.path.join(path_to_file, start_file_name+'_outlier_gene_list.p'), 'rb')
gene_list = pickle.load(fpgenelist)
fpgenelist.close()


df_by_gene = pd.DataFrame(by_cell, columns=gene_list, index=cell_list)
df_by_cell = pd.DataFrame(by_gene, columns=cell_list, index=gene_list)

def make_new_matrix(org_matrix_by_cell, gene_list_file):
    gene_df = pd.read_csv(os.path.join(path_to_file, gene_list_file), delimiter= '\t')
    gene_list = gene_df['GeneID'].tolist()
    group_list = gene_df['GroupID'].tolist()
    new_gmatrix_df = org_matrix_by_cell[gene_list]
    new_cmatrix_df = new_gmatrix_df.transpose()
    new_cmatrix_df.to_csv(os.path.join(path_to_file, 'goterms_monocle_count_matrix.txt'), sep = '\t', index_col=0)
    split_on='_'
    pos=0
    cat_name=['PNX', 'CTRL']
    score_tup =[]
    for gene, go_term in zip(gene_list, group_list):
        sorted_df = new_gmatrix_df.sort([gene])
        score_on = 'PNX'
        g_df = sorted_df[gene]
        ranked_cells = sorted_df.index.values
        ranked_cat = [x.split(split_on)[pos] for x in ranked_cells]
        div_by = int(len(ranked_cat)/len(cat_name))
        start = div_by *(len(cat_name)-1)
        score1 = len([x for x in ranked_cat[start:len(ranked_cat)] if x == score_on])
        tot = len([x for x in ranked_cat if x == score_on])
        res_score = float(score1)/float(tot)
        score = "%.3f" % res_score
        score_tup.append((gene, float(score), go_term))
    score_df = pd.DataFrame(score_tup, columns=['GeneID', 'Rankscore', 'GroupID'])
    score_df.to_csv(os.path.join(path_to_file, 'gene_feature_data.txt'), sep = '\t', index=False)
    sample_data = pd.read_csv(os.path.join(path_to_file, 'samples.table'), delimiter= '\t', index_col=0)
    by_sample = sample_data.transpose()
    map_data = pd.read_csv(os.path.join(path_to_file, 'results_spc2_all_align.txt'), delimiter= '\t', index_col=0)
    by_cell_map = map_data.transpose()
    loading_data = pd.read_csv(os.path.join(path_to_file, 'Cell_loading_SPC2_ctrl_pnx.txt'), delimiter= '\t', index_col=0)
    l_data = loading_data.transpose()
    new_cell_list = new_gmatrix_df.index.values
    cell_data = []
    for cell in new_cell_list:
        condition = cell.split(split_on)[pos]
        if condition == 'PNX':
            day = 4
        else:
            day = 0
        for s_cell in by_sample.columns.values:

            if cell.split(split_on)[1] == s_cell.split('_')[1]:
                tracking_id = cell
                total_mass = by_sample[s_cell][0]
                for m_cell in by_cell_map.columns.values:
                    if cell.split(split_on)[1] in m_cell:
                        print cell, m_cell, '2'
                        print by_cell_map[m_cell]
                        input_mass = by_cell_map[m_cell][0]
                        per_mapped = by_cell_map[m_cell][4]
                        for c_num in l_data.columns.values:
                            if int(c_num) == int(cell.split(split_on)[1].strip('C')):
                                if condition == 'PNX':
                                    loading = l_data[c_num][1]
                                else:
                                    loading = l_data[c_num][0]
                                if loading == '1':
                                    single_cell = 'yes'
                                else:
                                    single_cell = 'no'
                                c_data_tup = (tracking_id,total_mass,input_mass,per_mapped,condition,day,single_cell)
                                print c_data_tup
                                cell_data.append(c_data_tup)
    cell_data_df = pd.DataFrame(cell_data, columns=['tracking_id','total_mass','input_mass','per_mapped','condition','day','single_cell'])
    cell_data_df.to_csv(os.path.join(path_to_file, 'cell_feature_data.txt'), sep = '\t', index=False)

make_new_matrix(df_by_gene, gene_file_source)
