import os
import cPickle as pickle
import numpy as np
import pandas as pd
from collections import OrderedDict

def delete_cells(by_cell, cell_list, del_list):
  to_delete1 =[]
  for pos, cell_name in enumerate(cell_list):
    if cell_name in del_list:
      to_delete1.append(pos)
  to_delete = sorted(to_delete1, reverse=True)
  for pos in to_delete:
    print pos
    print 'Deleted specific cell '+cell_list[pos]
    del cell_list[pos]
  n_by_cell = np.delete(by_cell, to_delete, axis=0)
  return cell_list, n_by_cell


def filter_by_mapping(path_to_align, cutoff_per_map = 50000):
  c_to_del =[]
  with open(path_to_align, 'rb') as fp:
      a_data = pickle.load(fp)

      p_mapped = a_data['mapped_num']
      ind_list = p_mapped[p_mapped<cutoff_per_map]
      c_to_del = ind_list.index.values
  return c_to_del

def filter_cells_sd(by_cell, cell_list, sd=1.0):
  average_gene_exp = []

  for genes in by_cell:
    gen_exp = (genes >= 0).sum()
    average_gene_exp.append(gen_exp)
  np_av = np.array(average_gene_exp)
  averg = np.average(np_av)
  gene_sd = np.std(np_av)
  print averg, gene_sd
  to_delete= []
  for i1, exp_level in enumerate(np_av):
    if exp_level < (averg - (gene_sd*sd)) or exp_level > (averg + (gene_sd*sd)):
      to_delete.append(i1)
  to_delete1 = sorted(to_delete, reverse = True)
  print to_delete1
  for pos in to_delete1:
    print 'Deleted outlier '+cell_list[pos]
    del cell_list[pos]
  n_by_cell = np.delete(by_cell, to_delete1, axis=0)
  print "Number of cells remaining: "+str(len(cell_list))
  naverage_gene_exp = []
  for ngenes in n_by_cell:
    ngen_exp = (ngenes >= 0).sum()
    naverage_gene_exp.append(ngen_exp)
  nnp_av = np.array(naverage_gene_exp)
  naverg = np.average(nnp_av)
  ngene_sd = np.std(nnp_av)
  print "New", naverg, ngene_sd
  return cell_list, n_by_cell


def threshold_genes(by_gene, gen_list, number_expressed=2):
  g_todelete = []
  for g1, gene in enumerate(by_gene):
    cells_exp = (gene >= 0.5).sum()
    if cells_exp < number_expressed:
      g_todelete.append(g1)
  g1_todelete = sorted(g_todelete, reverse = True)
  print by_gene.shape
  for pos in g1_todelete:
    print 'Gene '+gen_list[pos]+' not expressed in '+str(number_expressed)+' cells.'
    del gen_list[pos]
  n_by_gene = np.delete(by_gene, g1_todelete, axis=0)
  print n_by_gene.shape
  return gen_list, n_by_gene

#given a pandas dataframe of gene expression split out ERCC and return seperate dataframes of each
def sep_ERCC(pd_by_gene, gen_list):
    ERCC_pos_list = []
    ERCC_list= []
    for i, gen in enumerate(gen_list):
        if 'ERCC' in gen:
            ERCC_list.append(gen_list.pop(i))
    pd_by_gene_no_ERCC = pd_by_gene[gen_list]
    pd_ERCC = pd_by_gene[ERCC_list]
    return pd_by_gene_no_ERCC.transpose(), pd_ERCC.transpose(), gen_list

path_to_file= '/Volumes/Seq_data/counts_sheppard_all/cpm_norm_liver_kidney_1'
file_name = 'sheppard_normalized_cpm_liver_kidney.txt'
name = 'liver_kidney'
with open(os.path.join(path_to_file,file_name.strip('.txt')+'.p'), 'rb') as fp:
  data = pickle.load(fp)
  gen_list = data.index.tolist()
  cell_list = [x.strip('_0') for x in list(data.columns.values)]
  path_to_align=os.path.join(path_to_file,'results_liver_all_align.p')
  del_list=filter_by_mapping(path_to_align)

  npdata = np.array(data.values, dtype='f')
  by_cell1 = npdata.transpose()
  rem_cell_list, rem_by_cell = delete_cells(by_cell1, cell_list, del_list)
  npdata2 = rem_by_cell.transpose()
  new_gene_list1, new_by_gene = threshold_genes(npdata2, gen_list)
  by_cell = new_by_gene.transpose()
  outlier_cell_list, outlier_by_cell = filter_cells_sd(by_cell, rem_cell_list)
  final_by_gene = outlier_by_cell.transpose()
  outlier_fpkm_dict = {}
  bulk_ctrl_dict = {}
  filter_on_lane = False
  if filter_on_lane:
    for i, l in enumerate(outlier_by_cell):
        split_cell_list = outlier_cell_list[i].split('_')
        if split_cell_list[0] == 'Lane5' and 'C' in split_cell_list[1] and split_cell_list[1] != 'pdgfra' or split_cell_list[0] == 'Lane6' and 'C' in split_cell_list[1] and split_cell_list[1] != 'pdgfra':
            cell_name = 'CTRL_'+split_cell_list[1]
            outlier_fpkm_dict[cell_name] = [float(lx) for lx in l]
        elif split_cell_list[0] == 'Lane7' and 'C' in split_cell_list[1] and split_cell_list[1] != 'pdgfra' or split_cell_list[0] == 'Lane8'and 'C' in split_cell_list[1] and split_cell_list[1] != 'pdgfra':
            cell_name = 'PNX_'+split_cell_list[1]
            outlier_fpkm_dict[cell_name] = [float(lx) for lx in l]
        else:
            bulk_ctrl_dict['_'.join(split_cell_list[1:])] = [float(lix) for lix in l]
  filter_controls = False
  if filter_controls:
      for i, l in enumerate(outlier_by_cell):
          cell_name = outlier_cell_list[i]
          if 'neg' in cell_name or '+' in cell_name or '-' in cell_name:
              bulk_ctrl_dict[cell_name] = [float(lix) for lix in l]
          else:
              outlier_fpkm_dict[cell_name] = [float(lx) for lx in l]
  for i, l in enumerate(outlier_by_cell):
    outlier_fpkm_dict[outlier_cell_list[i]] = l
  fpkm_df_outlier = pd.DataFrame(outlier_fpkm_dict, index = new_gene_list1)
  fpkm_df_outlier.to_csv(os.path.join(path_to_file, 'norm_cpm_outlier_filtered.txt'), sep = '\t')

  with open(os.path.join(path_to_file,name+'_norm_cpm_outlier_by_gene.p'), 'wb') as fp1:
    pickle.dump(pd.DataFrame(outlier_by_cell.transpose(), columns=outlier_cell_list, index=new_gene_list1), fp1)
  with open(os.path.join(path_to_file,name+'_norm_cpm_outlier_cell_list.p'), 'wb') as fp2:
    pickle.dump(outlier_cell_list, fp2)
  with open(os.path.join(path_to_file,name+'_norm_cpm_outlier_by_cell.p'), 'wb') as fp3:
    pickle.dump(pd.DataFrame(final_by_gene.transpose(), columns=new_gene_list1, index=outlier_cell_list), fp3)
  with open(os.path.join(path_to_file,name+'_norm_cpm_outlier_gene_list.p'), 'wb') as fp4:
    pickle.dump(new_gene_list1, fp4)
