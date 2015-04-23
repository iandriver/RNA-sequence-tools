import os
import cPickle as pickle
import numpy as np
import pandas as pd


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


def filter_by_mapping(path_to_align, cutoff_per_map = 50):
  c_to_del =[]
  with open(path_to_align, 'rb') as fp:
      a_data = pickle.load(fp)

      p_mapped = a_data['per_mapped']
      ind_list = p_mapped[p_mapped<cutoff_per_map]
      c_to_del = ind_list.index.values
  return c_to_del

def filter_cells_sd(by_cell, cell_list, sd=1.3):
  average_gene_exp = []

  for genes in by_cell:
    gen_exp = (genes >= 1).sum()
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
    ngen_exp = (ngenes >= 1).sum()
    naverage_gene_exp.append(ngen_exp)
  nnp_av = np.array(naverage_gene_exp)
  naverg = np.average(nnp_av)
  ngene_sd = np.std(nnp_av)
  print "New", naverg, ngene_sd
  return cell_list, n_by_cell


def threshold_genes(by_gene, gen_list, number_expressed=3):
  g_todelete = []
  for g1, gene in enumerate(by_gene):
    cells_exp = (gene >= 1.0).sum()
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
    w_gene_list = list(gen_list)
    ERCC_list= []
    pop_list =[]
    to_del = []
    for i, gen in enumerate(w_gene_list):
        if 'ERCC' in gen:
            pop_list.append(i)
        if '_' in gen:
            to_del.append(i)
    to_del = sorted(to_del, reverse=True)
    for d in to_del:
        del w_gene_list[d]
    pop_list = sorted(pop_list, reverse=True)
    for pos in pop_list:
        ERCC_list.append(w_gene_list.pop(pos))
    pd_by_gene_no_ERCC = pd_by_gene[w_gene_list]
    pd_ERCC = pd_by_gene[ERCC_list]
    return pd_by_gene_no_ERCC.transpose(), pd_ERCC.transpose(), w_gene_list

path_to_file = '/Volumes/Seq_data/Spc2_counts'
base_name ='count_spc2'
with open(os.path.join(path_to_file,'spc_count_raw.p'), 'rb') as fp:
  data = pickle.load(fp)
  gen_list = data.index.tolist()
  cell_list = [x.strip('_0') for x in list(data.columns.values)]
  path_to_align=os.path.join(path_to_file,'results_spc2_all_align.p')
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
  filter_controls = True
  if filter_controls:
      for i, l in enumerate(outlier_by_cell):
          cell_name = outlier_cell_list[i]
          if 'neg' in cell_name or '+' in cell_name or '-' in cell_name:
              bulk_ctrl_dict[cell_name] = [int(lix) for lix in l]
          else:
              outlier_fpkm_dict[cell_name] = [int(lx) for lx in l]
  df_bulk = pd.DataFrame(bulk_ctrl_dict, index = new_gene_list1)
  fpkm_df_outlier1 = pd.DataFrame(outlier_fpkm_dict, index = new_gene_list1)
  mod_gen_list = list(new_gene_list1)
  fpkm_df_outlier, df_ERCC, new_gene_list  = sep_ERCC(fpkm_df_outlier1.transpose(), mod_gen_list)
  mod_gen_list = list(new_gene_list1)
  bulk_ctrl_df, df_bulk_ERCC, bulk_gene_list = sep_ERCC(df_bulk.transpose(), mod_gen_list)
  outlier_cell_list = [x for x in list(fpkm_df_outlier.columns.values)]
  df_ERCC.to_csv(os.path.join(path_to_file, base_name+'_ERCC.txt'), sep = '\t')
  df_bulk_ERCC.to_csv(os.path.join(path_to_file, base_name+'_bulkonly_ERCC.txt'), sep = '\t')
  fpkm_df_outlier.to_csv(os.path.join(path_to_file, base_name+'_outlier_filtered.txt'), sep = '\t')
  bulk_ctrl_df.to_csv(os.path.join(path_to_file, base_name+'_bulk_ctrls.txt'), sep = '\t')
  with open(os.path.join(path_to_file,base_name+'_outlier_by_cell.p'), 'wb') as fp1:
    pickle.dump(fpkm_df_outlier.transpose(), fp1)
  with open(os.path.join(path_to_file,base_name+'_outlier_cell_list.p'), 'wb') as fp2:
    pickle.dump(outlier_cell_list, fp2)
  with open(os.path.join(path_to_file,base_name+'_outlier_by_gene.p'), 'wb') as fp3:
    pickle.dump(fpkm_df_outlier, fp3)
  with open(os.path.join(path_to_file,base_name+'_outlier_gene_list.p'), 'wb') as fp4:
    pickle.dump(new_gene_list, fp4)
