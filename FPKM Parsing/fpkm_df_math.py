
import cPickle as pickle
import numpy as np
import pandas as pd


def delete_cells(by_cell, cell_list, del_list):
  to_delete1 =[]
  for pos, num in enumerate(cell_list):
    if num in del_list:
      to_delete1.append(num)

  for pos in to_delete1:
    print 'Deleted specific cell '+cell_list[pos]
    del cell_list[pos]
  n_by_cell = np.delete(by_cell, to_delete1, axis=0)
  return cell_list, n_by_cell


def filter_by_mapping(path_to_align, cutoff_per_map = 50):
    c_to_del =[]
    with open(path_to_align, 'rb') as fp:
      a_data = pickle.load(fp)
      ind_list = a_data.index.tolist()
      cell_list = list(a_data.columns.values)
      for c in cell_list:
          if a_data[c]['per_mapped'] < cutoff_per_map:
              c_to_del.append(c)
    return c_to_del

def filter_cells_sd(by_cell, cell_list, sd=1.2):
  average_gene_exp = []

  for genes in by_cell:
    gen_exp = (genes >= 1.0).sum()
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
    ngen_exp = (ngenes >= 1.0).sum()
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



with open('Lib1and2_hg19_fpkm.p', 'rb') as fp:
  data = pickle.load(fp)
  gen_list = data.index.tolist()
  cell_list = list(data.columns.values)
  path_to_align='/Volumes/Seq_data/pdgfra_align.p'
  del_list=filter_by_mapping(path_to_align)

  npdata = np.array(data.values, dtype='f')
  by_cell1 = npdata.transpose()
  rem_cell_list, rem_by_cell = delete_cells(by_cell1, cell_list, del_list)
  npdata2 = rem_by_cell.transpose()
  new_gene_list, new_by_gene = threshold_genes(npdata2, gen_list)
  by_cell = new_by_gene.transpose()
  outlier_cell_list, outlier_by_cell = filter_cells_sd(by_cell, rem_cell_list)
  final_by_gene = outlier_by_cell.transpose()
  outlier_fpkm_dict = {}
  for i, l in enumerate(outlier_by_cell):
    outlier_fpkm_dict[outlier_cell_list[i]] = l
  fpkm_df_outlier = pd.DataFrame(outlier_fpkm_dict, index = new_gene_list)
  fpkm_df_outlier.to_csv('fpkm_matrix_genes_Lib1_2_hg19_outlier_filtered.txt', sep = '\t')
  with open('lib1_2_outlier_by_cell.p', 'wb') as fp1:
    pickle.dump(outlier_by_cell, fp1)
  with open('lib1_2_outlier_cell_list.p', 'wb') as fp2:
    pickle.dump(outlier_cell_list, fp2)
  with open('lib1_2_outlier_by_gene.p', 'wb') as fp3:
    pickle.dump(final_by_gene, fp3)
  with open('lib1_2_outlier_gene_list.p', 'wb') as fp4:
    pickle.dump(new_gene_list, fp4)
