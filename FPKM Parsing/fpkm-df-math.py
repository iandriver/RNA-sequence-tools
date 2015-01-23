
import cPickle as pickle
import numpy as np
import pandas as pd



def filter_cells_sd(by_cell, cell_list, sd=1.2):
  average_gene_exp = []

  for genes in by_cell:
    gen_exp = (genes >= 1.0).sum()
    average_gene_exp.append(gen_exp)
  np_av = np.array(average_gene_exp)
  averg = np.average(np_av)
  gene_sd = np.std(np_av)
  print averg, sd
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



with open('pdgfra_fpkm.p', 'rb') as fp:
  data = pickle.load(fp)
  gen_list = data.index.tolist()
  cell_list = list(data.columns.values)


  npdata = np.array(data.values, dtype='f')
  new_gene_list, new_by_gene = threshold_genes(npdata, gen_list)
  by_cell = new_by_gene.transpose()
  outlier_cell_list, outlier_by_cell = filter_cells_sd(by_cell, cell_list)
  final_by_gene = outlier_by_cell.transpose()
  with open('pdgfra_outlier_by_cell.p', 'wb') as fp1:
    pickle.dump(outlier_by_cell, fp1)
  with open('pdgfra_outlier_cell_list.p', 'wb') as fp2:
    pickle.dump(outlier_cell_list, fp2)
  with open('pdgfra_outlier_by_gene.p', 'wb') as fp3:
    pickle.dump(final_by_gene, fp3)
  with open('pdgfra_outlier_gene_list.p', 'wb') as fp4:
    pickle.dump(new_gene_list, fp4)
