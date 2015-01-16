import fnmatch
import os
import cPickle as pickle


fpkm_matrix_dict_g ={}

g = 0
for root, dirnames, filenames in os.walk('/Users/idriver/RockLab-files/test'):
  for filename in fnmatch.filter(filenames, '*.fpkm_tracking'):
    if '_genes' in filename:
      curr_cell_fpkm_g =[]
      curr_g_file = open(os.path.join(root, filename))
      g_cell_name = (root.split('/')[-1])

      if g == 0:
        curr_cell_genes = []
      for k, line in enumerate(curr_g_file):
        if k == 0:
          header = line.strip('\n').split('\t')
        if k > 0:
          curr_g_line = line.strip('\n').split('\t')
          #exclude RNA spike controls
          if g == 0 and curr_g_line[3][:4] != 'ERCC':
            curr_cell_genes.append(curr_g_line[0])
          if curr_g_line[3][:4] != 'ERCC':
            curr_cell_fpkm_g.append(curr_g_line[9])

      fpkm_matrix_dict_g[g_cell_name] = curr_cell_fpkm_g
      g=1
fpkm_matrix_dict_g['genes'] = curr_cell_genes

with open('spc_fpkm.p', 'wb') as fp:
  pickle.dump(fpkm_matrix_dict_g, fp)
