import fnmatch
import os
import pandas as pd

fpkm_matrix_dict_g ={}

g = 0
for root, dirnames, filenames in os.walk('/Users/idriver/RockLab-files/finished_spc_ensembl'):
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
          if g == 0 and 'mt-' not in curr_g_line[4] and curr_g_line[4][:2] != 'Gm':
            curr_cell_genes.append(curr_g_line[4]+'_'+curr_g_line[0][-5:])
          if 'mt-' not in curr_g_line[4] and curr_g_line[4][:2] != 'Gm':
            curr_cell_fpkm_g.append(curr_g_line[9])

      fpkm_matrix_dict_g[g_cell_name] = curr_cell_fpkm_g
      g=1
fpkm_df_g = pd.DataFrame(fpkm_matrix_dict_g, index = curr_cell_genes)
fpkm_df_g.to_csv('fpkm_matrix_genes_spc_ensembl_nomtorgm.txt', sep = '\t')
