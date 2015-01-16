import fnmatch
import os
import pandas as pd
import cPickle as pickle

fpkm_matrix_dict_g ={}
pats = ['/Volumes/Seq_data/results_pdgfra_mm10']
g = 0
for p in pats:
  for root, dirnames, filenames in os.walk(p):
    for filename in fnmatch.filter(filenames, '*.fpkm_tracking'):
      if '_genes' in filename:
        curr_cell_fpkm_g =[]
        curr_g_file = open(os.path.join(root, filename))
        if p == '/Volumes/Seq_data/results_spc_pnx_mm10':
          g_cell_name = (root.split('/')[-1])+'_pnx'
        elif p == '/Volumes/Seq_data/results_pdgfra_mm10':
          g_cell_name = (root.split('/')[-1])+'_ctrl'
        if g == 0:
          curr_cell_genes = []
        for k, line in enumerate(curr_g_file):
          if k == 0:
            header = line.strip('\n').split('\t')
          if k > 0:
            curr_g_line = line.strip('\n').split('\t')
            #exclude RNA spike controls
            if g == 0:
              curr_cell_genes.append(curr_g_line[4])

            curr_cell_fpkm_g.append(curr_g_line[9])

        fpkm_matrix_dict_g[g_cell_name] = curr_cell_fpkm_g
        g=1
fpkm_df_g = pd.DataFrame(fpkm_matrix_dict_g, index = curr_cell_genes)
fpkm_df_g.to_csv('fpkm_matrix_genes_pdgfra_mm10_combined_nomtorgm.txt', sep = '\t')

with open('spc_pdgfra_fpkm.p', 'wb') as fp:
  pickle.dump(fpkm_df_g, fp)
