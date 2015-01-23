import fnmatch
import os
import pandas as pd
import cPickle as pickle

fpkm_matrix_dict_g ={}
pats = ['/Volumes/Seq_data/results_cindy_rna']
g = 0
coun = 0
for p in pats:
  for root, dirnames, filenames in os.walk(p):
    for filename in fnmatch.filter(filenames, '*genes.fpkm_tracking'):
      if '_genes' in filename:
        curr_cell_fpkm_g =[]
        curr_g_file = open(os.path.join(root, filename))
        if p == '/Volumes/Seq_data/results_cindy_rna':
          g_cell_name = (root.split('/')[-1])
        elif p == '/Volumes/Seq_data/results_hg19_Lib1':
          g_cell_name = (root.split('/')[-1])+'_Lib1_hg19'
          print g_cell_name
          coun+=1
          print coun
        if g == 0:
          curr_cell_genes = []
        for k, line in enumerate(curr_g_file):
          if k == 0:
            header = line.strip('\n').split('\t')
          if k > 0:
            curr_g_line = line.strip('\n').split('\t')
            #exclude RNA spike controls
            if g == 0 and 'mt-' not in curr_g_line[4] and curr_g_line[4][:2] != 'Gm':
              curr_cell_genes.append(curr_g_line[4])
            if 'mt-' not in curr_g_line[4] and curr_g_line[4][:2] != 'Gm':
              curr_cell_fpkm_g.append(curr_g_line[9])
        fpkm_matrix_dict_g[g_cell_name] = curr_cell_fpkm_g
        g=1
fpkm_df_g = pd.DataFrame(fpkm_matrix_dict_g, index = curr_cell_genes)
fpkm_df_g.to_csv('fpkm_matrix_genes_cindy_rna.txt', sep = '\t')

with open('ck_rna_fpkm.p', 'wb') as fp:
  pickle.dump(fpkm_df_g, fp)
