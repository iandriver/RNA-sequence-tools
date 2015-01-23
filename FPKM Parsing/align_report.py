import fnmatch
import os
import pandas as pd
import cPickle as pickle


align_dict={}
cell_list =[]
for root, dirnames, filenames in os.walk('/Volumes/Seq_data/finished_pdgfra_ensembl'):
  for filename in fnmatch.filter(filenames, 'align_summary.txt'):
    cell_name = root.split('/')[-1]
    cell_list.append(cell_name)
    f = open(os.path.join(root,'align_summary.txt'), 'rw')
    for l in f:
        if 'Left' in l:
            side_s = 0
        elif 'Right' in l:
            side_s = 1
        if "Input" in l and side_s == 0:
            input_L_num = int(l.split(':')[-1])
        if "Mapped" in l and side_s == 0:
            mapped_L_1 = l.split(':')[-1]
            mapped_L_num = int(mapped_L_1.split('(')[0].strip())
        if "Input" in l and side_s == 1:
            input_R_num = int(l.split(':')[-1])
        if "Mapped" in l and side_s == 0:
            mapped_R_1 = l.split(':')[-1]
            mapped_R_num = int(mapped_R_1.split('(')[0].strip())
        if "overall read mapping rate." in l:
            per_mapped = float(l.split('%')[0])

    align_dict[cell_name] = [input_L_num, mapped_L_num, input_R_num, mapped_R_num, per_mapped]
    f.close()
align_df = pd.DataFrame(align_dict, index = ['input_L_num', 'mapped_L_num', 'input_R_num', 'mapped_R_num', 'per_mapped'])
align_df.to_csv('pdgfra_ensmbl_align.txt', sep = '\t')

with open('pdgfra_align.p', 'wb') as fp:
  pickle.dump(align_df, fp)
