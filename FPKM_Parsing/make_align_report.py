import fnmatch
import os
import pandas as pd
import cPickle as pickle
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

path = '/Volumes/Seq_data'
result_file_names = ['results_kidney_raw', 'results_liver_raw', 'results_lung_raw']
cell_list =[]
align_dict ={}
align_dict['input_L_num'] = []
align_dict['mapped_L_num'] = []
align_dict['input_R_num'] = []
align_dict['mapped_R_num'] = []
align_dict['per_mapped'] = []
for rf in result_file_names:
    path_to_file = os.path.join(path, rf)
    for root, dirnames, filenames in os.walk(path_to_file):
      for filename in fnmatch.filter(filenames, 'align_summary.txt'):
        g_cell_name = root.split('/')[-1]
        if len(g_cell_name.split('_')[-1]) == 2 and g_cell_name.split('_')[-1][0] == 'C':
          cell_name = g_cell_name[0:-1]+'0'+g_cell_name[-1]
        else:
          cell_name = g_cell_name
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

        align_dict['input_L_num'].append(input_L_num)
        align_dict['mapped_L_num'].append(mapped_L_num)
        align_dict['input_R_num'].append(input_R_num)
        align_dict['mapped_R_num'].append(mapped_R_num)
        align_dict['per_mapped'].append(per_mapped)
        f.close()
align_df = pd.DataFrame(align_dict, index = cell_list)
align_df.to_csv(os.path.join(path,'counts_sheppard_all','results_sheppard_all_align.txt'), sep = '\t')

plt.hist(align_df['mapped_L_num'])
plt.show()

with open(os.path.join(path,'counts_sheppard_all','results_sheppard_all_align.p'), 'wb') as fp:
  pickle.dump(align_df, fp)
