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
align_dict['input_num'] = []
align_dict['mapped_num'] = []
align_dict['mult_num'] = []
align_dict['per_mapped'] = []
align_dict['mult_mapped_per'] = []
for rf in result_file_names:
    path_to_file = os.path.join(path, rf)
    for root, dirnames, filenames in os.walk(path_to_file):
      for filename in fnmatch.filter(filenames, 'align_summary.txt'):
        g_cell_name = root.split('/')[-1]
        num = g_cell_name.split('_')[-1]
        if len(num) == 2:
          cell_name = '_'.join(g_cell_name.split('_')[0:-1])+'_'+'0'+num
        elif len(num) == 1:
          cell_name = '_'.join(g_cell_name.split('_')[0:-1])+'_'+'00'+num
        else:
          cell_name = g_cell_name
        cell_list.append(cell_name)
        f = open(os.path.join(root,'align_summary.txt'), 'rw')
        for l in f:
            if "Input" in l:
                input_num = int(l.split(':')[-1])
            if "Mapped" in l:
                mapped_1 = l.split(':')[-1]
                mapped_num = int(mapped_1.split('(')[0].strip())
                per_mapped_1 = mapped_1.split('(')[1]
                per_mapped = per_mapped_1.split('%')[0]
            if "of these" in l:
                mult_1 = l.split(':')[-1]
                mult_num = int(mult_1.split('(')[0].strip())
                mult_per_1 = mult_1.split('(')[1]
                mult_per = mult_per_1.split('%')[0]


        align_dict['input_num'].append(input_num)
        align_dict['mapped_num'].append(mapped_num)
        align_dict['mult_num'].append(mult_num)
        align_dict['per_mapped'].append(per_mapped)
        align_dict['mult_mapped_per'].append(mult_per)
        f.close()
align_df = pd.DataFrame(align_dict, index = cell_list)
align_df.to_csv(os.path.join(path,'counts_sheppard_all','results_sheppard_all_align.txt'), sep = '\t')

plt.hist(align_df['mapped_num'])
plt.show()

with open(os.path.join(path,'counts_sheppard_all','results_sheppard_all_align.p'), 'wb') as fp:
  pickle.dump(align_df, fp)
