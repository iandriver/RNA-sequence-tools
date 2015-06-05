import fnmatch
import os
import pandas as pd
import cPickle as pickle
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import OrderedDict

path = '/Volumes/Seq_data'
result_file_names = ['results_Human_NR']
basename = 'hu_NR'
cell_list =[]
align_dict =OrderedDict()
align_dict['input_L_num'] = []
align_dict['mapped_L_num'] = []
align_dict['input_R_num'] = []
align_dict['mapped_R_num'] = []
align_dict['per_mapped'] = []
for rf in result_file_names:
    path_to_file = os.path.join(path, rf)
    for root, dirnames, filenames in os.walk(path_to_file):
        for filename in fnmatch.filter(filenames, 'align_summary.txt'):
            cell_name = (root.split('/')[-1])
            if cell_name[-1] != '_' and '_' in cell_name:
                if cell_name.split('_')[-1][0] == 'b':
                    group_name = 'bulk'
            elif cell_name[-1] == '_':
                group_name = 'NR_repeat'
            else:
                group_name = 'NR'
            cell_name = group_name+'_'+cell_name.strip('C_')
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
align_df.to_csv(os.path.join(path,result_file_names[0],'results_'+basename+'_align.txt'), sep = '\t')

plt.hist(align_df['mapped_L_num'])
plt.show()

with open(os.path.join(path,result_file_names[0],'results_'+basename+'_align.p'), 'wb') as fp:
  pickle.dump(align_df, fp)
