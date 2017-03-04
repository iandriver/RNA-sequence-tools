import fnmatch
import os
import pandas as pd
from scipy import stats
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

path = '/Volumes/Drobo/Seq_data'
result_file_names = ['results_hisat_aged_epi']
basename = 'hisat2_aged_epi'
cell_list =[]
align_dict =OrderedDict()
align_dict['input_reads'] = []
align_dict['num_condcord_0'] = []
align_dict['per_condcord_0'] = []
align_dict['num_condcord_exactly1'] = []
align_dict['per_condcord_exactly1'] = []
align_dict['num_condcord_g1'] = []
align_dict['per_condcord_g1'] = []
align_dict['per_overall'] = []
for rf in result_file_names:
    path_to_file = os.path.join(path, rf)
    for root, dirnames, filenames in os.walk(path_to_file):
        for filename in fnmatch.filter(filenames, '*.o*'):
            check_file = False
            cell_name = (filename.split('.')[0])

            f = open(os.path.join(root,filename), 'r+')
            print(filename)
            check1 = True
            for l in f:
                if 'reads; of these:' in l:
                    total_reads = l.split(' ')[0]
                    check_file =True
                if 'aligned concordantly 0 times' in l and check1:
                    concord0_num = l.split(' ')[4]
                    concord0_per1 = l.split(')')[0]
                    concord0_per2 = concord0_per1.split('(')[-1].strip('%')
                    check1 = False
                if "aligned concordantly exactly 1 time" in l:
                    concord1_num = l.split(' ')[4]
                    concord1_per1 = l.split(')')[0]
                    concord1_per2 = concord1_per1.split('(')[-1].strip('%')
                if "aligned concordantly >1 times" in l:
                    concordg1_num = l.split(' ')[4]
                    concordg1_per1 = l.split(')')[0]
                    concordg1_per2 = concordg1_per1.split('(')[-1].strip('%')
                if "overall alignment rate" in l:
                    overall_per = l.split('%')[0]


            if check_file:
                cell_list.append(cell_name)
                align_dict['input_reads'].append(int(total_reads))
                align_dict['num_condcord_0'].append(int(concord0_num))
                align_dict['per_condcord_0'].append(float(concord0_per2))
                align_dict['num_condcord_exactly1'].append(int(concord1_num))
                align_dict['per_condcord_exactly1'].append(float(concord1_per2))
                align_dict['num_condcord_g1'].append(int(concordg1_num))
                align_dict['per_condcord_g1'].append(float(concordg1_per2))
                align_dict['per_overall'].append(float(overall_per))
                f.close()
align_df = pd.DataFrame(align_dict, index = cell_list)
align_df.to_csv(os.path.join(path,result_file_names[0],'results_'+basename+'_align.txt'), sep = '\t')

plt.hist(align_df['num_condcord_0'])
plt.show()
