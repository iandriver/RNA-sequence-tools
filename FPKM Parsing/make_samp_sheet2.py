import fnmatch
import os
import csv


samp_dict = {}
samp_dict['sample_name'] =[]
samp_dict['group'] = []
pats = ['/netapp/home/idriver/results_pdgfra1_ctrl_pnx']
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, '*.cxb'):
            g_cell_name = (root.split('/')[-1])
            if len(g_cell_name.split('_')[0]) == 2:
              cell_name = 'C0'+g_cell_name.split('_')[0][-1]
            else:
              cell_name = g_cell_name.split('_')[0]
            group_name = g_cell_name.split('_')[1]+'_'+cell_name
            samp_path = os.path.join(root,filename)
            samp_dict['sample_name'].append(samp_path)
            samp_dict['group'].append(group_name)
keys = sorted(samp_dict.keys(), reverse=True)
with open("/netapp/home/idriver/pdgfra_sample_sheet.txt", "wb") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(keys)
    writer.writerows(zip(*[samp_dict[key] for key in keys]))
