import os
import shutil
import pandas as pd
from subprocess import call

rename_tup = []
cell_name_done = []
path_to_file= '/Volumes/Seq_data/ihg-client.ucsf.edu/chapmanh'
outputdict = {'New_name':[], 'Old_name':[]}
orig_index =[]
barcode_map1 = pd.read_csv('/Volumes/Seq_data/AV_barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
map1_dict = {}
barcode_map2 = pd.read_csv('/Volumes/Seq_data/human_barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
map2_dict = {}
barcode_dict = {'TAAGGCGA':'N701','CGTACTAG':'N702','AGGCAGAA':'N703', 'TCCTGAGC':'N704',
                'GGACTCCT':'N705', 'TAGGCATG':'N706','CTCTCTAC':'N707', 'CAGAGAGG':'N708',
                'GCTACGCT':'N709', 'CGAGGCTG':'N710', 'AAGAGGCA':'N711', 'GTAGAGGA':'N712',
                'GCTCATGA':'N714', 'ATCTCAGG':'N715',
                'TAGATCGC':'S501', 'CTCTCTAT':'S502', 'TATCCTCT':'S503', 'AGAGTAGA':'S504',
                'GTAAGGAG':'S505', 'ACTGCATA':'S506', 'AAGGAGTA':'S507', 'CTAAGCCT':'S508',
                'CGTCTAAT':'S510', 'TCTCTCCG':'S511', 'TCGACTAG':'S513', 'GCGTAAGA':'S517'}
for n, c in zip(barcode_map1['name'].tolist(),barcode_map1['code'].tolist()):
    map1_dict[c]= n
for n1, c1 in zip(barcode_map2['name'].tolist(),barcode_map2['code'].tolist()):
    map2_dict[c1]= n1
new_name_list = []
old_name_list = []
search_dict = {}
results_to_rename = ['results_chapmanh-m-tooth-mesenchyme', 'results_chapmanh-m-doublecre', 'results_chapmanh-hu-IPF-HTII-280', 'results_chapmanh-mSCLib3', 'results_chapmanh-m-CC10-pos']
for root, dirs, files in os.walk(path_to_file):
    if root.split('/')[-1][0] == 'h' or root.split('/')[-1][0] == 'm':
        map_to_use = map1_dict
        for f in files:
            if f[-2:] == 'gz':
                name_split = f.split('_')
                code = [name_split[1][0:8],name_split[1][9:]]
                code_lookup = barcode_dict[code[1]]+'/'+barcode_dict[code[0]]
                name_start = '_'.join(root.split('/')[-1].split('-')[0:-1])
                new_cell_name = name_start+'_'+map_to_use[code_lookup]
                if new_cell_name not in cell_name_done:
                    outputdict['New_name'].append(new_cell_name)
                    cell_name_done.append(new_cell_name)
            elif f[-12:]== 'checksum.txt':
                print f
                num = f.split('-')[-1].split('_')[0]
                old_name = '_'.join(root.split('/')[-1].split('-')[0:-1])+'_'+num
                outputdict['Old_name'].append(old_name)
for n_name, o_name in zip(outputdict['New_name'], outputdict['Old_name']):
    search_dict[o_name] = n_name
name_df = pd.DataFrame.from_dict(outputdict)
for files in results_to_rename:
    pat = os.path.join('/Volumes/Seq_data', files)
    for root, dirs, files in os.walk(pat):
        if root.split('/')[-1] in name_df['Old_name'].values.tolist():
            n_name = search_dict[root.split('/')[-1]]
            shutil.move(root,os.path.join(os.path.dirname(root),n_name))
            print root, os.path.join(os.path.dirname(root),n_name)
