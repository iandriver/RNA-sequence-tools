import os
import shutil
import pandas as pd

rename_tup = []
path_to_file= '/Volumes/Seq_data/05202015'
file_name = 'Unknow_Name'
barcode_map1 = pd.read_csv('/Volumes/Seq_data/Spcd7-barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
map1_dict = {}
barcode_map2 = pd.read_csv('/Volumes/Seq_data/human_barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
map2_dict = {}
barcode_dict = {'TAAGGCGA':'N701','CGTACTAG':'N702','AGGCAGAA':'N703', 'TCCTGAGC':'N704',
                'GGACTCCT':'N705', 'TAGGCATG':'N706','CTCTCTAC':'N707', 'CAGAGAGG':'N708',
                'GCTACGCT':'N709', 'CGAGGCTG':'N710', 'AAGAGGCA':'N711', 'GTAGAGGA':'N712',
                'TAGATCGC':'S501', 'CTCTCTAT':'S502', 'TATCCTCT':'S503', 'AGAGTAGA':'S504',
                'GTAAGGAG':'S505', 'ACTGCATA':'S506', 'AAGGAGTA':'S507', 'CTAAGCCT':'S508'}
for n, c in zip(barcode_map1['name'].tolist(),barcode_map1['code'].tolist()):
    map1_dict[c]= n
for n1, c1 in zip(barcode_map2['name'].tolist(),barcode_map2['code'].tolist()):
    map2_dict[c1]= n1
for root, dirs, files in os.walk(path_to_file):
    if root.split('/')[-1] == file_name:
        if root.split('/')[-3] == 'Run754':
            map_to_use = map1_dict
        elif root.split('/')[-3] == 'Run753':
            map_to_use = map2_dict
        for d in dirs:
            code = d.split('_')
            code_lookup = barcode_dict[code[1]]+'/'+barcode_dict[code[0]]
            cell_name = map_to_use[code_lookup]+'_'
            rename_tup.append((os.path.join(root,d),cell_name))
for p, name in rename_tup:
    up_one = '/'.join(p.split('/')[0:-1])
    shutil.move(p,os.path.join(up_one,name))
