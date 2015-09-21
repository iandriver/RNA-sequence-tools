import os
import shutil
import pandas as pd
from subprocess import call

rename_tup = []
cell_name_done = []
path_to_file= '/Volumes/Seq_data/09142015'

barcode_map1 = pd.read_csv('/Volumes/Seq_data/v2-barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
map1_dict = {}
barcode_map2 = pd.read_csv('/Volumes/Seq_data/human_barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
map2_dict = {}
barcode_dict = {'TAAGGCGA':'N701','CGTACTAG':'N702','AGGCAGAA':'N703', 'TCCTGAGC':'N704',
                'GGACTCCT':'N705', 'TAGGCATG':'N706','CTCTCTAC':'N707', 'CAGAGAGG':'N708',
                'GCTACGCT':'N709', 'CGAGGCTG':'N710', 'AAGAGGCA':'N711', 'GTAGAGGA':'N712',
                'GCTCATGA':'N714', 'ATCTCAGG':'N715',
                'TAGATCGC':'S501', 'CTCTCTAT':'S502', 'TATCCTCT':'S503', 'AGAGTAGA':'S504',
                'GTAAGGAG':'S505', 'ACTGCATA':'S506', 'AAGGAGTA':'S507', 'CTAAGCCT':'S508',
                'CGTCTAAT':'S510', 'TCTCTCCG':'S511', 'TCGACTAG':'S513'}
for n, c in zip(barcode_map1['name'].tolist(),barcode_map1['code'].tolist()):
    map1_dict[c]= n
for n1, c1 in zip(barcode_map2['name'].tolist(),barcode_map2['code'].tolist()):
    map2_dict[c1]= n1
for root, dirs, files in os.walk(path_to_file):
    if root.split('/')[-1][0] == 'C':
        map_to_use = map1_dict
        for f in files:
            if f[-2:] == 'gz':
                name_split = f.split('_')
                code = [name_split[1][0:8],name_split[1][8:]]
                code_lookup = barcode_dict[code[1]]+'/'+barcode_dict[code[0]]
                cell_name = map_to_use[code_lookup]
                if root.split('/')[-2]+cell_name not in cell_name_done:
                    print code
                    cell_name_done.append(root.split('/')[-2]+cell_name)
                    rename_tup.append((root,cell_name))
                    print rename_tup[-1]
                else:
                    pass
for p, name in rename_tup:
    shutil.move(p,os.path.join(os.path.dirname(p),name))
    print p, os.path.join(os.path.dirname(p),name)
