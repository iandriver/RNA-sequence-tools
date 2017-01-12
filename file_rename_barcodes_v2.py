import subprocess
import os
import fnmatch
import pandas as pd


result_file_names = ['09142015']
path_to_file = '/Volumes/Drobo/Seq_data'
pats = []
new_filename = 'Renamed_ips17_Bu3'
os.mkdir(os.path.join(path_to_file,new_filename))


barcode_map1 = pd.read_csv('/Volumes/Seq_data/v2-barcodes.txt', sep='\t', header = None, names = ['name', 'code'])
map1_dict = {}
barcode_dict = {'TAAGGCGA':'N701','CGTACTAG':'N702','AGGCAGAA':'N703', 'TCCTGAGC':'N704',
                'GGACTCCT':'N705', 'TAGGCATG':'N706','CTCTCTAC':'N707', 'CAGAGAGG':'N708',
                'GCTACGCT':'N709', 'CGAGGCTG':'N710', 'AAGAGGCA':'N711', 'GTAGAGGA':'N712',
                'GCTCATGA':'N714', 'ATCTCAGG':'N715',
                'TAGATCGC':'S501', 'CTCTCTAT':'S502', 'TATCCTCT':'S503', 'AGAGTAGA':'S504',
                'GTAAGGAG':'S505', 'ACTGCATA':'S506', 'AAGGAGTA':'S507', 'CTAAGCCT':'S508',
                'CGTCTAAT':'S510', 'TCTCTCCG':'S511', 'TCGACTAG':'S513', 'GCGTAAGA':'S517'}

N_index1 = {'N701':'1', 'N702':'2', 'N703':'3', 'N704':'4', 'N705':'5', 'N706':'6', 'N707':'7', 'N710':'8',
            'N711':'9', 'N712':'10', 'N714':'11', 'N715':'12'}
S_index1 = {'S501':'A', 'S502':'B', 'S503':'C', 'S504':'D', 'S505':'E', 'S506':'F', 'S507':'G', 'S508':'H'}

N_index2 = {'N701':'1', 'N702':'2', 'N703':'3', 'N704':'4', 'N705':'5', 'N706':'6', 'N707':'7', 'N710':'8',
            'N711':'9', 'N712':'10', 'N714':'11', 'N715':'12'}
S_index2 = {'S502':'A', 'S503':'B', 'S505':'C', 'S506':'D', 'S507':'E', 'S508':'F', 'S510':'G', 'S511':'H'}

for n, c in zip(barcode_map1['name'].tolist(),barcode_map1['code'].tolist()):
    map1_dict[c]= n

for fname in result_file_names:
    pats.append(os.path.join(path_to_file, fname))
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        x=0
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
            if x==0:
                g_cell_title = filename.split('_')[1]
                barcode1 = g_cell_title[:8]
                barcode2 = g_cell_title[8:]
                code_lookup1 = barcode_dict[barcode1]
                code_lookup2 = barcode_dict[barcode2]
                code_lookup = code_lookup2+'/'+code_lookup1
                skip= False
                try:
                    cell_name = map1_dict[code_lookup]
                    nindex = N_index2[code_lookup1]
                    sindex = S_index2[code_lookup2]
                except KeyError:
                    print(code_lookup, 'file:',root, 'BarcodeError')
                    skip = True
                if not skip:
                    name = cell_name+'_'+sindex+nindex
                    
                    #uncomment when you actually want to move
                    #shutil.move(root,os.path.join(path_to_file,new_filename,name))
                    print(root, os.path.join(path_to_file,new_filename,name))
                    x+=1
