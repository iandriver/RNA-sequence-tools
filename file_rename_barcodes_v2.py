import subprocess
import os
import fnmatch


result_file_names = ['09142015_BU3_ips17_raw-data']
path_to_file = '/Volumes/Seq_data'
pats = []

barcode_dict = {'TAAGGCGA':'N701','CGTACTAG':'N702','AGGCAGAA':'N703', 'TCCTGAGC':'N704',
                'GGACTCCT':'N705', 'TAGGCATG':'N706','CTCTCTAC':'N707', 'CAGAGAGG':'N708',
                'GCTACGCT':'N709', 'CGAGGCTG':'N710', 'AAGAGGCA':'N711', 'GTAGAGGA':'N712',
                'GCTCATGA':'N714', 'ATCTCAGG':'N715',
                'TAGATCGC':'S501', 'CTCTCTAT':'S502', 'TATCCTCT':'S503', 'AGAGTAGA':'S504',
                'GTAAGGAG':'S505', 'ACTGCATA':'S506', 'AAGGAGTA':'S507', 'CTAAGCCT':'S508',
                'CGTCTAAT':'S510', 'TCTCTCCG':'S511', 'TCGACTAG':'S513', 'GCGTAAGA':'S517'}

N_index1 = {'N701':'1', 'N702':'2', 'N703':'3', 'N704':'4', 'N705':'5', 'N706':'6', 'N707':'7', 'N710':'8',
            'N711':'9', 'N712':'10', 'N714':'11', 'N715':'12'}
S_index1 = {'S502':'A', 'S503':'B', 'S505':'C', 'S506':'D', 'S507':'E', 'S508':'F', 'S510':'G', 'S511':'H'}
for fname in result_file_names:
    pats.append(os.path.join(path_to_file, fname))
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        x=0
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
            if x==0:
                g_cell_title = filename.split('_')[1]
                nindex = N_index1[barcode_dict[g_cell_title[:8]]]
                sindex = S_index1[barcode_dict[g_cell_title[8:]]]
                subprocess.call('mv '+root+' '+root+'_'+sindex+nindex, shell=True)
                print('mv '+root+' '+root+'_'+sindex+nindex)
                x+=1
