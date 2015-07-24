import subprocess
import os
import fnmatch


result_file_names = ['01272015']
pats = []
N_index1 = {'TAAGGCGA':'1', 'CGTACTAG':'2', 'AGGCAGAA':'3', 'TCCTGAGC':'4', 'GGACTCCT':'5', 'TAGGCATG':'6', 'CTCTCTAC':'7', 'CAGAGAGG':'8',
            'GCTACGCT':'9', 'CGAGGCTG':'10', 'AAGAGGCA':'11', 'GTAGAGGA':'12'}
S_index1 = {'GCGTAAGA':'A', 'CTCTCTAT':'B', 'TATCCTCT':'C', 'AGAGTAGA':'D', 'GTAAGGAG':'E', 'ACTGCATA':'F', 'AAGGAGTA':'G', 'CTAAGCCT':'H'}
for fname in result_file_names:
    pats.append(os.path.join('/Volumes/Seq_data', fname))
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
            g_cell_title = root.split('/')[-1]
            nindex = N_index1[g_cell_title[:8]]
            sindex = S_index1[g_cell_title[8:]]
            subprocess.call('mv '+root+' '+root+'_'+sindex+nindex, shell=True)

            if 'Sample' in g_cell_title:
                #if p.split('/')[-1][-1] == '1':
                    #subprocess.call('mv '+root+' '+'/'.join(root.split('/')[0:-1])+'/'+'1_lib1_'+g_cell_title.split('_')[-1], shell=True)
                #if p.split('/')[-1][-1] == '2':
                    #subprocess.call('mv '+root+' '+'/'.join(root.split('/')[0:-1])+'/'+'1_lib2_'+g_cell_title.split('_')[-1], shell=True)
                pass
