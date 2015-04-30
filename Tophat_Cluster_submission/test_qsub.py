#/usr/bin/env python
import commands
import os
from subprocess import call


#This is for adjusting parameters on your local system to generate proper file names for qsub
path = '/Volumes/C1_data/March2015/raw_data'
out= '${TMPDIR}'
annotation_file = '/netapp/home/idriver/mm10_ERCC/genes/genes.gtf'
index_gen_loc = '/netapp/home/idriver/mm10_ERCC/Bowtie2Index_mm10/mm10_ERCC'
result_file_name = 'results_pdgfra1_ctrl_pnx'

pathlist = []
for root, dirs, files in os.walk(path):
    if root.split('/')[-1]=='raw_data':
        for lane in dirs:
            result_name = 'results_'+lane
            if ' ' in result_name:
                result_s = result_name.split(' ')
                result_name = '_'.join(result_s)
            print result_name
    elif dirs == []:
        pathlist.append([root,files])
for p in pathlist:
    n = p[0].split('/')
    print n, 'n'
    name = n[-1].split('_')[0]
    if ' ' in name:
        name = name.split(' ')[0]
        print name, 'name'
    for f in p[1]:
        if 'fastq' in f and ".txt" not in f:
            f_split = f.split('.')
            f_name = name+'_'+f_split[0]
            final_files = f
            result_file = os.path.join(out,f_name)
            tophat_cmd = 'tophat2 -p 8 -r 50 -G '+annotation_file+' --transcriptome-index=/netapp/home/idriver/transcriptome_data_mm10/known -o '+result_file+' '+index_gen_loc+' '+final_files
            cufflinks_cmd = 'cufflinks -p 8 -G '+annotation_file+' -o '+result_file+' '+result_file+'/'+'accepted_hits.bam'
            cuffquant_cmd = 'cuffquant -p 8 -o '+result_file+' '+annotation_file+' '+result_file+'/'+'accepted_hits.bam'
            # Write script.
