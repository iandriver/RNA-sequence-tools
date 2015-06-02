#/usr/bin/env python
import commands
import os
from subprocess import call


#This is for adjusting parameters on your local system to generate proper file names for qsub
path = '/Volumes/Seq_data/chapmanh'
out= '${TMPDIR}'
annotation_file = '/netapp/home/idriver/genes_E_RS.gtf'
index_gen_loc = '/netapp/home/idriver/mm10_ERCC_RS_bt2/mm10_ERCC_RS/mm10_ERCC_RS'

#this next section parses the file names so that the paired end reads are in order and determines the name of the output file
#use test_qsub.py to test and modify this section locally to work for your file names
pathlist = []
for root, dirs, files in os.walk(path):
    if root.split('/')[-1]=='chapmanh':
        for lane in dirs:
            samp_file_name = lane.split('_')[-1].split('-')
            result_name = '_'.join(['results',samp_file_name[2],samp_file_name[3],samp_file_name[4]])
            print result_name
            #call('mkdir -p /netapp/home/idriver/%s' % result_name, shell=True)
    elif dirs == ['fastqc']:
        n = root.strip('/').split('/')
        out= '${TMPDIR}'
        name = '_'.join([n[-1].split('-')[1],n[-1].split('-')[2],n[-1].split('-')[-1]])
        print name
        data_file = root
        result_file = os.path.join(out,name)
        input_files=''
        r_num = []
        for f in files:
            if 'fastq' in f and ".txt" not in f:
                f_split = f.split('_')
                r_name = (f_split[3][1])
                en_split = f_split[4].split('.')
                p_num = en_split[0].strip('00')
                rank = r_name+p_num
                r_num.append(int(rank))
                input_files+=os.path.join(root,f)+' '
            in_split = input_files.split(' ')
            sort_num = [x for (y,x) in sorted(zip(r_num,in_split))]
        if len(in_split) > 2:
            name_build = ''
            for i, mul_f in enumerate(sort_num):
              if 'fastq' in mul_f:
                if i == len(in_split)-1:
                  name_build+=mul_f
                elif i < (len(in_split)/2)-1 or i > (len(in_split)/2)-1:
                  name_build+= mul_f+','
                elif i ==  (len(in_split)/2)-1:
                  name_build+= mul_f+' '
            final_files = name_build.strip(',')
        elif len(in_split) == 2:
            try:
                final_files = sort_num[0]+' '+sort_num[1].strip(',')
            except IndexError:
                print 'Incomplete File: '+name
        tophat_cmd = 'tophat2 -p 8 -r 230 -a 30 --read-realign-edit-dist 0 -G '+annotation_file+' --transcriptome-index=/netapp/home/idriver/transcriptome_data_mm10_RS/known_e_RS -o '+result_file+' '+index_gen_loc+' '+final_files
        samtools_cmd = 'samtools sort '+result_file+'/'+'accepted_hits.bam accepted_hits_sorted'
        cufflinks_cmd = 'cufflinks -p 8 --max-bundle-frags 10000000 -G '+annotation_file+' -o '+result_file+' '+result_file+'/'+'accepted_hits.bam'
        cuffquant_cmd = 'cuffquant -p 8 --max-bundle-frags 10000000 -o '+result_file+' '+annotation_file+' '+result_file+'/'+'accepted_hits.bam'
