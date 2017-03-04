import fnmatch
import os
import pandas as pd
import subprocess

'''This program takes accepted_hits.bam files from tophat and turns them into
counts and creates matrix file of cells/conditions and counts using htseq:
http://www-huber.embl.de/HTSeq/doc/overview.html

The files are sorted using samtools: http://www.htslib.org/

Paired end mates are fixed and RNA metrics collected using Picard tools:
http://broadinstitute.github.io/picard/
The 3' to 5' bias of each sample is collected as a matrix file for easy plotting.
'''
#list of file paths with mapped hits
pats = ['rsem_m38_Cfms-d7-Sham_insert_metrics']
#output path


#initialize dictonaries for collected output
fpkm_matrix_dict_g ={}
count_dict = {}
norm_read_dict = {}
picard_stats_dict = {}

#collect gene_list once since it the same between all samples
st = 1
gene_list = []

#loop through all files and sort, fix, count, collect metrics on each
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, 'insert_size_metrics.txt'):
            with open(os.path.join(root,'insert_size_metrics.txt'), mode='r') as f:
                next_line=-1
                cname=os.path.basename(root)
                for i, l in enumerate(f):
                     if l[0:18] == 'MEDIAN_INSERT_SIZE':
                             titles = l.split('\t')
                             next_line= i+1
                     if i == next_line:
                             metrics = l.split('\t')
                a = dict(zip(titles, metrics))
            picard_stats_dict[cname] = a
final_df = pd.DataFrame(picard_stats_dict)
final_df_t = final_df.transpose()
final_df_t.to_csv(os.path.basename(pats[0])+"_picard_insert_metrics.txt", sep='\t')
