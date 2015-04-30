import fnmatch
import os
import pandas as pd
import cPickle as pickle
import subprocess
import csv
'''This program takes accepted_hits.bam files from tophat and turns them into
counts and creates matrix file of cells/conditions and counts using htseq:
http://www-huber.embl.de/HTSeq/doc/overview.html

The files are sorted using samtools: http://www.htslib.org/

Paired end mates are fixed and RNA metrics collected using Picard tools:
http://broadinstitute.github.io/picard/
The 3' to 5' bias of each sample is collected as a matrix file for easy plotting.
'''
#list of file paths with mapped hits
pats = ['/Volumes/Seq_data/results_kidney_raw', '/Volumes/Seq_data/results_liver_raw', '/Volumes/Seq_data/results_lung_raw']

#output path
path = '/Volumes/Seq_data'

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
        for filename in fnmatch.filter(filenames, 'accepted_hits.bam'):
            #sort accepted_hits.bam using samtools
            cname = root.split('/')[-1]
            sort_out = os.path.join(root,cname+'_sorted')
            sam_sort_calln = 'samtools sort -n '+os.path.join(root,filename)+' '+sort_out
            sam_sort_call = 'samtools sort '+os.path.join(root,filename)+' '+sort_out
            print sam_sort_call
            #skip if file already exists
            if not os.path.isfile(sort_out+'.bam'):
                process = subprocess.Popen(sam_sort_call, stdout=subprocess.PIPE, shell=True)
                out, err = process.communicate()
                print(out)
            else:
                print sort_out+'.bam already exists'
                pass
            #format htseq-count command to generate raw counts from sorted accepted hits
            gf = '/Volumes/Seq_data/genes_E_RS.gtf'
            hts_out = os.path.join(root,cname+'_htseqcount.txt')
            if not os.path.isfile(hts_out):
                htseq_count_call = 'htseq-count -f bam '+sort_out+'.bam'+' '+gf+' > '+hts_out
                print htseq_count_call
                process = subprocess.Popen(htseq_count_call, stdout=subprocess.PIPE, shell=True)
                out, err = process.communicate()
                print(out)
            else:
                print hts_out+' already exists'
                pass
            #run picard_fixmate to clean up paired end reads in accepted_hits.bam (sorted)
            '''picard_fixmate_out = os.path.join(root,sort_out.strip('.bam')+'_FM.bam')
            if not os.path.isfile(picard_fixmate_out):
                picard_fixmate_call = 'java -Xmx3g -jar /Users/idriver/picard/dist/picard.jar FixMateInformation INPUT='+sort_out+'.bam OUTPUT='+picard_fixmate_out+' AS=true SORT_ORDER=coordinate'
                print picard_fixmate_call
                process = subprocess.Popen(picard_fixmate_call, stdout=subprocess.PIPE, shell=True)
                out, err = process.communicate()
                print(out)
            else:
                print picard_fixmate_out+' already exists'
            '''
            #run picard CollectRnaSeqMetrics (http://broadinstitute.github.io/picard/command-line-overview.html) and generate matrix of 3' to 5' bias (norm_read_dict)
            picard_rnaseqmetric_out = os.path.join(root,sort_out.strip('sorted.bam')+'RNA_metric.txt')
            picard_rnaseqchart_out = os.path.join(root,sort_out.strip('sorted.bam')+'RNA_metric.pdf')
            sorted_bam_in = os.path.join(root,sort_out+'.bam')
            if not os.path.isfile(picard_rnaseqchart_out):
                picard_seqmetric_call = 'java -Xmx3g -jar /Users/idriver/picard/dist/picard.jar CollectRnaSeqMetrics REF_FLAT=/Volumes/Seq_data/refFlat.txt.gz STRAND_SPECIFICITY=NONE MINIMUM_LENGTH=70 CHART_OUTPUT='+picard_rnaseqchart_out+' INPUT='+sorted_bam_in+' OUTPUT='+picard_rnaseqmetric_out
                print picard_seqmetric_call
                process = subprocess.Popen(picard_seqmetric_call, stdout=subprocess.PIPE, shell=True)
                out, err = process.communicate()
                print(out)
            else:
                print picard_rnaseqchart_out+' already exists'
            g_counts = []
            with open(hts_out, mode='r') as infile:
                hts_tab = csv.reader(infile, delimiter = '\t')
                print st
                for l in hts_tab:
                    if st == 1:
                        gene_list.append(l[0])
                    g_counts.append(l[1])
                st = 2
                print len(g_counts)
                print len(gene_list)
                count_dict[cname] = g_counts
            norm_read_dict[cname] = []
            index3 = []
            with open(picard_rnaseqmetric_out, mode='r') as infile:
                pic_tab = csv.reader(infile, delimiter = '\t')
                for i, l in enumerate(pic_tab):
                    if i == 6:
                        index1 = l
                    if i == 7:
                        num_stats = []
                        for n in l:
                            if n == '' or n == '?':
                                num_stats.append(0.0)
                            else:
                                num_stats.append(float(n))
                        picard_stats_dict[cname] = num_stats
                    if i == 10:
                        index2 = l
                    if i > 10 and i <= 111:
                        index3.append(int(l[0]))
                        norm_read_dict[cname].append(float(l[1]))
for k, v in norm_read_dict.items():
    if len(v) == 0:
        norm_read_dict[k] = [0 for x in range(101)]
        print norm_read_dict[k], len(norm_read_dict[k])
print index3

#form pandas dataframe of each and save as tab delimited file
count_df = pd.DataFrame(count_dict, index = gene_list)
count_df.to_csv(os.path.join(path,'sheppard_all_count_table.txt'), sep = '\t')
with open(os.path.join(path,'htseq_count_sheppard_all.p'), 'wb') as fp1:
    pickle.dump(count_df, fp1)
pic_stats_df = pd.DataFrame(picard_stats_dict, index = index1)
pic_stats_df.to_csv(os.path.join(path,'sheppard_picard_stats.txt'), sep = '\t')
norm_read_df = pd.DataFrame(norm_read_dict, index = index3)
norm_read_df.to_csv(os.path.join(path,'sheppard_read_bias.txt'), sep = '\t')
pd.DataFrame.plot(norm_read_df)
