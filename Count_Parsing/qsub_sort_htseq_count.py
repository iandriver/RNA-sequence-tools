import fnmatch
import os
import pandas as pd
import cPickle as pickle
import subprocess
import csv
from collections import OrderedDict

#list of file paths with mapped hits
pats = ['/Volumes/Seq_data/results_01272015', '/Volumes/Seq_data/results_spc2_n2']
#output path
path = '/Volumes/Seq_data'
base_name = 'combined_spc'

#initialize dictonaries for collected output
fpkm_matrix_dict_g =OrderedDict()
count_dict = OrderedDict()
norm_read_dict = OrderedDict()
picard_stats_dict = OrderedDict()

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

            #run picard_fixmate to clean up paired end reads in accepted_hits.bam (sorted)
            picard_fixmate_out = os.path.join(root,sort_out.strip('.bam')+'_FM.bam')
            picard_fixmate_call = 'java -Xmx3g -jar /netapp/home/idriver/broadinstitute-picard-fd8e773/dist/picard.jar FixMateInformation INPUT='+sort_out+'.bam OUTPUT='+picard_fixmate_out+' AS=true SORT_ORDER=coordinate'

            #format htseq-count command to generate raw counts from sorted accepted hits
            gf = '/Volumes/Seq_data/genes_E_RS.gtf'
            hts_out = os.path.join(root,cname+'_htseqcount.txt')
            htseq_count_call = 'htseq-count -f bam '+picard_fixmate_out+' '+gf+' > '+hts_out

            #run picard CollectRnaSeqMetrics (http://broadinstitute.github.io/picard/command-line-overview.html) and generate matrix of 3' to 5' bias (norm_read_dict)
            picard_rnaseqmetric_out = os.path.join(root,sort_out.strip('sorted.bam')+'RNA_metric.txt')
            picard_rnaseqchart_out = os.path.join(root,sort_out.strip('sorted.bam')+'RNA_metric.pdf')
            picard_seqmetric_call = 'java -Xmx3g -jar /netapp/home/idriver/broadinstitute-picard-fd8e773/dist/picard.jar CollectRnaSeqMetrics REF_FLAT=/Volumes/Seq_data/refFlat_mm10ERS.txt.gz STRAND_SPECIFICITY=NONE MINIMUM_LENGTH=70 CHART_OUTPUT='+picard_rnaseqchart_out+' INPUT='+picard_fixmate_out+' OUTPUT='+picard_rnaseqmetric_out

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
count_df.to_csv(os.path.join(path,'combined_spc_count_table.txt'), sep = '\t')
with open(os.path.join(path,'htseq_count_combined_spc.p'), 'wb') as fp1:
    pickle.dump(count_df, fp1)
pic_stats_df = pd.DataFrame(picard_stats_dict, index = index1)
pic_stats_df.to_csv(os.path.join(path,'combined_spc_picard_stats.txt'), sep = '\t')
norm_read_df = pd.DataFrame(norm_read_dict, index = index3)
norm_read_df.to_csv(os.path.join(path,'combined_spc_read_bias.txt'), sep = '\t')
pd.DataFrame.plot(norm_read_df)

            contents = """\
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o /netapp/home/idriver/%(result_file_name)s
#$ -e /netapp/home/idriver/error_spc
#$ -cwd
#$ -r y
#$ -j y
#$ -l netapp=10G,scratch=40G,mem_total=22G
#$ -pe smp 8
#$ -R yes
#$ -l h_rt=3:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.3
PATH=$PATH:/netapp/home/idriver/bin/samtools-0.1.19_2
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.0.13.Linux_x86_64
PATH=$PATH:/usr/bin/gunzip
export PATH
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
mkdir %(name)s
mkdir -p /netapp/home/idriver/%(result_file_name)s/%(name)s
%(tophat_cmd)s
%(cufflinks_cmd)s
%(cuffquant_cmd)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(name)s/* /netapp/home/idriver/%(result_file_name)s/%(name)s
rm -r %(name)s
date
""" % vars()
            if result_file_name == 'results_Lane8_data':
                filename = '%s.sh' % name
                write_file(filename, contents)
                print tophat_cmd
                print cufflinks_cmd
                print cuffquant_cmd
                jobid = qsub_submit(filename, fname = '%s' % name)
                print "Submitted. jobid = %d" % jobid
                # Write jobid to a file.
                import subprocess
                process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
                out, err = process.communicate()
                print(out)

path = '/Volumes/Seq_data'
count_df = pd.DataFrame(count_dict, index = gene_list)
count_df.to_csv(os.path.join(path,'pdgfra2_count_table.txt'), sep = '\t')
