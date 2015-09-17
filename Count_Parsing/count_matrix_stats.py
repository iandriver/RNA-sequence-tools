import fnmatch
import os
import pandas as pd
import cPickle as pickle
import csv
from collections import OrderedDict




#list of file paths with mapped hits
pats = ['/Volumes/Seq_data/results_01272015', '/Volumes/Seq_data/results_spc2_n2']
#output path
path = '/Volumes/Seq_data'
#base name for final output count matrix and picard metrics
base_name = 'combined_spc'

#initialize dictonaries for collected output
fpkm_matrix_dict_g = OderedDict()
count_dict = OderedDict()
norm_read_dict = OderedDict()
picard_stats_dict = OderedDict()

#collect gene_list once since it the same between all samples
st = 1
gene_list = []
g_counts = []

for p in pats:
    for root, dirnames, filenames in os.walk(os.path.join(path,p)):
        for filename in fnmatch.filter(filenames, 'accepted_hits.bam'):
            #sorted file path
            cname = root.split('/')[-1]
            sort_out = os.path.join(out, cname, cname+'_sorted')

            #fixmate file path
            picard_fixmate_out = sort_out.strip('.bam')+'_FM.bam'

            #format htseq-count command to generate raw counts from sorted accepted hits
            hts_out = os.path.join(out,cname+'_htseqcount.txt')

            #run picard CollectRnaSeqMetrics (http://broadinstitute.github.io/picard/command-line-overview.html) and generate matrix of 3' to 5' bias (norm_read_dict)
            picard_rnaseqmetric_out = sort_out.strip('sorted.bam')+'RNA_metric.txt'
            picard_rnaseqchart_out = sort_out.strip('sorted.bam')+'RNA_metric.pdf'

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
count_df.to_csv(os.path.join(path,base_name+'_count_table.txt'), sep = '\t')
with open(os.path.join(path,'htseq_count_'base_name'.p'), 'wb') as fp1:
pickle.dump(count_df, fp1)
pic_stats_df = pd.DataFrame(picard_stats_dict, index = index1)
pic_stats_df.to_csv(os.path.join(path,base_name+'_picard_stats.txt'), sep = '\t')
norm_read_df = pd.DataFrame(norm_read_dict, index = index3)
norm_read_df.to_csv(os.path.join(path,base_name'_read_bias.txt'), sep = '\t')
pd.DataFrame.plot(norm_read_df)
