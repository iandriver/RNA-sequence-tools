import fnmatch
import os
import pandas as pd
import cPickle as pickle
import subprocess
import csv

fpkm_matrix_dict_g ={}
pats = ['/Volumes/Seq_data/results_Lane5_data', '/Volumes/Seq_data/results_Lane6_data', '/Volumes/Seq_data/results_Lane7_data', '/Volumes/Seq_data/results_Lane8_data']
count_dict = {}
st = 1
gene_list = []
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, 'accepted_hits.bam'):
            cname = root.split('/')[-1]
            sort_out = os.path.join(root,cname+'_sorted')
            sam_sort_calln = 'samtools sort -n '+os.path.join(root,filename)+' '+sort_out
            sam_sort_call = 'samtools sort '+os.path.join(root,filename)+' '+sort_out
            print sam_sort_call
            if not os.path.isfile(sort_out+'.bam'):
                process = subprocess.Popen(sam_sort_call, stdout=subprocess.PIPE, shell=True)
                out, err = process.communicate()
                print(out)
            else:
                print sort_out+'.bam already exists'
                pass
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

path = '/Volumes/Seq_data'
count_df = pd.DataFrame(count_dict, index = gene_list)
count_df.to_csv(os.path.join(path,'pdgfra2_count_table.txt'), sep = '\t')
