import os
import pandas as pd
import cPickle as pickle
import subprocess
import csv

def make_count_matrix(dirs_in):
    count_dict = {}
    single_file_stop = True
    rows =[]
    headers = ['Gene_ID']
    for p in dirs_in:
        head_stop = True
        for root, dirnames, filenames in os.walk(p):
            cname = root.split('/')[-1]
            hts_out = os.path.join(root,cname+'_htseqcount.txt')
            for f in filenames:
                if cname+'_htseqcount.txt' == f:
                    with open(hts_out, mode='r') as infile:
                        hts_tab = csv.reader(infile, delimiter = '\t')
                        for i, l in enumerate(hts_tab):
                            if head_stop == 1:
                                headers.append(cname)
                            if single_file_stop:
                                rows.append({'Gene_ID':l[0], cname:l[1]})
                            else:
                                rows[i][cname] = l[1]
                        single_file_stop = False
        head_stop = False
    with open(os.path.join('/Volumes/Seq_data', 'spc2_counts.txt'), "wb") as outfile:
        writer = csv.DictWriter(outfile, headers, delimiter = "\t")
        writer.writeheader()
        writer.writerows(rows)

pats = ['/Volumes/Seq_data/results_Lane1_data', '/Volumes/Seq_data/results_Lane2_data', '/Volumes/Seq_data/results_Lane3_data', '/Volumes/Seq_data/results_Lane4_data']

make_count_matrix(pats)
