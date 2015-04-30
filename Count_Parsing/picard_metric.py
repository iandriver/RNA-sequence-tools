import fnmatch
import os
import pandas as pd
import cPickle as pickle
import subprocess
import csv

path = '/Volumes/Seq_data/results_Lane5_data/d1_neg'
picard_out = 'd1_negRNA_metric.txt'
with open(os.path.join(path, picard_out), mode='r') as infile:
    hts_tab = csv.reader(infile, delimiter = '\t')
    for i, l in enumerate(hts_tab):
        if i == 6:
            index1 = l
        if i == 7:
            pass
        if i == 10:
            index2 = l
        if i > 10 and i <= 111:
            print l[0], l[1]
