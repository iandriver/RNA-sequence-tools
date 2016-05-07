import fnmatch
import os
import subprocess


path_to_file= '/Volumes/Drobo/Seq_data/krasnow_sra_18_5'
sra_list = list(range(1033917,1033935))
for sra_num in sra_list:
        sra = 'SRR'+str(sra_num)
        subprocess.call('fastq-dump -O '+os.path.join(path_to_file,sra)+' --split-files '+sra, shell=True)
        subprocess.call('gzip '+os.path.join(path_to_file,sra,sra+'_1.fastq')+' '+os.path.join(path_to_file,sra,sra+'_2.fastq'), shell=True)
