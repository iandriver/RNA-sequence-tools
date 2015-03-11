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
            gf = '/Volumes/Seq_data/genes_E_RS.gtf'
            hts_out = os.path.join(root,cname+'_htseqcount.txt')
            htseq_count_call = 'htseq-count -f bam '+sort_out+'.bam'+' '+gf+' > '+hts_out
            print htseq_count_call
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
        keys = sorted(samp_dict.keys(), reverse=True)
        with open(os.path.join('/netapp/home/idriver', 'sample_sheet_'+sampl_sheet_name), "wb") as outfile:
            writer = csv.writer(outfile, delimiter = "\t")
            writer.writerow(keys)
            writer.writerows(zip(*[samp_dict[key] for key in keys]))

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
