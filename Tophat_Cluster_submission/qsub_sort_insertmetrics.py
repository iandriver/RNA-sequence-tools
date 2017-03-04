import fnmatch
import os
import pandas as pd
import cPickle as pickle
import subprocess
import csv
from collections import OrderedDict

def write_file(filename, contents):
    """Write the given contents to a text file.

    ARGUMENTS
       filename (string) - name of the file to write to, creating if it doesn't exist
       contents (string) - contents of the file to be written
    """

    # Open the file for writing
    file = open(filename, 'w')

    # Write the file contents
    file.write(contents)

    # Close the file
    file.close()

    return

def qsub_submit(command_filename, hold_jobid = None, fname = None):
    """Submit the given command filename to the queue.

    ARGUMENTS
      command_filename (string) - the name of the command file to submit

    OPTIONAL ARGUMENTS
      hold_jobid (int) - job id to hold on as a prerequisite for execution

    RETURNS
      jobid (integer) - the jobid
    """

    # Form command
    command = 'qsub'
    if name: command += ' -N %s' % fname
    if hold_jobid: command += ' -hold_jid %d' % hold_jobid
    command += ' %s' % command_filename

    # Submit the job and capture output.
    import subprocess
    print "> " + command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    print(out)

    # Match job id
    jobid = out.split(' ')[2]

    return int(jobid)

pats = ['/netapp/home/idriver/rsem_m38_Cfms-d7-Sham']
pt = pats[0]
#output path
out= '${TMPDIR}'
result_file_name = os.path.basename(pats[0])+'_insert_metrics'
#initialize dictonaries for collected output
fpkm_matrix_dict_g ={}
count_dict = {}
norm_read_dict = {}
picard_stats_dict = {}

#collect gene_list once since it the same between all samples
st = 1
gene_list = []
command_list = []
#loop through all files and sort, fix, count, collect metrics on each
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, '*.bam'):
            #sort accepted_hits.bam using samtools
            cname = root.split('/')[-1]
            sort_out = os.path.join(out, cname ,filename.strip(".bam")+"_sorted.bam")
            sort_call = "/netapp/home/idriver/bin/jre1.8.0_121/bin/java -Xmx3g -Djava.io.tmpdir="+os.path.join(out,cname)+" -jar /netapp/home/idriver/picard.jar SortSam I="+os.path.join(root,filename)+" O="+sort_out+" SORT_ORDER=coordinate TMP_DIR="+os.path.join(out,cname)
            print(sort_call)

            if not os.path.isfile(os.path.join(root,"insert_size_metrics.txt")):
                picard_sizemetric_call = '/netapp/home/idriver/bin/jre1.8.0_121/bin/java -Xmx3g -jar /netapp/home/idriver/picard.jar CollectInsertSizeMetrics I='+sort_out+' O='+os.path.join(out,cname,"insert_size_metrics.txt")+" H="+os.path.join(out,cname,'insert_size_histogram.pdf')
                print(picard_sizemetric_call)
            else:
                print(sort_out+' already sorted.')
            command_list.append([cname, sort_call, picard_sizemetric_call])


subprocess.call('mkdir -p /netapp/home/idriver/%s' % result_file_name, shell=True)


for i, calls in enumerate(command_list):
    name = calls[0]
    sorted_name = os.path.join(out,name,name+".transcript_sorted.bam")
    calls_zero = calls[1]
    calls_one = calls[2]
    contents = """\
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o /netapp/home/idriver/%(result_file_name)s
#$ -e /netapp/home/idriver/error_spc
#$ -cwd
#$ -r y
#$ -j y
#$ -l netapp=20G,scratch=25G,mem_total=5G
#$ -pe smp 2
#$ -R yes
#$ -l h_rt=10:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/bin/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.8
PATH=$PATH:/netapp/home/idriver/bin/samtools-1.3.1
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.1.0.Linux_x86_64
PATH=$PATH:/usr/bin/gunzip
export PATH
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
mkdir %(name)s
mkdir -p /netapp/home/idriver/%(result_file_name)s/%(name)s
cd %(name)s
%(calls_zero)s
%(calls_one)s
# Copy the results back to the project directory:
cd $TMPDIR
cp %(name)s/*.txt /netapp/home/idriver/%(result_file_name)s/%(name)s
cp %(name)s/*.pdf /netapp/home/idriver/%(result_file_name)s/%(name)s
rm %(sorted_name)s
rm -r %(pt)s/%(name)s
date
""" % vars()
    if name != 'A10_C52':
        filename = '%s.sh' % name
        write_file(filename, contents)
        print calls_zero
        print calls_one
        jobid = qsub_submit(filename, fname = '%s' % name)
        print "Submitted. jobid = %d" % jobid
        # Write jobid to a file.
        import subprocess
        process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
        out, err = process.communicate()
        print(out)
