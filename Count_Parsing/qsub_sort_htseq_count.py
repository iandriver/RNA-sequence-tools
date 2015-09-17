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

#list of file paths with mapped hits
pats = ['results_scler_ht280', 'results_chapmanh-hu-IPF-HTII-280', 'results_norm_alpha6']
#output path
path = '/netapp/home/idriver'
base_name = 'combined_spc'
result_file_name = 'count-picard_'+base_name
out= '${TMPDIR}'
genome = 'human'
if genome == 'human':
    annotation_file = '/netapp/home/idriver/hg19_ERCC_bt2/Annotation/hg19_ERCC.gtf'
    index_gen_loc = '/netapp/home/idriver/hg19_ERCC_bt2/hg19_ERCC/hg19_ERCC'
    refflat = 'refFlat.txt.gz'
if genome == 'mouse':
    annotation_file = '/netapp/home/idriver/genes_E_RS.gtf'
    index_gen_loc = '/netapp/home/idriver/mm10_ERCC_RS_bt2/mm10_ERCC_RS/mm10_ERCC_RS'
    refflat = 'refFlat_mm10ERS.txt.gz'

command_list = []

#loop through all files and sort, fix, count, collect metrics on each
for p in pats:
    for root, dirnames, filenames in os.walk(os.path.join(path,p)):
        for filename in fnmatch.filter(filenames, 'accepted_hits.bam'):
            #sort accepted_hits.bam using samtools
            cname = root.split('/')[-1]
            sort_out = os.path.join(out, cname, cname+'_sorted')
            sam_sort_calln = 'samtools sort -n '+os.path.join(root,filename)+' '+sort_out
            sam_sort_call = 'samtools sort '+os.path.join(root,filename)+' '+sort_out

            #run picard_fixmate to clean up paired end reads in accepted_hits.bam (sorted)
            picard_fixmate_out = sort_out.strip('.bam')+'_FM.bam'
            picard_fixmate_call = 'java -Xmx3g -jar /netapp/home/idriver/broadinstitute-picard-fd8e773/dist/picard.jar FixMateInformation INPUT='+sort_out+'.bam OUTPUT='+picard_fixmate_out+' AS=true SORT_ORDER=coordinate'

            #format htseq-count command to generate raw counts from sorted accepted hits
            hts_out = os.path.join(out,cname+'_htseqcount.txt')
            htseq_count_call = 'python -m HTSeq.scripts.count -f bam '+picard_fixmate_out+' '+annotation_file+' > '+hts_out

            #run picard CollectRnaSeqMetrics (http://broadinstitute.github.io/picard/command-line-overview.html) and generate matrix of 3' to 5' bias (norm_read_dict)
            picard_rnaseqmetric_out = sort_out.strip('sorted.bam')+'RNA_metric.txt'
            picard_rnaseqchart_out = sort_out.strip('sorted.bam')+'RNA_metric.pdf'
            picard_seqmetric_call = 'java -Xmx3g -jar /netapp/home/idriver/broadinstitute-picard-fd8e773/dist/picard.jar CollectRnaSeqMetrics REF_FLAT=/netapp/home/idriver/'+refflat+' STRAND_SPECIFICITY=NONE MINIMUM_LENGTH=70 CHART_OUTPUT='+picard_rnaseqchart_out+' INPUT='+picard_fixmate_out+' OUTPUT='+picard_rnaseqmetric_out
            command_list.append([cname, sam_sort_call, picard_fixmate_call, htseq_count_call, picard_seqmetric_call])

subprocess.call('mkdir -p /netapp/home/idriver/%s' % result_file_name, shell=True)
for i, calls in enumerate(command_list):
    name = calls[0]
    calls_zero = calls[1]
    calls_one = calls[2]
    calls_two = calls[3]
    calls_three = calls[4]
    contents = """\
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o /netapp/home/idriver/%(result_file_name)s
#$ -e /netapp/home/idriver/error_spc
#$ -cwd
#$ -r y
#$ -j y
#$ -l netapp=10G,scratch=10G,mem_total=12G
#$ -pe smp 8
#$ -R yes
#$ -l h_rt=1:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.3
PATH=$PATH:/netapp/home/idriver/bin/samtools-0.1.2
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.0.13.Linux_x86_64
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
%(calls_two)s
%(calls_three)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(name)s/* /netapp/home/idriver/%(result_file_name)s/%(name)s
rm -r %(name)s
date
""" % vars()
    if i != 0:
        filename = '%s.sh' % name+'_'+str(i)
        write_file(filename, contents)
        print calls_zero
        print calls_one
        print calls_two
        print calls_three
        jobid = qsub_submit(filename, fname = '%s' % name+'_'+str(i))
        print "Submitted. jobid = %d" % jobid
        # Write jobid to a file.
        import subprocess
        process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
        out, err = process.communicate()
        print(out)
