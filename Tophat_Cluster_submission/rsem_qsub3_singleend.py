#/usr/bin/env python
import commands
import os
from subprocess import call
import fnmatch


finished_cells =[]

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

#paths to raw reads and annotation and index of genome
path = '/netapp/home/idriver/liver_raw'
genome = 'mouse'
if genome == 'human2':
    annotation_file = '/netapp/home/idriver/hg19_all_bt2/Annotation/hg19_all.gtf'
    index_gen_loc = '/netapp/home/idriver/hg19_all_bt2/hg19_all/hg19_all'
    refflat = '/netapp/home/idriver/hg19_all_bt2/Annotation/refFlat_hg19_all.txt.gz'
    transcript_index = '/netapp/home/idriver/transcriptome_data_hg19_all/known_hg19_all'
if genome == 'human1':
    annotation_file = '/netapp/home/idriver/hg19_ERCC_bt2/Annotation/hg19_ERCC.gtf'
    index_gen_loc = '/netapp/home/idriver/hg19_ERCC_bt2/hg19_ERCC/hg19_ERCC'
    refflat = 'refFlat.txt.gz'
    transcript_index = '/netapp/home/idriver/transcriptome_data_hg19_ERCC_2'
if genome == 'mouse':
    annotation_file = '/netapp/home/idriver/genes_E_RS.gtf'
    index_gen_loc = '/netapp/home/idriver/mm10_ERCC_RS_bt2/mm10_ERCC_RS/mm10_ERCC_RS'
    refflat = 'refFlat_mm10ERS.txt.gz'
    transcript_index = '/netapp/home/idriver/transcriptome_data_mm10_RS/known_e_RS'
rsem_ref = '/netapp/home/idriver/rsem_m38/rsem_m38/GRCm38'
#this next section parses the file names so that the paired end reads are in order and determines the name of the output file
#use test_qsub.py to test and modify this section locally to work for your file names
pathlist = []
run_1 = True
result_file_name = 'DS_liver'
call('mkdir -p /netapp/home/idriver/%s' % result_file_name, shell=True)
for root, dirs, files in os.walk(path):
    if dirs == []:
        n = root.strip('/').split('/')

        for f in files:
            if 'fastq' in f and ".txt" not in f:
                name = f.split('.')[0]
                out= os.path.join('${TMPDIR}', name)
                print name
                data_file = root
                result_file = os.path.join(out,name)
                rank_tups =[]
                final_files = os.path.join(root,f)
            command = 'rsem-calculate-expression --bowtie2 -p 8 '+final_files+' '+rsem_ref+' '+out
            # Write bash script for qsub.
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
PATH=$PATH:/netapp/home/idriver/bin/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.8
PATH=$PATH:/netapp/home/idriver/bin/samtools-1.2
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.1.0.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/RSEM-1.2.28
PATH=$PATH:/netapp/home/idriver/bin/STAR-STAR_2.4.2a/source
PATH=$PATH:/usr/bin/gunzip
alias STAR="/netapp/home/idriver/bin/STAR-STAR_2.4.2a/source/STAR"
export PATH
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
mkdir -p /netapp/home/idriver/%(result_file_name)s/%(name)s
%(command)s
# Copy the results back to the project directory:
cd $TMPDIR
cp $TMPDIR/%(name)s* /netapp/home/idriver/%(result_file_name)s/%(name)s
date
""" % vars()
            if name != '':
                filename = '%s.sh' % name
                write_file(filename, contents)
                print command
                jobid = qsub_submit(filename, fname = '%s' % name)
                print "Submitted. jobid = %d" % jobid
                # Write jobid to a file.
                import subprocess
                process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
                out, err = process.communicate()
                print(out)
                run_1 = False
