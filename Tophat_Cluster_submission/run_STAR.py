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
path = '/netapp/home/idriver/chapmanh'
out= '${TMPDIR}'
annotation_file = '/netapp/home/idriver/hg19_ERCC_bt2/Annotation/hg19_ERCC.gtf'
index_gen_loc = '/netapp/home/idriver/hg19_ERCC_bt2/hg19_ERCC/hg19_ERCC'

result_name='hg19_ERCC_STAR'
name = 'hg19_star'
name2 = 'pass2'
name3 = 'hg19_2pass'
genomeDir = '/netapp/home/idriver/hg19_ERCC_STAR/hg19_star'
star_cmd1 ='STAR --runMode genomeGenerate --genomeDir '+name+' --genomeFastaFiles /netapp/home/idriver/hg19_ERCC_bt2/hg19_ERCC/hg19_ERCC.fa  --runThreadN 24 --genomeChrBinNbits 12 --genomeSAsparseD 2'
star_cmd2 = 'STAR --genomeDir '+genomeDir+' --readFilesIn /netapp/home/idriver/Sample_hu-IPF-HTII-280-66/hu-IPF-HTII-280-66_CGAGGCTG-TATCCTCT_L008_R1_001.fastq,/netapp/home/idriver/Sample_hu-IPF-HTII-280-66/hu-IPF-HTII-280-66_CGAGGCTG-TATCCTCT_L008_R1_002.fastq /netapp/home/idriver/Sample_hu-IPF-HTII-280-66/hu-IPF-HTII-280-66_CGAGGCTG-TATCCTCT_L008_R2_001.fastq,/netapp/home/idriver/Sample_hu-IPF-HTII-280-66/hu-IPF-HTII-280-66_CGAGGCTG-TATCCTCT_L008_R2_002.fastq --runThreadN 16'
star_cmd3 = 'STAR --runMode genomeGenerate --genomeDir '+genomeDir+' --genomeFastaFiles /netapp/home/idriver/hg19_ERCC_bt2/hg19_ERCC/hg19_ERCC.fa  --sjdbFileChrStartEnd /netapp/home/idriver/hg19_ERCC_STAR/hg19_pass1 --sjdbOverhang 75 --runThreadN 16 --genomeChrBinNbits 12 --genomeSAsparseD 2'
contents = """\
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o /netapp/home/idriver/%(result_name)s
#$ -e /netapp/home/idriver/error_spc
#$ -cwd
#$ -r y
#$ -j y
#$ -l netapp=30G,scratch=100G,mem_total=40G
#$ -pe smp 16
#$ -R yes
#$ -l h_rt=8:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.3
PATH=$PATH:/netapp/home/idriver/bin/samtools-0.1.19_2
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.0.13.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/STAR-STAR_2.4.2a/source
alias STAR="~/bin/STAR-STAR_2.4.2a/source/STAR"
PATH=$PATH:/usr/bin/gunzip
export PATH
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
mkdir %(name2)s
cd %(name2)s
mkdir -p /netapp/home/idriver/%(result_name)s/%(name2)s
%(star_cmd2)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(name2)s/* /netapp/home/idriver/%(result_name)s/%(name2)s
rm -r %(name2)s
date
""" % vars()
if name2 == 'pass2':
    filename = '%s.sh' % name2
    write_file(filename, contents)
    print star_cmd2
    jobid = qsub_submit(filename, fname = '%s' % name2)
    print "Submitted. jobid = %d" % jobid
    # Write jobid to a file.
    import subprocess
    process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
    out, err = process.communicate()
    print(out)
