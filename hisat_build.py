#/usr/bin/env python
import commands
import os
from subprocess import call


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
def qsub_submit(command_filename, hold_jobid = None, name = None):
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
	if name: command += ' -N %s' % name
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
out= '${TMPDIR}'

base_name = 'mm10_ERCC_hisat2'
path_ss= '/netapp/home/idriver/mm10_E_RS/mm10_splicesites.txt'
path_exon= '/netapp/home/idriver/mm10_E_RS/mm10_exons.txt'
path_fa='/netapp/home/idriver/mm10_E_RS/mm10_ERCC_RS.fa'

command = 'hisat2-build -p 6 --ss '+path_ss+' --exon '+path_exon+' '+path_fa+' '+base_name
print command
call('mkdir -p /netapp/home/idriver/'+base_name, shell=True)
contents = """\
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o /netapp/home/idriver/rsem_m38
#$ -e /netapp/home/idriver/error_spc
#$ -cwd
#$ -r y
#$ -j y
#$ -l netapp=40G,scratch=100G,mem_total=150G
#$ -pe smp 6
#$ -R yes
#$ -l h_rt=7:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/bin/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.8
PATH=$PATH:/netapp/home/idriver/bin/samtools-1.2
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.1.0.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/hisat2-2.0.4
PATH=$PATH:/netapp/home/idriver/bin/STAR-STAR_2.4.2a/source
PATH=$PATH:/usr/bin/gunzip
export PATH
alias STAR="/netapp/home/idriver/bin/STAR-STAR_2.4.2a/source/STAR"
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR

%(command)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r $TMPDIR /netapp/home/idriver/%(base_name)s
rm -r $TMPDIR
date
""" % vars()

filename = 'hisat2_build.sh'
write_file(filename, contents)
jobid = qsub_submit(filename, name = 'hisat_build')
print "Submitted. jobid = %d" % jobid
# Write jobid to a file.
import subprocess
process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
out, err = process.communicate()
print(out)
