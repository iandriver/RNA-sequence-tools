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

call_files = ''
p= '/netapp/home/idriver/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes'
for root, dir, files in os.walk(p):
	for f in files:
		call_files += os.path.join(p, f)+','

call_file = call_files.strip(',')

command = 'bowtie2-build -f '+call_file+' hg19_ERCC'
print command

contents = """\
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o /netapp/home/idriver/hg19_ERCC_bt2
#$ -e /netapp/home/idriver/error_spc
#$ -cwd
#$ -r y
#$ -j y
#$ -l netapp=50G,scratch=100G,mem_total=60G
#$ -pe smp 8
#$ -R yes
#$ -l h_rt=5:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.3
PATH=$PATH:/netapp/home/idriver/bin/samtools-0.1.19_2
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.0.13.Linux_x86_64
export PATH
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
mkdir $TMPDIR/hg19_ERCC
cd $TMPDIR/hg19_ERCC
mkdir -p /netapp/home/idriver/hg19_ERCC_bt2
%(command)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r hg19_ERCC /netapp/home/idriver/hg19_ERCC_bt2
rm -r hg19_ERCC
date
""" % vars()

filename = 'b2_build.sh'
write_file(filename, contents)
jobid = qsub_submit(filename, name = 'b2_build')
print "Submitted. jobid = %d" % jobid
# Write jobid to a file.
import subprocess
process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
out, err = process.communicate()
print(out)	
			