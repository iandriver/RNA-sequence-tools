import fnmatch
import os
import csv
import subprocess

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
    if fname: command += ' -N %s' % fname
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


trim_file_name = 'krasnow_sra_18_5'

path = os.path.join('/netapp/home/idriver', trim_file_name)

cuffdiff_files_1 = ''
cuffdiff_files_rest = ''
file_out= '${TMPDIR}'
trim_fasta_path = '/netapp/home/idriver/Trimmomatic-0.36/adapters/NexteraPE-PE.fa'
for root, dirnames, filenames in os.walk(path):
    if dirnames == []:
        for files in filenames:
            num = files.strip('.fastq.gz').split('_')[-1]
            if num == '1':
                pair_1 = os.path.join(root,files)
            if num == '2':
                pair_2 = os.path.join(root,files)
        srr_num = os.path.basename(root)
        file_name_end_list = ['_forward_paired.fq.gz','_forward_unpaired.fq.gz','_reverse_paired.fq.gz','_reverse_unpaired.fq.gz']
        output_text = ' '
        for extension in file_name_end_list:
            output_text = output_text+file_out+'/'+srr_num+extension+' '
        trim_cmd = 'java -jar /netapp/home/idriver/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 '+pair_1+' '+pair_2+output_text+'ILLUMINACLIP:${TMPDIR}/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

        command = """\
#!/bin/sh
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o %(path)s
#$ -e /netapp/home/idriver/error_spc
#$ -cwd
#$ -r y
#$ -j y
#$ -l netapp=10G,scratch=10G,mem_total=12G
#$ -pe smp 6
#$ -R yes
#$ -l h_rt=1:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/bin/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.8
PATH=$PATH:/netapp/home/idriver/bin/samtools-1.3
PATH=$PATH:/netapp/home/idriver/bin/tophat-2.1.0.Linux_x86_64
PATH=$PATH:/usr/bin/gunzip
export PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
cp %(trim_fasta_path)s .
%(trim_cmd)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(srr_num)s* %(path)s/%(srr_num)s/
""" % vars()
        if srr_num != '':
            print(trim_cmd)
            filename = 'trim_'+srr_num+'.sh'
            write_file(filename, command)
            jobid = qsub_submit(filename, fname=srr_num)
            print "Submitted. jobid = %d" % jobid
            # Write jobid to a file.
            import subprocess
            process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
            out, err = process.communicate()
            print(out)
