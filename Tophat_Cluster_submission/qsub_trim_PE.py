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


trim_file_name = 'pengt_2nd'

path = os.path.join('/netapp/home/idriver', trim_file_name)


file_out= '${TMPDIR}'

for root, dirnames, filenames in os.walk(path):
    if dirnames == []:
        fastq_files = [files for files in filenames if '.fastq.gz' in files]

        name = os.path.basename(root)

        #subprocess.call('mkdir -p '+path+'/%s' % name, shell=True)
        if len(fastq_files) > 2:
            r1_files = [r1 for r1 in fastq_files if r1.strip('.fastq.gz').split('_')[-2] == 'R1']
            r1_files.sort()
            r2_files = [r2 for r2 in fastq_files if r2.strip('.fastq.gz').split('_')[-2] == 'R2']
            r2_files.sort()
            r1_cat_cmd = 'cat '
            r2_cat_cmd = 'cat '
            for r1,r2 in zip(r1_files,r2_files):
                r1_cat_cmd+= os.path.join(root,r1)+' '
                r2_cat_cmd+= os.path.join(root,r2)+' '
            r1_cat_cmd+= '> '+os.path.join(root,r1_files[0].strip('.fastq.gz')+'_new'+'.fastq.gz')
            r2_cat_cmd+= '> '+os.path.join(root,r2_files[0].strip('.fastq.gz')+'_new'+'.fastq.gz')
            print(r1_cat_cmd)
            if os.path.isfile(os.path.join(root,r1_files[0].strip('.fastq.gz')+'_new'+'.fastq.gz')):
                print(os.path.join(root,r1_files[0].strip('.fastq.gz')+'_new'+'.fastq.gz')+' exists.')
                pass
            else:
                process = subprocess.Popen(r1_cat_cmd, stdout=subprocess.PIPE, shell=True)
                out1, err1 = process.communicate()
                print(out1)
                print(r2_cat_cmd)
                process = subprocess.Popen(r2_cat_cmd, stdout=subprocess.PIPE, shell=True)
                out2, err2 = process.communicate()
                print(out2)
            pair_1 = os.path.join(root,r1_files[0].strip('.fastq.gz')+'_new'+'.fastq.gz')
            pair_2 = os.path.join(root,r2_files[0].strip('.fastq.gz')+'_new'+'.fastq.gz')
        if len(fastq_files) == 2:
            for files in fastq_files:
                num = files.strip('.fastq.gz').split('_')[-2]
                if num == 'R1':
                    pair_1 = os.path.join(root,files)
                if num == 'R2':
                    pair_2 = os.path.join(root,files)
        if len(fastq_files) == 1:
            trim_fasta_path = '/netapp/home/idriver/Trimmomatic-0.36/adapters/TruSeq3-SE.fa'
            output_text = ' '+file_out+'/'+name+'_trimed_fastq.gz '
            input_file = os.path.join(root,fastq_files[0])
            trim_cmd = 'java -jar /netapp/home/idriver/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 '+input_file+output_text+'ILLUMINACLIP:${TMPDIR}/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75'
        else:
            trim_fasta_path = '/netapp/home/idriver/Trimmomatic-0.36/adapters/NexteraPE-PE.fa'
            file_name_end_list = ['_forward_paired.fq.gz','_forward_unpaired.fq.gz','_reverse_paired.fq.gz','_reverse_unpaired.fq.gz']
            output_text = ' '
            for extension in file_name_end_list:
                output_text = output_text+file_out+'/'+name+extension+' '
            trim_cmd = 'java -jar /netapp/home/idriver/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 '+pair_1+' '+pair_2+output_text+'ILLUMINACLIP:${TMPDIR}/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70'

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
#$ -l netapp=5G,scratch=8G,mem_total=10G
#$ -pe smp 4
#$ -R yes
#$ -l h_rt=1:00:00
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
cp -r %(name)s_* %(path)s/%(name)s/
""" % vars()
        if name != '': #['BU3_C36_F9','BU3_C76_E12','ips17_C13_C3','ips17_C26_E2','ips17_C27_E3','ips17_C25_E1']:
            print(trim_cmd)
            filename = 'trim_'+name+'.sh'
            write_file(filename, command)
            jobid = qsub_submit(filename, fname=name)
            print "Submitted. jobid = %d" % jobid
            # Write jobid to a file.
            import subprocess
            process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
            out, err = process.communicate()
            print(out)
