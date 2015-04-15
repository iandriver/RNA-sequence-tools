#/usr/bin/env python
import commands
import os
from subprocess import call
import fnmatch

do_cells =['199', '200', '201', '202', '204', '184', '183', '132', '108', '107', '71', '69', '67']

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

path = '/netapp/home/idriver/raw_data'
out= '${TMPDIR}'
annotation_file = '/netapp/home/idriver/genes_E_RS.gtf'
index_gen_loc = '/netapp/home/idriver/mm10_ERCC_RS_bt2/mm10_ERCC_RS/mm10_ERCC_RS'


pathlist = []
result_name_list = []
for root, dirs, files in os.walk(path):
    if root.split('/')[-1] == 'raw_data':
        for lane in dirs:
            result_name = 'results_'+lane
            if ' ' in result_name:
                result_s = result_name.split(' ')
                result_name = '_'.join(result_s)
            call('mkdir -p /netapp/home/idriver/%s' % result_name, shell=True)
            result_name_list.append(result_name)
    elif root.split('/')[-2]=='raw_data' and dirs == []:
        for f in files:
            n = root.split('/')[-1]
            if ' ' in n:
                sample_name = n.split(' ')[0]
            elif '_' in n:
                sample_name = n.split('_')[0]
            for rn in result_name_list:
                if sample_name in rn:
                    result_file_name = rn
            if 'fastq' in f and ".txt" not in f:
                f_split = f.split('.')
                f_name = sample_name+'_'+f_split[0]
                root_s = root.split('/')
                if ' ' in root_s[-1]:
                    new_file = '\''+root_s[-1]+'\''
                    new_root = '/'.join(root_s[0:-1])+'/'+new_file
                    final_files = os.path.join(new_root,f)
                else:
                    final_files = os.path.join(root,f)
                print final_files
                result_file = os.path.join(out,f_name)
                tophat_cmd = 'tophat2 -p 8 --read-realign-edit-dist 0 -G '+annotation_file+' --transcriptome-index=/netapp/home/idriver/transcriptome_data_mm10_RS/known_e_RS -o '+result_file+' '+index_gen_loc+' '+final_files
                samtools_cmd = 'samtools sort '+result_file+'/'+'accepted_hits.bam accepted_hits_sorted'
                cufflinks_cmd = 'cufflinks -p 8 --max-bundle-frags 10000000 -G '+annotation_file+' -o '+result_file+' '+result_file+'/'+'accepted_hits.bam'
                cuffquant_cmd = 'cuffquant -p 8 --max-bundle-frags 10000000 -o '+result_file+' '+annotation_file+' '+result_file+'/'+'accepted_hits.bam'
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
#$ -l netapp=10G,scratch=30G,mem_total=15G
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
mkdir %(f_name)s
mkdir -p /netapp/home/idriver/%(result_file_name)s/%(f_name)s
%(tophat_cmd)s
%(cufflinks_cmd)s
%(cuffquant_cmd)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(f_name)s/* /netapp/home/idriver/%(result_file_name)s/%(f_name)s
rm -r %(f_name)s
date
""" % vars()

            if 'kidney' in f_name and f_name.split('_')[-1] in do_cells:
                call('rm -r /netapp/home/idriver/results_kidney_raw/%s*' % f_name, shell=True)
                filename = '%s.sh' % f_name
                write_file(filename, contents)
                print tophat_cmd
                print cufflinks_cmd
                print cuffquant_cmd
                jobid = qsub_submit(filename, fname = '%s' % f_name)
                print "Submitted. jobid = %d" % jobid
                # Write jobid to a file.
                import subprocess
                process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
                out, err = process.communicate()
                print(out)
