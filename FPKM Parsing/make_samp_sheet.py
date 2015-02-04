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

samp_dict = {}
samp_dict['sample_name'] =[]
samp_dict['group'] = []
pats = ['/netapp/home/idriver/results_pdgfra1_ctrl_pnx']
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, '*.cxb'):
            g_cell_name = (root.split('/')[-1])
            if len(g_cell_name.split('_')[0]) == 2:
              cell_name = 'C0'+g_cell_name.split('_')[0][-1]
            else:
              cell_name = g_cell_name.split('_')[0]
            group_name = g_cell_name.split('_')[1]+'_'+cell_name
            samp_path = os.path.join(root,filename)
            samp_dict['sample_name'].append(samp_path)
            samp_dict['group'].append(group_name)
keys = sorted(samp_dict.keys(), reverse=True)
with open("/netapp/home/idriver/pdgfra_sample_sheet.txt", "wb") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(keys)
    writer.writerows(zip(*[samp_dict[key] for key in keys]))

result_file_name = 'results_pdgfra1_ctrl_pnx'
cuffnorm_cmd = 'cuffnorm --use-sample-sheet -p 8 -o sc_expr_out /netapp/home/idriver/mm10_ERCC/genes/genes.gtf /netapp/home/idriver/pdgfra_sample_sheet.txt'
name ='sc_expr_out'
command = """\
#!/bin/sh
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
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
mkdir %(name)s
mkdir -p /netapp/home/idriver/%(result_file_name)s/%(name)s
%(cuffnorm_cmd)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(name)s/* /netapp/home/idriver/%(result_file_name)s/%(name)s
""" % vars()

filename = 'cuffnorm.sh'
write_file(filename, command)
jobid = qsub_submit(filename, fname='cuffnorm_pdgfra')
print "Submitted. jobid = %d" % jobid
# Write jobid to a file.
import subprocess
process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
out, err = process.communicate()
print(out)
