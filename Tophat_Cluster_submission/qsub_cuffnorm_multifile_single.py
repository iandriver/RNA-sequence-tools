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
sampl_sheet_name = 'sheppard_all_RS'
result_file_names = ['results_lung_raw', 'results_liver_raw', 'results_kidney_raw']

samp_dict = {}
samp_dict['sample_name'] =[]
samp_dict['group'] = []
pats = []
for fname in result_file_names:
    p = os.path.join('/netapp/home/idriver', fname)
    file_name = fname.split('_')[1]
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, '*.cxb'):
            g_cell_name = root.split('/')[-1]
            num = g_cell_name.split('_')[-1]
            if len(num) == 2:
              cell_name = '_'.join(g_cell_name.split('_')[0:-1])+'_'+'0'+num
            elif len(num) == 1:
              cell_name = '_'.join(g_cell_name.split('_')[0:-1])+'_'+'00'+num
            else:
              cell_name = g_cell_name
            samp_path = os.path.join(root,filename)
            samp_dict['sample_name'].append(samp_path)
            samp_dict['group'].append(cell_name)
    keys = sorted(samp_dict.keys(), reverse=True)
with open(os.path.join('/netapp/home/idriver', 'sample_sheet_'+sampl_sheet_name), "wb") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(keys)
    writer.writerows(zip(*[samp_dict[key] for key in keys]))
annotation_file = '/netapp/home/idriver/genes_E_RS.gtf'
index_gen_loc = '/netapp/home/idriver/mm10_ERCC_RS_bt2/mm10_ERCC_RS/mm10_ERCC_RS'
cuff_name = 'cuffnorm_'+sampl_sheet_name
cuffnorm_cmd = 'cuffnorm --use-sample-sheet -p 8 -o '+cuff_name+' '+annotation_file+' '+os.path.join('/netapp/home/idriver', 'sample_sheet_'+sampl_sheet_name)
mk_dir = 'mkdir -p '+os.path.join('/netapp/home/idriver', cuff_name)
subprocess.call(mk_dir, shell=True)

command = """\
#!/bin/sh
#!/bin/sh
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -o /netapp/home/idriver/%(cuff_name)s
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
mkdir %(cuff_name)s

%(cuffnorm_cmd)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(cuff_name)s/* /netapp/home/idriver/%(cuff_name)s
""" % vars()
filename = 'cuffnorm_'+sampl_sheet_name+'.sh'
write_file(filename, command)
jobid = qsub_submit(filename, fname=cuff_name)
print "Submitted. jobid = %d" % jobid
# Write jobid to a file.
import subprocess
process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
out, err = process.communicate()
print(out)
