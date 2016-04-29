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
path = '/netapp/home/idriver/160310_D00328_0472_AC9CR6ANXX_Project_chapmanh-YX'
out= '${TMPDIR}'
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
    hisat_ref = '/netapp/home/idriver/mm10_ERCC_RS_bt2/mm10_ERCC_RS/mm10_ERS_hisat'
    splice_file = '/netapp/home/idriver/mm10_splicesites.txt'

#this next section parses the file names so that the paired end reads are in order and determines the name of the output file
#use test_qsub.py to test and modify this section locally to work for your file names
pathlist = []
run_1 = True
result_file_name = 'hisat_chapmanh-YX'
call('mkdir -p /netapp/home/idriver/%s' % result_file_name, shell=True)
for root, dirs, files in os.walk(path):
    if dirs == ['fastqc']:
        n = root.strip('/').split('/')
        out= '${TMPDIR}'
        name1 = n[-1].split('_')
        name = name1[1]
        data_file = root
        result_file = os.path.join(out,name)
        rank_tups =[]
        for f in files:
            if 'fastq' in f and ".txt" not in f:
                f_split = f.split('_')
                r_name = int(f_split[3][1])
                en_split = f_split[4].split('.')
                p_num = int(en_split[0])
                rank_tups.append((os.path.join(root,f), r_name,p_num))
        files_sorted = sorted(rank_tups, key=lambda element: (element[1], element[2]))
        files1 = ''
        files2 = ''
        for fs in files_sorted:
            if fs[1] == 1:
                files1 += fs[0]+','
            if fs[1] == 2:
                files2 += fs[0]+','
        final_files = '-1 '+ files1.strip(',')+' -2 '+files2.strip(',')
        hisat_cmd = 'hisat2 -p 14 --known-splicesite-infile '+splice_file+' -x '+hisat_ref+' --trim3 --dta-cufflinks '+final_files+' -S '+name
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
#$ -l netapp=25G,scratch=120G,mem_total=200G
#$ -pe smp 14
#$ -R yes
#$ -l h_rt=36:59:00
set echo on
date
hostname
pwd
export PATH=$PATH:${HOME}/bin
PATH=$PATH:/netapp/home/idriver/bin/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.8
PATH=$PATH:/netapp/home/idriver/bin/samtools-1.2
PATH=$PATH:/netapp/home/idriver/hisat2-2.0.3-beta
PATH=$PATH:/usr/bin/gunzip
export PATH
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
mkdir %(name)s
cd %(name)s
mkdir -p /netapp/home/idriver/%(result_file_name)s/%(name)s
%(hisat_cmd)s

# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(name)s/* /netapp/home/idriver/%(result_file_name)s/%(name)s
rm -r %(name)s
date
""" % vars()
        if name == 'K5-1':
            filename = '%s.sh' % name
            write_file(filename, contents)
            print hisat_cmd
            jobid = qsub_submit(filename, fname = '%s' % name)
            print "Submitted. jobid = %d" % jobid
            # Write jobid to a file.
            import subprocess
            process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
            out, err = process.communicate()
            print(out)
            run_1 = False
