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
path = '/netapp/home/idriver/MD_SPC_RBPJ_PNX_fastq'
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
    hisat_ref = '/netapp/home/idriver/hg19_ERCC_hisat2/hg19_ERCC_hisat2'
    refflat = 'refFlat.txt.gz'
    transcript_index = '/netapp/home/idriver/transcriptome_data_hg19_ERCC_2'
    splice_file = '/netapp/home/idriver/hg19_ERCC_splicesites.txt'
if genome == 'mouse':
    annotation_file = '/netapp/home/idriver/mm10_E_RS/genes_E_RS.gtf'
    index_gen_loc = '/netapp/home/idriver/mm10_ERCC_RS_bt2/mm10_ERCC_RS/mm10_ERCC_RS'
    refflat = 'refFlat_mm10ERS.txt.gz'
    transcript_index = '/netapp/home/idriver/transcriptome_data_mm10_RS/known_e_RS'
    hisat_ref = '/netapp/home/idriver/mm10_ERCC_hisat2/mm10_ERCC_hisat2'
    splice_file = '/netapp/home/idriver/mm10_E_RS/mm10_splicesites.txt'

#this next section parses the file names so that the paired end reads are in order and determines the name of the output file
#use test_qsub.py to test and modify this section locally to work for your file names
pathlist = []
run_1 = True
result_file_name = 'results_hisat_md_spc_pnx_rbpj_bulk'
call('mkdir -p /netapp/home/idriver/%s' % result_file_name, shell=True)
for root, dirs, files in os.walk(path):
    if dirs == [] or dirs ==['fastqc']:
        out= '${TMPDIR}'
        name = os.path.basename(root)
        data_file = root
        result_file = os.path.join(out,name)
        final_files = ''
        for f in files:
            if '_trimed_fastq' in f:
                final_files = os.path.join(root,f)
                single_end = True
            elif 'forward_paired' in f:
                single_end = False
                final_files+=os.path.join(root,f)+' '
            elif 'reverse_paired' in f:
                final_files+=os.path.join(root,f)
        #if pairs exist format paired end read command
        if not single_end:
            hisat_cmd = 'hisat2 -p 10 --phred33 --rna-strandness FR --known-splicesite-infile '+splice_file+' -x '+hisat_ref+' --dta-cufflinks '+final_files+' -S '+name+'.sam --met-file '+result_file+'/'+name+'_metrics.txt'
        #otherwise format for single end reads
        else:
            hisat_cmd = 'hisat2 -p 10 --phred33 --rna-strandness F --known-splicesite-infile '+splice_file+' -x '+hisat_ref+' --dta-cufflinks -U '+final_files+' -S '+name+'.sam --met-file '+result_file+'/'+name+'_metrics.txt'

        samtools_cmd = '/netapp/home/idriver/bin/samtools-1.3.1/samtools view -u '+name+'.sam | /netapp/home/idriver/bin/samtools-1.3.1/samtools sort -o '+name+'.sorted.bam'
        stringtie_cmd = 'stringtie -eB '+result_file+'/'+name+'.sorted.bam -o '+result_file+'/'+name+'.gff -p 10 -G '+annotation_file
        cufflinks_cmd = 'cufflinks -p 10 --max-bundle-frags 10000000 -G '+annotation_file+' -o '+result_file+' '+name+'.sorted.bam'
        cuffquant_cmd = 'cuffquant -p 10 --max-bundle-frags 10000000 -o '+result_file+' '+annotation_file+' '+name+'.sorted.bam'
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
#$ -l netapp=10G,scratch=20G,mem_total=25G
#$ -pe smp 10
#$ -R yes
#$ -l h_rt=8:59:00
set echo on
date
hostname
pwd
PATH=$PATH:/netapp/home/idriver/bin/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.8
PATH=$PATH:/netapp/home/idriver/bin/samtools-1.3.1
PATH=$PATH:/netapp/home/idriver/hisat2-2.0.4
PATH=$PATH:/netapp/home/idriver/bin/stringtie-1.2.2.Linux_x86_64
PATH=$PATH:/usr/bin/gunzip
alias samtools='/netapp/home/idriver/bin/samtools-1.3.1/samtools'
PATH=$PATH:${HOME}/bin
export PATH
echo $PATH
export TMPDIR=/scratch
echo $TMPDIR
cd $TMPDIR
mkdir %(name)s
cd %(name)s
mkdir -p /netapp/home/idriver/%(result_file_name)s/%(name)s
%(hisat_cmd)s
%(samtools_cmd)s
%(stringtie_cmd)s
%(cufflinks_cmd)s
%(cuffquant_cmd)s
# Copy the results back to the project directory:
cd $TMPDIR
cp -r %(name)s/* /netapp/home/idriver/%(result_file_name)s/%(name)s
rm -r %(name)s
date
""" % vars()
        if name != 'Sample_C1':
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
