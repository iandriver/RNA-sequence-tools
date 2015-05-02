#/usr/bin/env python
import commands
import os
from subprocess import call

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

path = '/netapp/home/idriver/08182014'
out= '${TMPDIR}'
annotation_file = '/netapp/home/idriver/Mus_musculus_Ensembl_GRCm38/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf'
index_gen_loc = '/netapp/home/idriver/Mus_musculus_Ensembl_GRCm38/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome'

pathlist = []
for root, dirs, files in os.walk(path):
  if 'fastq' in root:
    pathlist.append([root,files])
for p in pathlist:
  n = p[0].strip('/').split('_')
  name = n[0].split('/')[-1]
  data_file = p[0]
  result_file = os.path.join(out,name)
  input_files=''
  r_num = []
  for f in p[1]:
    if 'fastq' in f and 'qz' not in f:
      f_split = f.split('_')
      r_name = (f_split[3][1])
      en_split = f_split[4].split('.')
      p_num = en_split[0].strip('00')
      rank = r_name+p_num
      r_num.append(int(rank))
      input_files+=os.path.join(p[0],f)+' '
  in_split = input_files.split(' ')
  sort_num = [x for (y,x) in sorted(zip(r_num,in_split))]
  if len(in_split) > 2:
    name_build = ''
    for i, mul_f in enumerate(sort_num):
      if 'fastq' in mul_f:
        if i == len(in_split)-1:
          name_build+=mul_f
        elif i < (len(in_split)/2)-1 or i > (len(in_split)/2)-1:
          name_build+= mul_f+','
        elif i ==  (len(in_split)/2)-1:
          name_build+= mul_f+' '
    final_files = name_build.strip(',')
  elif len(in_split) == 2:
    final_files = sort_num[0]+' '+sort_num[1].strip(',')
  cell_number = int(name.strip('C'))
  tophat_cmd = 'tophat2 -p 8 -r 50 --no-coverage-search --keep-tmp -G '+annotation_file+' --transcriptome-index=/netapp/home/idriver/transcriptome_data/known -o '+result_file+' '+index_gen_loc+' '+final_files
  cufflinks_cmd = 'cufflinks -p 4 -G '+annotation_file+' -o '+result_file+' '+result_file+'/'+'accepted_hits.bam'
  cuffquant_cmd = 'cuffquant -p 4 -o '+result_file+' '+annotation_file+' '+result_file+'/'+'accepted_hits.bam'
  print cell_number