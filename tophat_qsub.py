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
  out = process.stdout.readline()
  print(out)

  # Match job id
  import re
  matches = re.match('Your job-array (\d+).(\d+)-(\d+):(\d+)', out)
  jobid = matches.group(1)

  return int(jobid)

path = '/netapp/home/idriver/06062014'
out= '/netapp/home/idriver/spc_results'
annotation_file = '/netapp/home/idriver/genes.gtf'
index_gen_loc = '/netapp/home/idriver/mm10/mm10'
pathlist = []
for root, dirs, files in os.walk(path):
  if 'fastq' in root:
    pathlist.append([root,files])
for p in pathlist:
  n = p[0].strip('./').split('_')
  name = n[0]
  result_file = os.path.join(out, name)
  try:
    os.mkdir(result_file)
  except OSError:
    pass
  input_files=''
  for f in p[1]:
    if 'gz' in f:
      f_split = f.split('_')
      r_name = f_split[3]
      en_split = f_split[4].split('.')
      r_num = en_split[0]
      input_files+=os.path.join(os.getcwd(),p[0].strip('./'),f)+' '
  if len(p[1]) > 2:
    in_split = input_files.split(' ')
    name_build = ''
    for i, mul_f in enumerate(in_split):
      if 'gz' in mul_f:
        if i == len(p[1])-1:
          name_build+=mul_f
        elif i < (len(p[1])/2)-1 or i > (len(p[1])/2)-1:
          name_build+= mul_f+','
        elif i ==  (len(p[1])/2)-1:
          name_build+= mul_f+' '
  if len(p[1]) >2:
    final_files = name_build
  else:
    final_files = input_files
  for h, k, l in os.walk(result_file):
    g_cell_name = (h.split('/')[-1])
    cell_number = int(g_cell_name.strip('C'))
    tophat_cmd = 'tophat -p 4 -r 50 -G '+annotation_file+' -o '+result_file+' '+index_gen_loc+' '+final_files
    cufflinks_cmd = 'cufflinks -p 4 -G '+annotation_file+' -o '+result_file+'/'+g_cell_name+'_cufflinks_out '+result_file+'/'+'accepted_hits.bam'
    cuffquant_cmd = 'cuffquant -p 4 -o '+result_file+'/'+g_cell_name+'/'+g_cell_name+'_cuffquant_out '+annotation_file+' '+result_file+'/'+g_cell_name+'/'+'accepted_hits.bam'
    # Write script.
    contents = """\
    #!/bin/bash
    #$ -l arch=linux-x64
    #$ -S /bin/bash
    #$ -o /netapp/home/idriver/results_spc
    #$ -e /netapp/home/idriver/error_spc
    #$ -cwd
    #$ -r y
    #$ -j y
    #$ -l netapp=2G,scratch=1G,mem_total=4G
    #$ -pe smp 4
    #$ -R yes
    #$ -l h_rt=3:59:00


    set echo on

    date
    hostname
    pwd

    export PATH=$PATH:${HOME}/bin
    PATH=$PATH:/netapp/home/idriver/bin/cufflinks-2.1.1.Linux_x86_64
    PATH=$PATH:/netapp/home/idriver/bin/bowtie2-2.2.3
    PATH=$PATH:/netapp/home/idriver/bin/samtools-0.1.19_2
    export PATH
    echo $PATH
    tophat %(tophat_cmd)d
    tophat_pid = $!
    wait $tophat_pid
    cufflinks %(cufflinks_cmd)d
    cuffquant %(cuffquant_cmd)d

    date
    """ % vars()
    filename = 'SPC-C%d.cmd' % cell_number
    write_file(filename, contents)
    jobid = qsub_submit(filename, name = 'C%d' % cell_number)
    print "Submitted. jobid = %d" % jobid
    # Write jobid to a file.
    import subprocess
    process = subprocess.Popen('echo %d > jobids' % jobid, stdout=subprocess.PIPE, shell = True)
    out, err = process.stdout.readline()
    print(out)
    break
