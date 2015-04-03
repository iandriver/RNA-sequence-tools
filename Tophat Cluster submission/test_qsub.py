#/usr/bin/env python
import commands
import os
from subprocess import call


#This is for adjusting parameters on your local system to generate proper file names for tophat_qsub
path = '/Volumes/Seq_data/02262015'
out= '${TMPDIR}'
annotation_file = '/netapp/home/idriver/mm10_ERCC/genes/genes.gtf'
index_gen_loc = '/netapp/home/idriver/mm10_ERCC/Bowtie2Index_mm10/mm10_ERCC'
result_file_name = 'results_pdgfra1_ctrl_pnx'

pathlist = []
for root, dirs, files in os.walk(path):

  if dirs == []:
      print root, files
      pathlist.append([root,files])
for p in pathlist:
  n = p[0].strip('/').split('/')
  print n, 'n'
  name = n[-1]
  print name, 'name'
  data_file = p[0]
  result_file = os.path.join(out,name)
  input_files=''
  r_num = []
  for f in p[1]:
    if 'fastq' in f and ".txt" not in f:
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
  print final_files
  tophat_cmd = 'tophat2 -p 8 -r 50 -G '+annotation_file+' --transcriptome-index=/netapp/home/idriver/transcriptome_data_mm10/known -o '+result_file+' '+index_gen_loc+' '+final_files
  cufflinks_cmd = 'cufflinks -p 8 -G '+annotation_file+' -o '+result_file+' '+result_file+'/'+'accepted_hits.bam'
  cuffquant_cmd = 'cuffquant -p 8 -o '+result_file+' '+annotation_file+' '+result_file+'/'+'accepted_hits.bam'
  # Write script.
