import os
from subprocess import call

print 'y'
path = '/Users/idriver/RockLab-files/07022014'
#path to gene annotation file
annotation_file = '/Users/idriver/Downloads/Mus_musculus_Ensembl_GRCm38/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf'
#path to gene idex file
index_gen_loc = '/Users/idriver/Downloads/Mus_musculus_Ensembl_GRCm38/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome'
out ='/Users/idriver/RockLab-files/pdgfra_test'
pathlist = []
#make a list of paths, files that contain fastq sequence data
for root, dirs, files in os.walk(path):
  if 'fastq' in root:
    pathlist.append([root,files])
#extract cell name and create filename for each cell sequence
for p in pathlist:
  n = p[0].split('/')
  n1 = n[-1].split('_')
  name = n1[0]
  if name != 'C10':
    print name
    result_file = os.path.join(out, name)
    print result_file
    try:
      os.mkdir(result_file)
    except OSError:
      pass
    input_files=''
    #parse inpout file path
    r_num = []
    for f in p[1]:
      if 'fastq' in f:
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
      print sort_num,'s'
      final_files = sort_num[0]+' '+sort_num[1].strip(', ')
      print sort_num
    cell_number = int(name.strip('C'))
    tophat_cmd = 'tophat2 -p 4 -r 50 -G '+annotation_file+' -o '+result_file+' '+index_gen_loc+' '+final_files
    cufflinks_cmd = 'cufflinks -p 4 -G '+annotation_file+' -o '+result_file+' '+result_file+'/'+'accepted_hits.bam'
    cuffquant_cmd = 'cuffquant -p 4 -o '+result_file+' '+annotation_file+' '+result_file+'/'+'accepted_hits.bam'
    print tophat_cmd
    call(tophat_cmd, shell=True)
    call(cufflinks_cmd, shell=True)
    call(cuffquant_cmd, shell=True)
