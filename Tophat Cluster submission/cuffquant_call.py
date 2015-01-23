import os
from subprocess import call



annotation_file = '/Users/idriver/Downloads/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2014-05-23-16-05-10/Genes/genes.gtf'
result_file = '/Users/idriver/RockLab-files/test'
for h, k, l in os.walk(result_file):
  g_cell_name = (h.split('/')[-1])
  if g_cell_name != 'logs' and 'cuffquant' not in g_cell_name and g_cell_name != 'tmp' and g_cell_name != 'test':
    to_call = 'cuffquant -p 4 -o '+result_file+'/'+g_cell_name+'/'+g_cell_name+'_cuffquant_out '+annotation_file+' '+result_file+'/'+g_cell_name+'/'+'accepted_hits.bam'
    print to_call

    call(to_call, shell=True)
