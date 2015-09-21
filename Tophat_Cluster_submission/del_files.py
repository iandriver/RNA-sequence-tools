import os
from subprocess import call

for path, d, files in os.walk('/netapp/home/idriver/results_scler_ht280'):
  for f in files:
    if f != 'accepted_hits.bam':
      call("rm "+os.path.join(path,f), shell=True)
