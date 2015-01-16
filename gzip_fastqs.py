import os
from subprocess import call


pathlist = []
for root, dirs, files in os.walk('/netapp/home/idriver/10242014_macs'):
  if 'fastq' in root:
    pathlist.append([root,files])
for p in pathlist:
  for f in p[1]:
    if 'gz' in f:
      gzip_call = 'gunzip -f '+p[0]+'/'+f
      print gzip_call
      call(gzip_call, shell=True)
