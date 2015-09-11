import os
from subprocess import call


pathlist = []
for root, dirs, files in os.walk('/netapp/home/idriver/Sample_hu-IPF-HTII-280-66'):
    pathlist.append([root,files])
for p in pathlist:
  for f in p[1]:
    if 'gz' in f:
      gzip_call = 'gunzip -f '+p[0]+'/'+f
      print gzip_call
      call(gzip_call, shell=True)
