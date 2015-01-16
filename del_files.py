import os
from subprocess import call

for path, d, files in os.walk('/Volumes/seq_data/Lung_sabre'):
  for f in files:
    if "index" in f:
      call("rm \""+os.path.join(path,f)+"\"", shell=True)
