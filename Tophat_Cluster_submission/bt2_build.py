import os
from subprocess import call
bt2_str=''
p= '/netapp/home/idriver/all_fa_ERCC'
for root, dirnames, filenames in os.walk(p):
  bt2_str= ' '.join(filenames)
bt2_str1 = bt2_str.strip(',')
print bt2_str1

call('cat '+bt2_str1+' > mm10_E.fa', shell=True)
