import fnmatch
import os
import pandas as pd

table=[]
for root, dirnames, filenames in os.walk('/Users/idriver/RockLab-files/finished'):
  for filename in fnmatch.filter(filenames, '*.cxb'):
      g_cell_name = (root.split('/')[-1])
      samp_path = os.path.join(root,filename)
      table.append([samp_path, g_cell_name])

header = ['sample_name', 'group']
df = pd.DataFrame(table, columns=header)
df.to_csv('sample_sheet.txt', sep = '\t')
