import fnmatch
import os
import cPickle as pickle


for root, dirnames, filenames in os.walk('/Users/idriver/RockLab-files/test'):
  for filename in filenames:
    if filename == 'accepted_hits.bam':
      g_cell_name = (root.split('/')[-1])
      os.rename(root + os.sep + filename, root + os.sep + g_cell_name+'_'+str(filename))
