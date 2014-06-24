import os
import fnmatch
import sys, csv ,operator

for root, dirnames, filenames in os.walk('/Users/idriver/RockLab-files/test'):
  for filename in fnmatch.filter(filenames, '*.fpkm_tracking'):
    if filename =='isoforms.fpkm_tracking':
      data = csv.reader(open(os.path.join(root, filename), 'rU'),delimiter='\t')
      header = next(data, None)  # returns the headers or `None` if the input is empty
      sortedlist = sorted(data, key=operator.itemgetter(0))
      #now write the sorte result into new CSV file
      with open(root+'/'+root.split('/')[-1]+'_isoforms.fpkm_tracking', "wb") as f:
          fileWriter = csv.writer(f, delimiter='\t')
          fileWriter.writerow(header)
          for row in sortedlist:
              fileWriter.writerow(row)
