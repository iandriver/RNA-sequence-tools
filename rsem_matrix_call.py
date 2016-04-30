import fnmatch
import os
import subprocess


path_to_file= '/Volumes/Drobo/YX-chapman_rsem'

gene_list = []
iso_list = []
for root, dirnames, filenames in os.walk(path_to_file):
    for f in filenames:
        if '.genes.results' in f:
            gpat = os.path.join(root,f)
            gene_list.append(gpat)
        if '.isoforms.results' in f:
            ipat = os.path.join(root,f)
            iso_list.append(ipat)

subprocess.call('rsem-generate-data-matrix '+' '.join(gene_list)+' > '+os.path.basename(path_to_file)+'.gene.count.matrix',shell=True)
subprocess.call('rsem-generate-data-matrix '+' '.join(iso_list)+' > '+os.path.basename(path_to_file)+'.iso.count.matrix',shell=True)
