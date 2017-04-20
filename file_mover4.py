import subprocess
import os
import fnmatch


result_file_names = ['170414_D00108_0693_AHJG57BCXY_rockj-Tcf21-Ctr', '170414_D00108_0693_AHJG57BCXY_rockj-Tcf21-Bleo']
path_to_file = '/netapp/home/idriver'
pats = []

for fname in result_file_names:
    pats.append(os.path.join(path_to_file, fname))
cells_seen = []
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
                if 'ctrl' in filename:
                    well_num = filename.split('_')[2]
                    current_filepath= os.path.join(root, 'Tcf21_ctrl_'+well_num)
                    print(current_filepath)
                    try:
                        os.mkdir(current_filepath)
                    except OSError:
                        print('Subdirectory already exists.')
                    subprocess.call('mv '+root+'/'+filename+' '+current_filepath+'/', shell=True)
                    print('mv '+root+'/'+filename+' '+current_filepath+'/')
                if 'Bleo' in filename:
                    well_num = filename.split('_')[2]
                    current_filepath= os.path.join(root, 'Tcf21_bleo_'+well_num)
                    print(current_filepath)
                    try:
                        os.mkdir(current_filepath)
                    except OSError:
                        print('Subdirectory already exists.')
                    subprocess.call('mv '+root+'/'+filename+' '+current_filepath+'/', shell=True)
                    print('mv '+root+'/'+filename+' '+current_filepath+'/')
