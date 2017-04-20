import subprocess
import os
import fnmatch


result_file_names = ['pengt']
path_to_file = '/netapp/home/idriver'
pats = []
new_filename1 = os.path.join(path_to_file,'pengt_1st')
new_filename2 = os.path.join(path_to_file,'pengt_2nd')
print(new_filename1)
try:
    os.mkdir(new_filename1)
    os.mkdir(new_filename2)
except OSError:
    print('File already exists.')
for fname in result_file_names:
    pats.append(os.path.join(path_to_file, fname))
cells_seen = []
for p in pats:
    for root, dirnames, filenames in os.walk(p):
        for filename in fnmatch.filter(filenames, '*.fastq.gz'):
                if '1st' in filename:
                    well_num = filename.split('_')[5]
                    current_filepath= os.path.join(new_filename1,well_num+'_1st')
                    print(current_filepath)
                    try:
                        os.mkdir(current_filepath)
                    except OSError:
                        print('Subdirectory already exists.')
                    subprocess.call('mv '+root+'/'+filename+' '+current_filepath+'/', shell=True)
                    print('mv '+root+'/'+filename+' '+current_filepath+'/')
                if '2nd' in filename:
                    well_num = filename.split('_')[5]
                    current_filepath= os.path.join(new_filename2,well_num+'_2nd')
                    print(current_filepath)
                    try:
                        os.mkdir(current_filepath)
                    except OSError:
                        print('Subdirectory already exists.')
                    subprocess.call('mv '+root+'/'+filename+' '+current_filepath+'/', shell=True)
                    print('mv '+root+'/'+filename+' '+current_filepath+'/')
