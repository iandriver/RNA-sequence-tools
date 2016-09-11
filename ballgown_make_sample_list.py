import os
import pandas as pd



pats = ['/Volumes/Drobo/macs_all_+pnx']

sample_list = []
sample_folder = []
for path_to_file in pats:
    for root, dirnames, filenames in os.walk(path_to_file):
        for f in filenames:
            if 'e_data.ctab' in f:
                sample_list.append(os.path.basename(root))
                sample_folder.append(root.split('/')[-2])
sample_df = pd.DataFrame.from_dict({'run':sample_list, 'folder':sample_folder})


sample_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_samples.txt'), sep='\t')
