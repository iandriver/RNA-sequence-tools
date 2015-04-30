import os
import cPickle as pickle
import pandas as pd



path_to_file= '/Volumes/Seq_data/counts_sheppard_all/cpm_norm_liver_kidney_1'
file_name = 'sheppard_normalized_cpm_liver_kidney.txt'
cuff_df = pd.DataFrame.from_csv(os.path.join(path_to_file,file_name), sep='\t')

with open(os.path.join(path_to_file,file_name.strip('.txt')+'.p'), 'wb') as fp:
  pickle.dump(cuff_df, fp)
