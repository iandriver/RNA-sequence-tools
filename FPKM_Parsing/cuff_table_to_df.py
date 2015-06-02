import os
import cPickle as pickle
import pandas as pd



path_to_file= '/Volumes/Seq_data/cuffnorm_chapman_hu_all'
file_name = 'genes.fpkm_table'
cuff_df = pd.DataFrame.from_csv(os.path.join(path_to_file,file_name), sep='\t')

with open(os.path.join(path_to_file,'chapman_hu_all.p'), 'wb') as fp:
  pickle.dump(cuff_df, fp)
