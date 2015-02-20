import os
import cPickle as pickle
import pandas as pd

cuff_df = pd.DataFrame.from_csv('/Volumes/Seq_data/results_macs_pnx/sc_expr_out_macs/genes.fpkm_table', sep='\t')

with open('/Volumes/Seq_data/results_macs_pnx/cuff_fpkm_table.p', 'wb') as fp:
  pickle.dump(cuff_df, fp)
