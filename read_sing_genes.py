import pandas as pd

def file_gene_list(path='/Users/idriver/RockLab-files/spc_ensembl/'):
  file = 'new-spc-highPC1-2.txt'
  gene_df = pd.read_csv(path+file, delimiter= '\t')
  gene_list = gene_df['GeneID']
  return gene_list
