import os
import fnmatch
import pandas as pd



paths = ['/Volumes/Drobo/Seq_data/rsem_pdgfra_and_pteng']
gene_dict_tpm = {}
gene_dict_fpkm = {}
iso_dict_tpm = {}
iso_dict_fpkm = {}
sample_list = []
sample_folder = []
for path_to_file in paths:
    for root, dirnames, filenames in os.walk(path_to_file):
        for f in filenames:
            if '.genes.results' in f:
                name = os.path.basename(root)
                print(name)
                gene_df = pd.read_csv(os.path.join(root,f), sep=None, engine='python')
                index_gene = gene_df['gene_id']
                gene_dict_tpm[name] = gene_df['TPM']
                gene_dict_fpkm[name] = gene_df['FPKM']
                sample_list.append(os.path.basename(root))
                sample_folder.append(root)
            if '.isoforms.results' in f:
                name = os.path.basename(root)
                iso_df = pd.read_csv(os.path.join(root,f), sep=None, engine='python')
                index_iso = iso_df['transcript_id']
                iso_dict_tpm[name] = iso_df['TPM']
                iso_dict_fpkm[name] = iso_df['FPKM']
sample_df = pd.DataFrame.from_dict({'run':sample_list, 'folder':sample_folder})
gene_tpm_df = pd.DataFrame.from_dict(gene_dict_tpm)
gene_fpkm_df = pd.DataFrame.from_dict(gene_dict_fpkm)
iso_tpm_df = pd.DataFrame.from_dict(iso_dict_tpm)
iso_fpkm_df = pd.DataFrame.from_dict(iso_dict_fpkm)

gene_tpm_df.set_index(index_gene, inplace=True)
gene_fpkm_df.set_index(index_gene, inplace=True)
iso_tpm_df.set_index(index_iso, inplace=True)
iso_fpkm_df.set_index(index_iso, inplace=True)

sample_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_samples.txt'), sep='\t')
gene_tpm_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_gene_tpm_matrix.txt'), sep='\t')
gene_fpkm_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_gene_fpkm_matrix.txt'), sep='\t')
iso_tpm_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_isoform_tpm_matrix.txt'), sep='\t')
iso_fpkm_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_isoform_fpkm_matrix.txt'), sep='\t')
