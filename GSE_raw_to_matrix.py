import os
import fnmatch
import pandas as pd
import subprocess
import sys
import mygene
import collections
import numpy as np

#path to GSE RAW file series or single cell data matrices to merge
files_and_matrix_list = ['/Users/idriver/Downloads/GSE71318_RAW','/Users/idriver/Downloads/GSE73727_RAW','/Users/idriver/Downloads/GSE63818_RAW','/Users/idriver/Downloads/GSE66507_human_blastocyst_rnaseq_counts.txt', '/Users/idriver/Downloads/GSE94820_raw.expMatrix_deeper.characterization.set.submission.txt','/Users/idriver/Downloads/GSE72056_melanoma_single_cell_revised_v2_2.txt']

def ensembl_to_sym2(df):
    from pyensembl import EnsemblRelease
    data = EnsemblRelease(87)
    ensml_ids = df.index.tolist()
    sym_gene_list = []
    ensembl_to_sym_dict = {}
    no_ensembl_list = []
    for e in ensml_ids:
        try:
            g = data.gene_by_id(e)
        except ValueError:
            no_ensembl_list.append(e)
        if g.gene_name.strip() not in sym_gene_list:
            ensembl_to_sym_dict[e] = g.gene_name.strip()
            sym_gene_list.append(g.gene_name.strip())
    df.drop(no_ensembl_list, inplace=True)
    renamed_df = df.rename(index=ensembl_to_sym_dict)
    renamed_df_grouped = renamed_df.groupby(renamed_df.index).first()
    numeric_df = renamed_df_grouped.convert_objects(convert_numeric=True)
    return(numeric_df)

def merge_datasets(list_of_paths):
    cell_data_dict = collections.OrderedDict()
    count_files = 0
    for path_to_file in list_of_paths:
        if os.path.isdir(path_to_file):
            if os.path.basename(path_to_file) == 'GSE63818_RAW':
                file_df_current = make_matrix_from_dir(path_to_file, count_or_RPKM='expression')
            else:
                file_df_current = make_matrix_from_dir(path_to_file)
            cell_set = os.path.basename(path_to_file).split('_')[0]
            for cols in file_df_current.columns.tolist():
                cell_data_dict[cols] = cell_set
        else:
            file_df_current = pd.read_csv(path_to_file, sep=None, engine='python', index_col=0)
            if isinstance(file_df_current.index.tolist()[0], str):
                if file_df_current.index.tolist()[0][0:4] =='ENSG':
                    file_df_current = ensembl_to_sym2(file_df_current)
            cell_set = os.path.basename(path_to_file).split('_')[0]
            for cols in file_df_current.columns.tolist():
                cell_data_dict[cols] = cell_set
        file_df_current = file_df_current.groupby(file_df_current.index).first()
        if count_files == 0:
            index1 = file_df_current.index.tolist()
            merge_df1 = file_df_current.copy()
            count_files += 1
        elif count_files > 0:
            index2 = file_df_current.index.tolist()
            if index1 != index2:
                print(len(index1), len(index2))
                merge_index = list(set(index1)&set(index2))
                merge_index.sort()
                print(len(merge_index))
                merge_df_new = merge_df1.loc[merge_index, :].merge(file_df_current.loc[merge_index, :], how='left', left_index=True, right_index=True)
            else:
                merge_df_new = merge_df1.merge(file_df_current, how='left', left_index=True, right_index=True)
            index1 = merge_index.copy()
            merge_df1 = merge_df_new.copy()
            merge_df1 = merge_df1.groupby(merge_df1.index).first()
            count_files += 1
        print(path_to_file, merge_df1.shape)

    cell_data_df = pd.DataFrame.from_dict(cell_data_dict,orient='index')
    merge_df1.fillna(value=0, inplace=True)
    try:
        merge_df1_to_save = merge_df1.astype(np.float64)
    except ValueError:
        merge_df1.replace({"x1f|\x02|\x12|\x0f|\x16|\x03(?:\d{1,2}(?:,\d{1,2})?)?": ''}, regex=True, inplace=True)
    merge_df1_to_save = merge_df1.astype(np.float64)
    merge_df1_to_save.to_csv('merged_gene_matrix_all_v1.txt', sep='\t')
    cell_data_df.to_csv('merged_cell_data_source_all_v1.txt', sep='\t')

def make_matrix_from_dir(path_to_dir, count_or_RPKM='RPKM'):
    gene_dict_tpm = {}
    index_gene_old = []
    index_gene_new = []
    ensembl=False
    for root, dirnames, filenames in os.walk(path_to_dir):
        for f in filenames:
            if '.gz' in f and count_or_RPKM in f:
                gz =True
                f_name =  f[0:-3]
                name = f[0:-7]
                gzip_command = 'gunzip '+os.path.join(root,f)
                subprocess.call(gzip_command, shell=True)
                file_found = True
            elif count_or_RPKM in f and 'gz' not in f and 'matrix' not in f:
                name = f[0:-4]
                f_name =  f
                file_found = True
            else:
                file_found = False
            if file_found:
                gene_df = pd.read_csv(os.path.join(root,f_name), sep=None, engine='python')
                cols_to_check = gene_df.columns.tolist()
                #check for headerless data (will have a gene value as a header)
                try:
                    if isinstance(float(cols_to_check[1]), float):
                        gene_df = pd.read_csv(os.path.join(root,f_name), sep=None, engine='python', header = None)
                        gene_df.rename(columns={0:'tracking_id', 1:name}, inplace=True)
                except ValueError:
                    pass
                cols_to_check = gene_df.columns.tolist()
                for col in cols_to_check:
                    if isinstance(gene_df[col][0], str):
                        if gene_df[col][0][0:4] == 'ENSG':
                            ensembl = True
                            ensembl_gen_colname = col
                            col_name = col

                        elif 'EEF1A1' in gene_df[col].tolist() or 'GAPDH' in gene_df[col].tolist():
                            ensembl = False
                            col_name = col

                try:
                    gene_df['tracking_id'] = gene_df['tracking_id'].str.strip()
                    index_gene_new = gene_df['tracking_id'].copy()
                except KeyError:
                    index_gene_new = gene_df['ensG'].copy()
                    ensembl= True
                if len(index_gene_old) != 0:
                    if index_gene_new.all() == index_gene_old.all():
                        pass
                    else:
                        sys.exit('Gene Indices do not match')
                else:
                    index_gene_old = index_gene_new
                col_tries = [c for c in [name, count_or_RPKM, cols_to_check[1]] if c not in [col_name, 'ensT']]
                stop = False
                for cols in col_tries:
                    try:
                        gene_dict_tpm[name] = gene_df[cols]
                        stop= True
                    except KeyError:
                        stop = False
                        pass
                    if stop:
                        break


    gene_tpm_df = pd.DataFrame.from_dict(gene_dict_tpm)
    if len(index_gene_new) != 0:
        gene_tpm_df.set_index(index_gene_new, inplace=True)


    if ensembl:
        gene_df_collapsed = gene_tpm_df.groupby('ensG').agg(['sum'])
        gene_df_collapsed.columns = gene_df_collapsed.columns.get_level_values(0)
        gene_df_collapsed.to_csv(os.path.join(path_to_dir,os.path.basename(path_to_dir)+'_ensembl_gene_'+count_or_RPKM+'_matrix.txt'), sep='\t')
        gene_tpm_df_new = ensembl_to_sym2(gene_df_collapsed)
        gene_tpm_df_new.sort_index(inplace=True)
        gene_tpm_df_new.to_csv(os.path.join(path_to_dir,os.path.basename(path_to_dir)+'_gene_'+count_or_RPKM+'_matrix.txt'), sep='\t')
        final_df = gene_tpm_df_new.copy()
    else:
        gene_tpm_df.sort_index(inplace=True)
        gene_tpm_df.to_csv(os.path.join(path_to_dir,os.path.basename(path_to_dir)+'_gene_'+count_or_RPKM+'_matrix.txt'), sep='\t')
        final_df = gene_tpm_df.copy()

    return(final_df)

merge_datasets(files_and_matrix_list)
