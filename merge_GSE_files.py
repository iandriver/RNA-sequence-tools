import os
import fnmatch
import pandas as pd
import subprocess
import sys
import mygene



def ensembl_to_gene(df, index_col_name, species='human', ensembl_type='gene'):
    #ensembl_type: either 'transcript' or 'gene'
    #species: either 'mouse' or 'human'

    mg = mygene.MyGeneInfo()
    ensml_ids = df[index_col_name].tolist()
    if ensembl_type =='transcript':
        scope= 'ensembltranscript'
    elif ensembl_type == 'gene':
        scope = 'ensemblgene'
    sym_dict = mg.querymany(ensml_ids, scopes=scope, fields='symbol,name,alias,entrezgene,refseq', species=species, returnall=True, verbose=True, email='ian.driver@ucsf.edu')

    new_columns = [x for x in iso_de_df.columns.values]
    iso_de_df.columns = new_columns
    symbol_list = ['' for x in ensml_ids]
    alias_list = ['' for x in ensml_ids]
    id_num = [0 for x in ensml_ids]
    name_list = ['' for x in ensml_ids]
    entrez_list = ['' for x in ensml_ids]
    refseq_list = ['' for x in ensml_ids]

    for e_dicts in sym_dict['out']:
        if e_dicts['query'] in ensml_ids:
            index_num = ensml_ids.index(e_dicts['query'])
            try:
                id_num[index_num] = e_dicts['_id']
                good = True
            except KeyError:
                good =False
            if good:
                try:
                    symbol_list[index_num] = e_dicts['symbol']
                except KeyError:
                    symbol_list[index_num] = 'No_sym'
                try:
                    entrez_list[index_num] = e_dicts['entrezgene']
                except KeyError:
                    entrez_list[index_num] = ''
                try:
                    refseq_list[index_num] = e_dicts['refseq']
                except KeyError:
                    refseq_list[index_num] = ''
                try:
                    alias_list[index_num] = e_dicts['alias']
                except KeyError:
                    alias_list[index_num] = ''

                try:
                    name_list[index_num] = e_dicts['name']
                except KeyError:
                    name_list[index_num] = ''

    new_index = []
    new_name_dict = {}
    new_upper = {}
    dup = [item for item, count in collections.Counter(gene_df[index_col_name]).items() if count > 1]
    ensembl_sym_dict = dict(zip(ensml_ids, symbol_list))
    for tid, sym in zip(ensml_ids, symbol_list):
        new_index.append(sym+'_'+tid)
        if ensembl_type == 'transcript':
            new_name_dict[tid] = sym+'_'+tid
        elif ensembl_type == 'gene':
            if sym == '':
                new_name_dict[tid] = tid
            else:
                new_name_dict[tid] = sym
                new_upper[tid] = sym.upper()
    new_name_df = pd.DataFrame({'new_name': new_index, 'symbol':symbol_list,'alias':alias_list,'id_num':id_num, 'name':name_list, 'entrez gene': entrez_list, 'refseq': refseq_list}, index=ensml_ids)
    renamed_df = iso_de_df.rename(index=new_name_dict)
    renamed_upper_df = iso_de_df.rename(index=new_upper)
    renamed_df.to_csv(os.path.join(os.path.dirname(path_to_file),'sym_'+os.path.basename(path_to_file)), sep='\t')
    renamed_upper_df.to_csv(os.path.join(os.path.dirname(path_to_file),'sym_uppercase_genes_'+os.path.basename(path_to_file)), sep='\t')
    merge_df = pd.concat([new_name_df,iso_de_df], axis=1)
    new_name_df.to_csv(os.path.join(os.path.dirname(path_to_file),'ensml_to_sym_'+os.path.basename(path_to_file)), sep='\t')
    merge_df.to_csv(os.path.join(os.path.dirname(path_to_file),'ensml_to_sym_all_'+os.path.basename(path_to_file)), sep='\t')

def merge_datasets(list_of_paths):
    gene_dict_tpm = {}
    index_gene_old = []
    for path_to_file in paths:
        if os.path.isdir(path_to_file):
            make_matrix_from_dir(path_to_file)
        else:
            
def make_matrix_from_dir(path_to_file):
    ensembl_gene = False
    ensembl_transcript = False
    names_already_seen = []
    for root, dirnames, filenames in os.walk(path_to_file):
        gz_file=False
        for f in filenames:
            if '.gz' in f:
                gz_file =True
            if 'txt' in f or 'tsv' in f or 'csv' in f:
                if gz_file:
                    gzip_command = 'gunzip '+os.path.join(root,f)
                    subprocess.call(gzip_command, shell=True)
                    name = f[0:-7]
                    if name not in names_already_seen:
                        gene_df = pd.read_csv(os.path.join(root,f[0:-3]), sep=None, engine='python')
                        names_already_seen.append(name)
                        print(name)
                else:
                    if name not in names_already_seen:
                        gene_df = pd.read_csv(os.path.join(root,f), sep=None, engine='python')
                        names_already_seen.append(name)
                        print(name)
                cols_to_check = gene_df.columns.tolist()
                #check for headerless data (will have a gene value as a header)
                for col in cols_to_check:
                    try:
                        if isinstance(float(col), float):
                            gene_df = pd.read_csv(os.path.join(root,f), sep=None, engine='python', header = None)
                            gene_df.rename(columns={0:'tracking_id', 1:name}, inplace=True)
                    except ValueError:
                        pass
                cols_to_check = gene_df.columns.tolist()
                for col in cols_to_check:
                    if gene_df[col][0][0:4] == 'ENSG':
                        ensembl_gene = True
                        ensembl_gen_colname = col
                    elif gene_df[col][0][0:4] == 'ENST':
                        ensembl_transcript = True
                        ensembl_transcript_colname = col
                    elif 'EEF1A1' in gene_df[col].tolist() or 'GAPDH' in gene_df[col].tolist():

                index_gene_new = gene_df['tracking_id']
                if len(index_gene_old) != 0:
                    if index_gene_new.all() == index_gene_old.all():
                        pass
                    else:
                        sys.exit('Gene Indices do not match')
                else:
                    index_gene_old = index_gene_new
                gene_dict_tpm[name] = gene_df[name]
            elif 'txt' in f and 'gz' not in f:
                name = '_'.join(f.split('_')[1:6])
                print(name)
                gene_df = pd.read_csv(os.path.join(root,f), sep=None, engine='python')
                index_gene_new = gene_df['tracking_id']
                if len(index_gene_old) != 0:
                    if index_gene_new.all() == index_gene_old.all():
                        pass
                    else:
                        sys.exit('Gene Indices do not match')
                else:
                    index_gene_old = index_gene_new
                gene_dict_tpm[name] = gene_df[name]

gene_tpm_df = pd.DataFrame.from_dict(gene_dict_tpm)


gene_tpm_df.set_index(index_gene_new, inplace=True)



gene_tpm_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_gene_matrix.txt'), sep='\t')
