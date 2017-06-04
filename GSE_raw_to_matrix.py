import os
import fnmatch
import pandas as pd
import subprocess
import sys
import mygene
import collections


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
        ensembl_to_sym_dict[e] = g.gene_name
        sym_gene_list.append(g.gene_name)
    df.drop(no_ensembl_list, inplace=True)
    renamed_df = df.rename(index=ensembl_to_sym_dict)
    return(renamed_df)

def ensembl_to_sym(iso_de_df, ensembl_type ='gene',species ='human'):
    #ensembl_type either 'transcript' or 'gene'
    #species 'mouse' or 'human'


    mg = mygene.MyGeneInfo()
    ensml_ids = iso_de_df.index.tolist()
    print(len(ensml_ids))
    if ensembl_type =='transcript':
        scope= 'ensembltranscript'
    elif ensembl_type == 'gene':
        scope = 'ensemblgene'
    sym_dict = mg.querymany(ensml_ids, scopes=scope, fields='symbol,name,alias,entrezgene,refseq', species=species, returnall=True, verbose=True, email='ian.driver@ucsf.edu')
    print(len(sym_dict))
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

    enseml_sym_dict = dict(zip(ensml_ids, symbol_list))
    new_index = []
    new_name_dict = {}
    for tid, sym in zip(ensml_ids, symbol_list):
        new_index.append(sym+'_'+tid)
        if ensembl_type == 'transcript':
            new_name_dict[tid] = sym+'_'+tid
        elif ensembl_type == 'gene':
            if sym == '':
                new_name_dict[tid] = tid
            else:
                new_name_dict[tid] = sym
    renamed_df = iso_de_df.rename(index=new_name_dict)
    return(renamed_df)

def merge_raw(paths, count_or_RPKM='RPKM'):
    gene_dict_tpm = {}
    index_gene_old = []
    for path_to_file in paths:
        for root, dirnames, filenames in os.walk(path_to_file):
            for f in filenames:
                if '.gz' in f and count_or_RPKM in f and 'matrix' not in f:
                    gzip_command = 'gunzip '+os.path.join(root,f)
                    name = f[0:-7]
                    subprocess.call(gzip_command, shell=True)
                    print(name)
                    gene_df = pd.read_csv(os.path.join(root,f[0:-3]), sep=None, engine='python')
                    try:
                        index_gene_new = gene_df['tracking_id']
                    except KeyError:
                        index_gene_new = gene_df['ensG']
                        ensembl= True
                    if len(index_gene_old) != 0:
                        if index_gene_new.all() == index_gene_old.all():
                            pass
                        else:
                            sys.exit('Gene Indices do not match')
                    else:
                        index_gene_old = index_gene_new
                    try:
                        gene_dict_tpm[name] = gene_df[name]
                    except KeyError:
                        gene_dict_tpm[name] = gene_df[count_or_RPKM]
                elif count_or_RPKM in f and 'gz' not in f and 'matrix' not in f:
                    name = f[0:-4]
                    print(name)
                    gene_df = pd.read_csv(os.path.join(root,f), sep=None, engine='python')
                    try:
                        index_gene_new = gene_df['tracking_id']
                    except KeyError:
                        index_gene_new = gene_df['ensG']
                        ensembl= True
                    if len(index_gene_old) != 0:
                        if index_gene_new.all() == index_gene_old.all():
                            pass
                        else:
                            sys.exit('Gene Indices do not match')
                    else:
                        index_gene_old = index_gene_new
                    try:
                        gene_dict_tpm[name] = gene_df[name]
                    except KeyError:
                        gene_dict_tpm[name] = gene_df[count_or_RPKM]


    gene_tpm_df = pd.DataFrame.from_dict(gene_dict_tpm)

    gene_tpm_df.set_index(index_gene_new, inplace=True)

    if ensembl:
        gene_df_collapsed = gene_tpm_df.groupby('ensG').agg(['sum'])
        gene_df_collapsed.columns = gene_df_collapsed.columns.get_level_values(0)
        gene_df_collapsed.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_ensembl_gene_'+count_or_RPKM+'_matrix.txt'), sep='\t')
        gene_tpm_df_new = ensembl_to_sym2(gene_df_collapsed)
        gene_tpm_df_new.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_gene_'+count_or_RPKM+'_matrix.txt'), sep='\t')
    else:
        gene_tpm_df.to_csv(os.path.join(path_to_file,os.path.basename(path_to_file)+'_gene_'+count_or_RPKM+'_matrix.txt'), sep='\t')


merge_raw(['/Users/iandriver/Downloads/GSE73727_RAW'])
