import mygene
import pandas as pd
import os



path_to_file = 'DESeq_all_ctrl_vs_d4pnx_stats.txt'

#either 'transcript' or 'gene'
ensembl_type ='gene'

#species 'mouse' or 'human'
species ='mouse'

mg = mygene.MyGeneInfo()
iso_de_df = pd.read_csv(path_to_file, index_col=0, delimiter= '\t')
ensml_ids = iso_de_df.index.tolist()
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
for tid, sym in zip(ensml_ids, symbol_list):
    new_index.append(sym+'_'+tid)
    if ensembl_type == 'transcript':
        new_name_dict[tid] = sym+'_'+tid
    elif ensembl_type == 'gene':
        if sym == '':
            new_name_dict[tid] = tid
        else:
            new_name_dict[tid] = sym
new_name_df = pd.DataFrame({'new_name': new_index, 'symbol':symbol_list,'alias':alias_list,'id_num':id_num, 'name':name_list, 'entrez gene': entrez_list, 'refseq': refseq_list}, index=ensml_ids)
renamed_df = iso_de_df.rename(index=new_name_dict)
renamed_df.to_csv(os.path.join(os.path.dirname(path_to_file),'sym_'+os.path.basename(path_to_file)), sep='\t')
merge_df = pd.concat([new_name_df,iso_de_df], axis=1)
new_name_df.to_csv(os.path.join(os.path.dirname(path_to_file),'ensml_to_sym_'+os.path.basename(path_to_file)), sep='\t')
merge_df.to_csv(os.path.join(os.path.dirname(path_to_file),'ensml_to_sym_all_'+os.path.basename(path_to_file)), sep='\t')
