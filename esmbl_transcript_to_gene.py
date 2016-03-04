import mygene
import pandas as pd

mg = mygene.MyGeneInfo()
iso_de_df = pd.read_csv('pdgfra_gene_matrix.txt', index_col=0, delimiter= '\t')
ensml_ids = iso_de_df.index.tolist()
sym_dict = mg.querymany(ensml_ids, scopes='ensemblgene', fields='symbol,name,alias,entrezgene,refseq', species='mouse', returnall=True)
new_columns = [x.split('/')[-2] for x in iso_de_df.columns.values]
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
        symbol_list[index_num] = e_dicts['symbol']
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
            print(e_dicts['symbol']+' has no alias')
        id_num[index_num] = e_dicts['_id']
        try:
            name_list[index_num] = e_dicts['name']
        except KeyError:
            name_list[index_num] = ''
            print(e_dicts['symbol']+' has no names')

new_index = []
for tid, sym in zip(ensml_ids, symbol_list):
    new_index.append(sym+'_'+tid)
new_name_df = pd.DataFrame({'new_name': new_index, 'symbol':symbol_list,'alias':alias_list,'id_num':id_num, 'name':name_list, 'entrez gene': entrez_list, 'refseq': refseq_list}, index=ensml_ids)

merge_df = pd.concat([new_name_df,iso_de_df], axis=1)
new_name_df.to_csv('ensml_to_sym.txt', sep='\t')
merge_df.to_csv('ensml_to_sym_de_gene_pdgfra.txt', sep='\t')
