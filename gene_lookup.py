import os
from Bio import SeqIO
import sys
from Bio import Entrez
from sets import Set
import fnmatch
import collections
import cPickle as pickle
import pandas as pd

# *Always* tell NCBI who you are
Entrez.email = "ian.driver@ucsf.edu"
p_name = 'gene_annotation.p'
def make_pickle(file_to_dump, p_name='gene_annotation.p'):
  with open(p_name, 'wb') as fp:
    pickle.dump(file_to_dump, fp)

def load_pickle(p_name='gene_annotation.p'):
  with open(p_name, 'rb') as fp:
    data = pickle.load(fp)
  return data

def get_gene(gene_list):
  new_gene_id_list = []
  arch_gene_id = []
  archive = load_pickle(p_name='gene_list.p')
  for g_term in gene_list:
    if g_term not in archive:
      archive.append(g_term)
      g_handle = Entrez.esearch(db ='gene', term="Mus musculus[Orgn] AND "+g_term+'[sym]')
      g_record = Entrez.read(g_handle)
      gene_id = g_record["IdList"]
      if len(gene_id)>0:
       new_gene_id_list.append(gene_id[0])
    else:
      arch_gene_id.append(archive[g_term])
  make_pickle(archive, p_name='gene_list.p')
  return new_gene_id_list, arch_gene_id

def flatten(d, parent_key=''):
  items = []
  for k, v in d.items():
    new_key = parent_key + '_' + k if parent_key else k
    if isinstance(v, collections.MutableMapping):
      items.extend(flatten(v, new_key).items())
    else:
      items.append((new_key, v))
  return dict(items)

def retrieve_annotation(id_list):
  new_genes, archived_genes = get_gene(gene_list)
  goterms = {}
  handle = Entrez.efetch(db='gene', id=",".join(new_genes), retype='gb', retmode='xml')
  all_records = Entrez.parse(handle, 'genebank')
  for record in all_records:
    for rec in record['Entrezgene_properties']:
      for i, k in rec.items():
        if i == 'Gene-commentary_properties':
          for n in k:
            sym_check = False
            for ni1, n1 in n.items():
              if n1 == 'Official Symbol':
                sym_check = True
              if ni1 =='Gene-commentary_text' and sym_check:
                #Finds the gene name (symbol or short version) and adds an empty list
                #in the goterms dictionary.  Changes sym_check to false so only
                #the one gene name is created for that gene
                gene_name = n1
                goterms[gene_name]= []
                sym_check = False
        if i == 'Gene-commentary_comment':
          for k2 in k:
            #there are three gene ontology types in NCBI gene database:
            #Function, Process, and Component.  The gotype dict will
            #store each go term under its type
            gotype= {}
            for i3, k3 in k2.items():
              if i3 == 'Gene-commentary_label':
                gotype_name = k3
                gotype[gotype_name] = []
                goterms[gene_name].append(gotype)
              if i3 == 'Gene-commentary_comment':
                for k4 in k3:
                  for i5, k5 in k4.items():
                    if i5 == 'Gene-commentary_source':
                      for k6 in k5:
                        for i7, k7 in k6.items():
                          if i7 == 'Other-source_anchor':
                            gotype[gotype_name].append(k7)
            #reduce the Go terms to a set to remove dublicate entries
            set_type = list(Set(gotype[gotype_name]))
            gotype[gotype_name]= set_type
            #only add the gotype once
            if gotype not in goterms[gene_name]:
              goterms[gene_name].append(gotype)
  #this function returns a dictionary of gene names each with a list of dicts
  #with GO terms by type example {gene1: [{Function_type1: [GOterm1, GOterm2]},
  # {Function_type2: [GOterm3, GOterm4]}]}
  return goterms


def chop_genelist(path='/Users/idriver/RockLab-files/test'):
  g = 0
  for root, dirnames, filenames in os.walk(path):
    for filename in fnmatch.filter(filenames, '*.fpkm_tracking'):
      if '_genes' in filename:
        curr_cell_fpkm_g =[]
        curr_g_file = open(os.path.join(root, filename))
        g_cell_name = (root.split('/')[-1])
        if g == 0:
          curr_cell_genes = []
        for k, line in enumerate(curr_g_file):
          if k == 0:
            header = line.strip('\n').split('\t')
          if k > 0:
            curr_g_line = line.strip('\n').split('\t')
            #exclude RNA spike controls
            if g == 0 and curr_g_line[3][:4] != 'ERCC':
              curr_cell_genes.append(curr_g_line[0])
        g = 1
        break
  sp_gene_list = []
  x_last =0
  for x in range(999, len(curr_cell_genes), 999):
    sp_gene_list.append(curr_cell_genes[x_last:x])
    x_last = x
  return sp_gene_list



def file_gene_list(path='/Users/idriver/RockLab-files/spc_ensembl/'):
  file = 'new-spc-highPC1-2.txt'
  gene_df = pd.read_csv(path+file, delimiter= '\t')
  gene_list = gene_df['GeneID'].tolist()
  final_list = []
  for g in gene_list:
    final_list.append(g.split('_')[0])
  return final_list

def print_annotation(ann_list, gene_name, sub_part='all'):
  for s in ann_list.items():
    if s[0] in gene_name:
      print s[0]
      for s2 in s[1]:
        if sub_part=='all':
          print s2
        for k, v in s2.items():
          if sub_part=k:
            print k, v

def get_p_gene(ann_list):
  genes = []
  for s in ann_list.items():
    genes.append(s[0])

gene_list = file_gene_list()
answer = retrieve_annotation(get_gene(gene_list))
get_annotation(answer,gene_list)
