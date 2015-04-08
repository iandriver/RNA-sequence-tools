import os
from Bio import SeqIO
import sys
from Bio import Entrez
from sets import Set
import json
from pprint import pprint

#Basic script for fetching gene ontology from Entrez
#Starting with common gene names in a list 1) get_gene takes the list and returns gene ids that are suitable
#for Entrez. 2) retrieve_annotation takes a gene id list and returns the gene ontology in a dict of dicts

# *Always* tell NCBI who you are
Entrez.email = "whoever@where.edu"

path_to_gojson ="/path/to/go_json_file

def file_gene_list(path='/Users/idriver/RockLab-files/pdgfra/', file = 'pdgfra-ly6a-dcn-group.txt'):
  gene_df = pd.read_csv(path+file, delimiter= '\t')
  gene_list = gene_df['GeneID'].tolist()
  final_list = []
  for g in gene_list:
    final_list.append(g.split('_')[0])
  return final_list

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

def check_ids(id_list, filename='gene_gos.json'):
    already_known = []
    with open(filename, 'r') as gg2:
        go_json = json.load(gg2)
    for gg_entry in go_json:
        if gg_entry in id_list:
            already_known.append(id_list.pop(id_list.index(gg_entry)))
    new_id_list = id_list
    return new_id_list, already_known

#get_gene takes a list of gene names as symbols (set for mouse here) and returns a list of gene IDs
#that NCBI can use
def get_gene(gene_list):
  new_gene_list, already_known = check_ids(gene_list)
  gene_id_list = []
  for g_term in new_gene_list:
    g_handle = Entrez.esearch(db ='gene', term="Mus musculus[Orgn] AND "+g_term+'[sym]')
    g_record = Entrez.read(g_handle)
    gene_id = g_record["IdList"]
    gene_id_list.append(gene_id[0])
  return gene_id_list, already_known
            
#retrieve_annotation takes a list of gene IDs and returns all of the associated
#Gene Ontology terms as a dictionary
def retrieve_annotation(id_list):
  goterms = {}
  handle = Entrez.efetch(db='gene', id=",".join(id_list), retype='gb', retmode='xml')
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

def update_json(gene_list, filename=os.path.join(path_to_gojson, 'gene_gos.json'):
    id_list, already_known = get_gene(gene_list)
    GO_json = json.dumps(retrieve_annotation(id_list))
    gos = json.loads(GO_json)
    if os.path.isfile(filename):
        with open(filename,'r') as f:
            data = json.load(f)

        data.update(gos)

        with open(filename, 'w') as f2:
            json.dump(data, f2)
    else:
        print "Creating new go JSON: "+filename
        with open(filename, 'w') as f2:
            json.dump(data, f2)


gene_list = ['Dcn']

update_json(gene_list)

with open('gene_gos.json', 'r') as gg2:
    go_json = json.load(gg2)
pprint(go_json)