import os
from Bio import SeqIO
import sys
from Bio import Entrez
from sets import Set

# *Always* tell NCBI who you are
Entrez.email = "your_email@whatever"

#get_gene takes a list of gene names as symbols (set for mouse here) and returns a list of gene IDs
#that NCBI can use 
def get_gene(gene_list):
  gene_id_list = []
  for g_term in gene_list:
    g_handle = Entrez.esearch(db ='gene', term="Mus musculus[Orgn] AND "+g_term+'[sym]')
    g_record = Entrez.read(g_handle)
    gene_id = g_record["IdList"]
    gene_id_list.append(gene_id[0])
  return gene_id_list

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



gene_list = ['Notch3', 'Fgf10']

id_list = get_gene(gene_list)
print retrieve_annotation(id_list)
