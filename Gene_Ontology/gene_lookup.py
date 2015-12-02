import os
from Bio import SeqIO
import cPickle as pickle
import sys
from Bio import Entrez
from sets import Set
import json
from pprint import pprint
import pandas as pd
from collections import OrderedDict
from itertools import izip
#Basic script for fetching gene ontology from Entrez
#Starting with common gene names in a list 1) get_gene takes the list and returns gene ids that are suitable
#for Entrez. 2) retrieve_annotation takes a gene id list and returns the gene ontology in a dict of dicts

# *Always* tell NCBI who you are
Entrez.email = "ian.driver@ucsf.edu"

#the file path where gene list will be and where new list will output
path_to_file = '/Volumes/Seq_data/cuffnorm_pdgfra_ctrl_only'
species = 'mouse'
if species == 'human':
    entrez_species = 'Homo sapiens'
    go_json='hu_gene_go.json'
elif species == 'mouse':
    entrez_species = 'Mus musculus'
    go_json='gene_gos.json'
#where the json gene ontology database is stored
path_to_gojson ="/Volumes/Seq_data"
#name of file containing gene
#the file path where gene list will be and where new list will output
use_gene_file = False
#name of file containing gene
gene_file_source = 'go_search_genes_lung_all.txt'
#if you want to update the database change to True
update = False
base_name = 'pdgfra_ctrl_only'

#load file gene
by_cell = pd.DataFrame.from_csv(os.path.join(path_to_file, base_name+'_outlier_filtered.txt'), sep='\t')
by_gene = by_cell.transpose()

def return_pickle(path_to_file, file):
    fpgenelist = open(os.path.join(path_to_file,file), 'rb')
    gene_list = pickle.load(fpgenelist)
    return gene_list
    fpgenelist.close()


def singular_gene_list(path_to_file, file):
  gene_df = pd.read_csv(os.path.join(path_to_file,file), delimiter= '\t')
  gene_list = gene_df['GeneID'].tolist()
  final_list = []
  for g in gene_list:
    final_list.append(g.split('_')[0])
  return final_list

#create list of genes
if use_gene_file:
    #will only take pickle files and .txt singular files (has [GeneID] column)
    if gene_file_source[-1] != 'p':
        g_list = return_pickle(path_to_file, gene_file_source)
    else:
        g_list = singular_gene_list(path_to_file, gene_file_source)
else:
    g_list = by_cell.index.tolist()

def check_ids(id_list, filename=os.path.join(path_to_gojson, go_json)):
    already_known = []
    if os.path.isfile(filename):
        with open(filename, 'r') as gg2:
            go_json = json.load(gg2)
        for gg_entry in go_json:
            if gg_entry in id_list:
                already_known.append(id_list.pop(id_list.index(gg_entry)))
        new_id_list = id_list
    else:
        new_id_list = id_list
        already_known =[]
    return new_id_list, already_known

#get_gene takes a list of gene names as symbols (set for mouse here) and returns a list of gene IDs
#that NCBI can use
def get_gene(gene_list):
    new_gene_list, already_known = check_ids(gene_list)
    gene_id_list = []
    for g_term in new_gene_list:
        g_handle = Entrez.esearch(db ='gene', term=entrez_species+"[Orgn] AND "+g_term+'[sym]')
        g_record = Entrez.read(g_handle)
        gene_id = g_record["IdList"]
        if gene_id != []:
            gene_id_list.append(gene_id[0])
        else:
            print g_term+' not in Entrez'
    return gene_id_list, already_known

#retrieve_annotation takes a list of gene IDs and returns all of the associated
#Gene Ontology terms as a dictionary
def retrieve_annotation(id_list):
  goterms = OrderedDict()
  handle = Entrez.efetch(db='gene', id=",".join(id_list), retype='gb', retmode='xml')
  all_records = Entrez.parse(handle, 'genebank')
  for record in all_records:
    try:
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
                gotype= OrderedDict()
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
    except KeyError:
        pass
  return goterms

def update_json(gene_list, filename=os.path.join(path_to_gojson, go_json)):
    id_list, already_known = get_gene(gene_list)
    if id_list != []:
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
            data = gos
            with open(filename, 'w') as f2:
                json.dump(data, f2)

def return_json(gene_list, filename=os.path.join(path_to_gojson, go_json)):
    id_list, already_known = get_gene(gene_list)
    if id_list != []:
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
            data = gos
            with open(filename, 'w') as f2:
                json.dump(data, f2)
    else:
        with open(filename, 'rw') as gg2:
            go_json_file = json.load(gg2)
        for g in already_known:
            pprint(go_json_file[g])

def compress(data, selectors):
    # compress('ABCDEF', [1,0,1,0,1,1]) --> A C E F
    return [d for d, s in izip(data, selectors) if s]

if update:
    update_json(g_list)
with open(os.path.join(path_to_gojson, go_json), 'rw') as gg2:
    go_json = json.load(gg2)
go_search_term =[('Process', 'Type II pneumocyte differentiation'),
                ('Process', 'lung alveolus development'),
                ('Process', 'lung development'),
                ('Process', 'Clara cell differentiation'),
                ('Process', 'lung saccule development'),
                ('Process', 'lung epithelium development'),
                ('Process', 'mesenchymal-epithelial cell signaling'),
                ('Process', 'establishment of planar polarity'),
                ('Process', 'keratinocyte differentiation'),
                ('Process', 'negative regulation of cellular senescence'),
                ('Process', 'epithelial cell development'),
                ('Process', 'positive regulation of mesenchymal cell proliferation'),
                ('Process', 'Notch signaling pathway'),
                ('Process', 'Wnt signaling pathway'),
                ('Process', 'negative regulation of fibroblast growth factor receptor signaling pathway'),
                ('Process', 'negative regulation of epithelial cell proliferation',),
                ('Process', 'non-canonical Wnt signaling pathway via JNK cascade'),
                ('Process', 'canonical Wnt signaling pathway'),
                ('Process', 'Wnt signaling pathway, calcium modulating pathway'),
                ('Process', 'activation of JUN kinase activity'),
                ('Process', 'negative regulation of epithelial to mesenchymal transition'),
                ('Process', 'negative regulation of transforming growth factor beta receptor signaling pathway'),
                ('Process', 'negative regulation of apoptotic process'),
                ('Process', 'proximal/distal pattern formation'),
                ('Component', 'alveolar lamellar body'),
                ('Component', 'alveolar lamellar body membrane'),
                ('Function', 'cytokine activity'),
                ('Function', 'NF-kappaB binding'),
                ('Process', 'cilium assembly'),
                ('Component', 'integral component of membrane'),
                ('Component', 'external side of plasma membrane'),
                ('Component', 'keratin filament'),
                ('Function', 'sequence-specific DNA binding transcription factor activity'),
                ('Function', 'transcription regulatory region sequence-specific DNA binding'),
                ('Function', 'transcription factor activity, sequence-specific DNA binding'),
                ('Process','positive regulation of fibroblast migration'),
                ('Process', 'positive regulation of cell migration'),
                ('Process','tumor necrosis factor-mediated signaling pathway'),
                ('Process', 'negative regulation of inflammatory response'),
                ('Process', 'regulation of inflammatory response'),
                ('Process', 'activation of MAPK activity'),
                ('Process', 'cell fate commitment'),
                ('Process', 'midgut development'),
                ('Process', 'positive regulation of type I interferon-mediated signaling pathway',),
                ('Process', 'negative regulation of Notch signaling pathway'),
                ('Process', 'positive regulation of interleukin-6 production'),
                ('Process', 'positive regulation of NF-kappaB transcription factor activity'),
                ('Process', 'positive regulation of cytokine secretion involved in immune response'),
                ('Process', 'establishment of planar polarity',),
                ('Component', 'extracellular matrix'),
                ('Process', 'lipid storage'),
                ('Process', 'liver development'),
                ('Process', 'fibroblast growth factor receptor signaling pathway'),
                ('Process', 'kidney development'),
                ('Function','extracellular matrix structural constituent'),
                ('Process', 'regulation of MAPK cascade'),
                ('Function', 'calcium channel activity'),
                ('Process', 'negative regulation of apoptotic process'),
                ('Process', 'ion transmembrane transport'),
                ('Function', 'SMAD binding'),
                ('Function', 'transcription corepressor activity'),
                ('Process', 'negative regulation of neuron apoptotic process'),
                ('Process', 'negative regulation of JAK-STAT cascade'),
                ('Process', 'negative regulation of inflammatory response'),
                ('Process', 'cell migration'),
                ('Process', 'response to hyperoxia'),
                ('Function', 'phospholipase activity'),
                ('Process', 'lung saccule development'),
                ('Process', 'Clara cell differentiation'),
                ('Component', 'integral component of plasma membrane'),
                ('Function', 'RNA binding'),
                ('Function', 'protein binding',),
                ('Function', 'chromatin binding'),
                ('Process', 'gene expression'),
                ('Component', 'cytosol'),
                ('Component', 'cytoplasm'),
                ('Component', 'extracellular space'),
                ('Component', 'extracellular region'),
                ('Component', 'nucleus'),
                ('Component', 'plasma membrane'),
                ('Component', 'membrane')]
term_index =['Function', 'Process', 'Component']
search_term_dict =OrderedDict()
search_term_list = []
search_term_dict['GeneID'] = []
search_term_dict['GroupID'] = []

for go_term in go_search_term:
    for g in g_list:
        if g == 'Ly6a':
            print go_json[g][0]['Component'], len(go_json[g])
        gene_tries = []
        gene_mask = []
        try:
            for parts in go_json[g]:
                try:
                    gene_test = parts[go_term[0]]
                    gene_tries.append(gene_test)
                    gene_mask.append(True)
                except:
                    gene_tries.append([])
                    gene_mask.append(False)
                    pass
            if True in gene_mask:
                print gene_mask, gene_tries
                gene = compress(gene_tries, gene_mask)[0]
                print gene
                if go_term[1] in gene:
                    if g not in search_term_dict['GeneID']:
                        search_term_list.append(g)
                        search_term_dict['GeneID'].append(g)
                        search_term_dict['GroupID'].append(go_term[1])
                    else:
                        pass
        except KeyError:
            pass

searches_df = pd.DataFrame(search_term_dict)
searches_df.to_csv(os.path.join(path_to_file, 'go_search_genes_lung_all.txt'), sep = '\t', index=False)
