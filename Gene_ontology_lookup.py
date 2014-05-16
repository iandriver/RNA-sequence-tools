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

'''Sample Output:
{'Fgf10': [{'Function': ['chemoattractant activity', 'type 2 fibroblast growth factor receptor binding', 'receptor binding', 'protein binding', 'heparin binding', 'fibroblast growth factor receptor binding', 'growth factor activity']}, {'Process': ['semicircular canal fusion', 'lung morphogenesis', 'lung development', 'branching morphogenesis of an epithelial tube', 'otic vesicle formation', 'negative regulation of cell proliferation', 'bud elongation involved in lung branching', 'positive regulation of mitotic cell cycle', 'embryonic camera-type eye development', 'regulation of smoothened signaling pathway', 'epithelial cell proliferation', 'regulation of activin receptor signaling pathway', 'embryonic genitalia morphogenesis', 'branch elongation involved in salivary gland morphogenesis', 'white fat cell differentiation', 'mesenchymal-epithelial cell signaling involved in lung development', 'inner ear morphogenesis', 'bronchiole morphogenesis', 'positive regulation of Wnt signaling pathway', 'negative regulation of extrinsic apoptotic signaling pathway in absence of ligand', 'epithelial cell migration', 'embryonic digestive tract development', 'lung alveolus development', 'positive regulation of canonical Wnt signaling pathway', 'positive regulation of epithelial cell proliferation', 'epithelial cell differentiation', 'mammary gland specification', 'lung proximal/distal axis specification', 'positive regulation of fibroblast proliferation', 'odontogenesis of dentin-containing tooth', 'lung saccule development', 'mammary gland bud formation', 'positive regulation of vascular endothelial growth factor receptor signaling pathway', 'somatic stem cell maintenance', 'ERK1 and ERK2 cascade', 'salivary gland development', 'embryonic digestive tract morphogenesis', 'prostatic bud formation', 'organ formation', 'positive regulation of ATPase activity', 'cell-cell signaling', 'thyroid gland development', 'organ morphogenesis', 'positive regulation of MAPK cascade', 'positive regulation of epithelial cell migration', 'spleen development', 'positive regulation of hair follicle cell proliferation', 'metanephros morphogenesis', 'lacrimal gland development', 'negative regulation of cell differentiation', 'activation of MAPK activity', 'urothelial cell proliferation', 'submandibular salivary gland formation', 'positive regulation of transcription from RNA polymerase II promoter', 'organ growth', 'keratinocyte proliferation', 'positive regulation of white fat cell proliferation', 'fibroblast growth factor receptor signaling pathway', 'positive regulation of DNA replication', 'positive regulation of keratinocyte proliferation', 'digestive tract development', 'blood vessel remodeling', 'female genitalia morphogenesis', 'induction of positive chemotaxis', 'positive regulation of transcription, DNA-templated', 'fibroblast growth factor receptor signaling pathway involved in mammary gland specification', 'epithelial cell proliferation involved in salivary gland morphogenesis', 'male genitalia morphogenesis', 'epidermis morphogenesis', 'wound healing', 'mesenchymal cell differentiation involved in lung development', 'determination of left/right symmetry', 'limb bud formation', 'positive regulation of urothelial cell proliferation', 'positive regulation of lymphocyte proliferation', 'respiratory system development', 'pancreas development', 'positive regulation of peptidyl-tyrosine phosphorylation', 'radial glial cell differentiation', 'positive regulation of cyclin-dependent protein serine/threonine kinase activity involved in G1/S transition of mitotic cell cycle', 'regulation of gene expression', 'Harderian gland development', 'pituitary gland development', 'positive regulation of Notch signaling pathway', 'positive regulation of cell proliferation', 'muscle cell fate commitment', 'tear secretion', 'limb development', 'hair follicle morphogenesis', 'branching involved in salivary gland morphogenesis', 'actin cytoskeleton reorganization', 'positive regulation of keratinocyte migration', 'limb morphogenesis', 'secretion by lung epithelial cell involved in lung growth', 'salivary gland morphogenesis', 'epithelial tube branching involved in lung morphogenesis', 'semicircular canal morphogenesis', 'organ induction', 'establishment of mitotic spindle orientation', 'regulation of epithelial cell proliferation', 'negative regulation of cell cycle arrest', 'regulation of branching involved in salivary gland morphogenesis by mesenchymal-epithelial signaling', 'positive chemotaxis', 'protein localization to cell surface', 'embryonic pattern specification', 'blood vessel morphogenesis', 'bud outgrowth involved in lung branching', 'smooth muscle cell differentiation', 'positive regulation of Ras protein signal transduction', 'positive regulation of DNA repair', 'lung epithelium development', 'Type II pneumocyte differentiation', 'thymus development', 'chemotaxis', 'epidermis development', 'positive regulation of ERK1 and ERK2 cascade', 'regulation of saliva secretion']}, {'Component': ['extracellular space', 'extracellular matrix', 'extracellular region', 'cell', 'nucleus', 'cell surface', 'plasma membrane']}], 'Notch3': [{'Function': ['calcium ion binding', 'enzyme binding', 'protein binding']}, {'Process': ['multicellular organismal development', 'negative regulation of cell differentiation', 'forebrain development', 'regulation of transcription, DNA-templated', 'cell differentiation', 'negative regulation of neuron differentiation', 'neuron fate commitment', 'Notch signaling pathway', 'transcription, DNA-templated', 'positive regulation of smooth muscle cell proliferation', 'regulation of developmental process']}, {'Component': ['integral component of plasma membrane', 'receptor complex', 'nucleus', 'membrane', 'integral component of membrane', 'plasma membrane']}]}'''
