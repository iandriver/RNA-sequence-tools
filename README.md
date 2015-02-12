RNA-sequence-tools
==================

Tools for RNA seq and gene annotation

Tophat Cluster submission contains scripts for processing raw RNA-seq files for submission to a linux cluster (qb3 at UCSF specifically).  Additionally it contains some scripts for testing command formatting and managing files on a cluster.

FPKM parsing contains scripts for turning multiple tophat/cufflinks output files into FPKM matrices suitable for analysis in Fluidigm Singular R package.

RNA Seq analysis contains scripts for clustering and pca analysis of RNA-seq data

Gene Ontology contains files for fetching and organizing Entrez Gene ontology information from lists of genes

Sample Work Flow:                
1) Use tophat_qsub.py to submit sequencing to cluster -> Output: tophat and cufflinks fpkm files         
2) Use cuffnorm_qsub to create sample sheet and normalize sequecing reads with cuffnorm -> Output: cuffnorm gene fpkm table                
3) Use align_report2.py to create and alignment report on mapped reads -> Output: alignment report for filtering                
4) Filter samples using fpkm_df_math.py to delete cells or genes not expressed at threshold -> Output: filtered fpkm usable in monocle (https://github.com/cole-trapnell-lab/monocle-release) or Singular (Fluidigm)                            
5)Use Pca or corr_search to analyze data -> Work in progress
