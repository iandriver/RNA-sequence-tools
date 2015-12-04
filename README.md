RNA-sequence-tools
==================

Tools for RNA seq and gene annotation

Tophat Cluster submission contains scripts for processing raw RNA-seq files for submission to a linux cluster (qb3 at UCSF specifically).  Additionally it contains some scripts for testing command formatting and managing files on a cluster.

FPKM parsing contains scripts for turning multiple tophat/cufflinks output files into FPKM matrices suitable for analysis in Fluidigm Singular R package.

RNA Seq analysis contains scripts for clustering and pca analysis of RNA-seq data

Gene Ontology contains files for fetching and organizing Entrez Gene ontology information from lists of genes

Sample FPKM Work Flow:                
1) Use tophat_qsub.py to submit sequencing to cluster -> Output: tophat and cufflinks fpkm files         
2) Use cuffnorm_qsub to create sample sheet and normalize sequecing reads with cuffnorm -> Output: cuffnorm gene fpkm table               

Sample Count based Work Flow:         

1) Use tophat_qsub.py to submit sequencing to cluster -> Output: tophat and cufflinks fpkm files 
2) Use sort_htseq_count.py to clean up accepted hits and generate htseq counts and picard metric information (3' to 5' bias, CG, etc.)
3) Use R-scripts: DESeq or edgeR to process raw counts.

Data Analysis Tools:         

1) Use filter_outliers to filter on mapping rates,genes expressed, or other rule based metrics. -> outlier_filtered matrix           
2) Run cluster.py to do a broad unbiased clustering and subclustering search. Produces hierarchical clustering and pca and correlation groups for all cells and cell subgroups (down to a defined threshold for minimal number of cells in a subgroup)                          
3)Use cluster1.py for a more targeted search uses selected cells or gene files to establish starting point and does targeted significance searching based a single gene.
4)corr_search.py allows searching for correlation with any one gene for a whole matrix.
5)make_monocle_data.py only works with cuffnorm output (per monocles requirements: http://cole-trapnell-lab.github.io/monocle-release/) and produces the 3 files necassary for running monocle (gene and cell feature sheets). It requires gene_lookup.py, which fetches GO terms to populate the gene feature sheet. Also required is a table with fluidigm capture data (single clean cell capture or not.) Terms can be added to classify populations by gene expression. 
6) run_monocle.R runs takes the make_monocle_data.py output and runs through a basic monocle workflow (making spanning tree and plotting genes in pseudotime) and then finds the clustergroups and plots representative groups in psuedotime. 
7)make_scLVM.R is a script to do the ERCC normalization and cell cycle analysis using the scLVM package (https://github.com/PMBio/scLVM).





