pca_from_go <- function(exp){
	exp_spc2_n3 <- updateGeneListFromFile(exp)
	pca_n3 <- PCA(exp_spc2_n3)
	pca_topgene_list3 <- getTopPCAGenes (pca_n3, top_gene_num=150)
	exp_spc2_n4 <- updateGeneListFromList(exp_spc2_n3, pca_topgene_list3)
	hc_n4 <- HC(exp_spc2_n4)
	pca_n4 <- PCA(exp_spc2_n4)
	hc_sample_list_n4 <- identifySampleClusters (hc_n4)
	exp_spc2_n5 <- updateSampleListFromList(exp_spc2_n4, hc_sample_list_n4)
	pca_n5 <- PCA(exp_spc2_n5)
	hc_n5 <- HC(exp_spc2_n5)
	pc_n5_list <- displayPCALoading(pca_n5)
}