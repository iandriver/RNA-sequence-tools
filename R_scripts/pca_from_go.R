pca_from_go <- function(exp){
	exp_1 <- updateGeneListFromFile(exp)
	pca_1 <- PCA(exp_1)
	pca_topgene_1 <- getTopPCAGenes (pca_1, top_gene_num=150)
	write.table(pca_topgene_1, file="pca_goterms_allPCA_list.txt", sep="\t", col.names=F)
	exp_2 <- updateGeneListFromList(exp_1, pca_topgene_1)
	hc_2 <- HC(exp_2)
	pca_2 <- PCA(exp_2)
	hc_sample_list_2 <- identifySampleClusters (hc_2)
	exp_3 <- updateSampleListFromList(exp_2, hc_sample_list_2)
	saveRDS(exp_3, file="exp_fso_gene_sample_clusters.rds")
	pca_3 <- PCA(exp_3)
	hc_3 <- HC(exp_3)
	pc_3_list <- displayPCALoading(pca_3)
	write.table(pc_3_list, file="pca_goterms_topPCALoading_selected.txt", sep="\t", col.names=F)
	exp_4 <- resetGeneList(exp_3)
	saveRDS(exp_4, file="exp_fso_sample_clusters_resetgenes.rds")
	anova_go <- ANOVA(exp_4)
	saveRDS(anova_go, file="anova_goterm_reset.rds")
	go_levels = levels(anova_go$"sample_list"$"GroupID")
	perms <- permutations(length(go_levels),2,v=l_levels)
	for (p in 1:nrow(perms)){
		term1 = perms[p,1]
		term2 = perms[p,2]
		fc_list <- foldChangeAnalysis(anova_go, term1, term2, foldchange_threshold = 1, pvalue_threshold = 1, display_plot = FALSE)
		write.table(fc_list, file=paste("Anova_", term1, term2,"folldchange_all.txt", sep = "_"), sep="\t")
	}
}
