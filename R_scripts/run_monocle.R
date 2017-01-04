run_monocle <- function(directory){
	library(monocle)
	library(reshape2)
	library(plyr)
	subdir = paste(directory,'/','monocle_2', sep='')
	dir.create(subdir)
	setwd(subdir)
	fpkm_matrix <- read.delim(paste(directory,"/goterms_monocle_count_matrix.txt", sep=''), row.names=1)
	sample_sheet <- read.delim(paste(directory,"/cell_feature_data.txt", sep=''))
	feature_sheet <- read.delim(paste(directory,"/gene_feature_data.txt", sep=''))
	row.names(sample_sheet) <- sample_sheet$tracking_id
	row.names(feature_sheet) <- feature_sheet$GeneID
  fpkm2 <- fpkm_matrix[,match(row.names(sample_sheet), colnames(fpkm_matrix))]
	pd <- new("AnnotatedDataFrame", data = sample_sheet)
	fd <- new("AnnotatedDataFrame", data = feature_sheet)
	my_data <- newCellDataSet(as.matrix(fpkm2), phenoData = pd, featureData=fd)
	my_data <- detectGenes(my_data, min_expr = 1)
	expressed_genes <- row.names(subset(fData(my_data), num_cells_expressed >= 5))

	valid_cells <- row.names(subset(pData(my_data), total_mass>=100000 & single_cell=='yes' & per_mapped >=50))
	my_data <- my_data[,valid_cells]
	saveRDS(my_data, "monocle_data_pdgfra1.rds")
	L <- log(exprs(my_data[expressed_genes,]))
	melted_dens_df <- melt(t(scale(t(L))))
	pdf("qplot_density.pdf")
	qp <- qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab("Standardized log(FPKM)") + ylab("Density")
	print(qp)
	dev.off()
	marker_genes <- row.names(subset(fData(my_data), GeneID %in% c('SFTPB', 'SFTPC', 'NKX2-1', 'KRT5', 'KRT7','ABCA3', 'COL1A1', 'MYOF', 'NOTCH3', 'ITGA2', 'HES1', 'ITGAV', 'PDGFRA', 'DUOX1', 'HOPX', 'DPP4', 'BMP2', 'FOXA1', 'FOXJ1', 'KRT15')))
	pdgfra_marker_genes <- row.names(subset(fData(my_data), GeneID %in% c('Dcn', 'Serpinf1', 'Lmna', 'Serping1', 'Cxcl1', 'Lum', 'Sod1', 'Apoe', 'Thbs2', 'Has1', 'Has2', 'Notch2', 'Col5a3', 'Dpt', 'Tcf21', 'Acta2', 'Ly6e', 'Ly6a', 'Pdgfra', 'Jun')))
	notch_genes<- row.names(subset(fData(my_data), GeneID %in% c('Sftpc', 'Notch3', 'Hes1', 'Hey1', 'Hey2', 'Hes4', 'Hes6', 'Notch1', 'Notch2', 'Notch4', 'Dll1', 'Jag1', 'Jag2', 'Adam17', 'Numb', 'Psen1', 'Psen2', 'Psenen')))
	diff_test_res_cond <- differentialGeneTest(my_data[pdgfra_marker_genes,], fullModelFormulaStr="expression~hi_low")
	sig_genes_cond <- subset(diff_test_res_cond, qval < 0.1)
	print(sig_genes_cond)

	sig_genes_cond <- merge(fData(my_data), sig_genes_cond, by="row.names")

	sig_genes_cond[,c("GeneID", "pval", "qval")]

	diff_test_res <- differentialGeneTest(my_data[expressed_genes,], fullModelFormulaStr="expression~hi_low")
	ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
	ordering_genes <- intersect(ordering_genes, expressed_genes)
	my_data <- setOrderingFilter(my_data, ordering_genes)
	my_data <- reduceDimension(my_data, use_irlba=FALSE)
	my_data <- orderCells(my_data, num_paths=3, reverse=FALSE)
	pdf('spanning_tree_plot.pdf', width=12, height=12)
	span_tree <- plot_spanning_tree(my_data, x = 1, y = 2, color_by = "State", show_tree = TRUE, show_backbone = TRUE, backbone_color = "black", markers = NULL, show_cell_names = TRUE, cell_name_size = 2)
	print(span_tree)
	dev.off()
	pdf('marker_genes_psuedotime.pdf', width=14, height=9)
	pstime <- plot_genes_in_pseudotime(my_data[pdgfra_marker_genes,], min_expr = NULL, cell_size = 0.75, nrow = NULL, ncol = 4, panel_order = NULL, color_by = "State", trend_formula = "adjusted_expression ~ sm.ns(Pseudotime, df=3)", label_by_short_name = TRUE)
	print(pstime)
	dev.off()
	saveRDS(my_data, "monocle_data_pdgfra2.rds")
	#find cluster groups and plot out a subset of each group in pseudotime
	full_model_fits <- fitModel(my_data, modelFormulaStr="expression~sm.ns(Pseudotime, df=3)")
	expression_curve_matrix <- responseMatrix(full_model_fits)
	clusters <- clusterGenes(expression_curve_matrix, k=4)
	group1 <- c()
	group2 <- c()
	group3 <- c()
	group4 <- c()
	for(p in names(clusters$clustering)){
		for(group in clusters$clustering[p]){
			if(group == 1){
				group1 <- append(group1, p)
			}
			if(group == 2){
				group2 <- append(group2, p)
			}
			if(group == 3){
				group3 <- append(group3, p)
			}
			if(group == 4){
				group4 <- append(group4, p)

			}
		}
	}
	write.table(group1, file="gene_cluster_1.txt", sep="\t", col.names='group1_genes')
	write.table(group2, file="gene_cluster_2.txt", sep="\t", col.names='group2_genes')
	write.table(group3, file="gene_cluster_3.txt", sep="\t", col.names='group3_genes')
	write.table(group4, file="gene_cluster_4.txt", sep="\t", col.names='group4_genes')
	pdf(file='group1_cluster_genes_in_pseudotime.pdf', width=14, height=9)
	group1_genes <- row.names(subset(fData(my_data), GeneID %in% group1[1:24]))
	pp1 <- plot_genes_in_pseudotime(my_data[group1_genes,], min_expr = NULL, cell_size = 0.75, nrow = NULL, ncol = 4, panel_order = NULL, color_by = "State", trend_formula = "adjusted_expression ~ sm.ns(Pseudotime, df=3)", label_by_short_name = TRUE)
	print(pp1)
	dev.off()
	pdf(file='group2_cluster_genes_in_pseudotime.pdf', width=14, height=9)
	group2_genes <- row.names(subset(fData(my_data), GeneID %in% group2[1:24]))
	pp2 <- plot_genes_in_pseudotime(my_data[group2_genes,], min_expr = NULL, cell_size = 0.75, nrow = NULL, ncol = 4, panel_order = NULL, color_by = "State", trend_formula = "adjusted_expression ~ sm.ns(Pseudotime, df=3)", label_by_short_name = TRUE)
	print(pp2)
	dev.off()
	pdf(file='group3_cluster_genes_in_pseudotime.pdf', width=14, height=9)
	group3_genes <- row.names(subset(fData(my_data), GeneID %in% group3[1:24]))
	pp3 <- plot_genes_in_pseudotime(my_data[group3_genes,], min_expr = NULL, cell_size = 0.75, nrow = NULL, ncol = 4, panel_order = NULL, color_by = "State", trend_formula = "adjusted_expression ~ sm.ns(Pseudotime, df=3)", label_by_short_name = TRUE)
	print(pp3)
	dev.off()
	pdf(file='group4_cluster_genes_in_pseudotime.pdf', width=14, height=9)
	group4_genes <- row.names(subset(fData(my_data), GeneID %in% group4[1:24]))
	pp4 <- plot_genes_in_pseudotime(my_data[group4_genes,], min_expr = NULL, cell_size = 0.75, nrow = NULL, ncol = 4, panel_order = NULL, color_by = "State", trend_formula = "adjusted_expression ~ sm.ns(Pseudotime, df=3)", label_by_short_name = TRUE)
	print(pp4)
	dev.off()
	}
