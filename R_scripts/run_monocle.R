run_monocle <- function(dir){
	library(monocle)
	library(reshape2)
	library(plyr)
	setwd(dir)
	fpkm_matrix <- read.delim("goterms_monocle_count_matrix.txt", row.names=1)
	sample_sheet <- read.delim("cell_feature_data.txt", row.names=1)
	feature_sheet <- read.delim("gene_feature_data.txt", row.names=1)
<<<<<<< Updated upstream
=======
	head(fpkm_matrix)
	head(sample_sheet)
	head(feature_sheet)
>>>>>>> Stashed changes
	colnames(fpkm_matrix) <- row.names(sample_sheet)
	row.names(fpkm_matrix) <- row.names(feature_sheet)
	pd <- new("AnnotatedDataFrame", data = sample_sheet)
	fd <- new("AnnotatedDataFrame", data = feature_sheet)
	my_data <- newCellDataSet(as.matrix(fpkm_matrix), phenoData = pd, featureData=fd)
	my_data <- detectGenes(my_data, min_expr = 1)
	expressed_genes <- row.names(subset(fData(my_data), num_cells_expressed >= 50))
<<<<<<< Updated upstream
	valid_cells <- row.names(subset(pData(my_data), total_mass>=100000 & single_cell == 'yes' ))
	my_data <- my_data[,valid_cells]
	L <- log(exprs(my_data[expressed_genes,]))
	melted_dens_df <- melt(t(scale(t(L))))
	pdf("monocle.pdf")
	qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab("Standardized log(FPKM)") + ylab("Density")
	marker_genes <- row.names(subset(fData(my_data), GeneID %in% c('Scnn1a', 'Pdpn', 'Hopx', 'Aqp5', 'Ager')))
=======
	valid_cells <- row.names(subset(pData(my_data), total_mass>=100000 & single_cell=='yes' & per_mapped >=50))
	my_data <- my_data[,valid_cells]
	L <- log(exprs(my_data[expressed_genes,]))
	melted_dens_df <- melt(t(scale(t(L))))
	pdf("qplot_density.pdf")
	qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab("Standardized log(FPKM)") + ylab("Density")
	marker_genes <- row.names(subset(fData(my_data), tracking_id %in% c('Scnn1a', 'Pdpn', 'Hopx', 'Aqp5', 'Ager')))
>>>>>>> Stashed changes
	diff_test_res_cond <- differentialGeneTest(my_data[marker_genes,], fullModelFormulaStr="expression~condition")
	diff_test_res_day <- differentialGeneTest(my_data[marker_genes,], fullModelFormulaStr="expression~day")
	sig_genes_cond <- subset(diff_test_res_cond, qval < 0.1)
	sig_genes_day <- subset(diff_test_res_day, qval < 0.1)
	sig_genes_cond <- merge(fData(my_data), sig_genes_cond, by="row.names")
	sig_genes_day <- merge(fData(my_data), sig_genes_day, by="row.names")
<<<<<<< Updated upstream
	sig_genes_cond[,c("GeneID", "pval", "qval")]
	sig_genes_day[,c("GeneID", "pval", "qval")]
	type1 <- my_data[row.names(subset(fData(my_data),GeneID  %in% c("Ager", "Aqp5","Pdpn"))),]
	plot_genes_jitter(type1, grouping="condition", ncol=3)
	diff_test_res <- differentialGeneTest(my_data[expressed_genes,], fullModelFormulaStr="expression~condition")
	ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
	ordering_genes <- intersect(ordering_genes, expressed_genes)
	my_data <- setOrderingFilter(my_data, ordering_genes)
	my_data <- reduceDimension(my_data, use_irlba=FALSE)
	my_data <- orderCells(my_data, num_paths=2, reverse=TRUE)
	plot_spanning_tree(my_data)
	dev.off()
	saveRDS(my_data, "monocle_data_spc2_n2.rds")
=======
	sig_genes_cond[,c("gene_short_name", "pval", "qval")]
	sig_genes_day[,c("gene_short_name", "pval", "qval")]
	type1 <- my_data[row.names(subset(fData(my_data),tracking_id %in% c("Ager", "Aqp5","Pdpn"))),]
	plot_genes_jitter(type1, grouping="condition", ncol=3)
	dev.off()
>>>>>>> Stashed changes
	}
