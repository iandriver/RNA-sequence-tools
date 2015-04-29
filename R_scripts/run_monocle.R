run_monocle <- function(dir){
	library(monocle)
	library(reshape2)
	library(plyr)
	setwd(dir)
	fpkm_matrix <- read.delim("genes.fpkm_table", row.names=1)
	sample_sheet <- read.delim("samples.table", row.names=1)
	feature_sheet <- read.delim("genes.attr_table", row.names=1)
	colnames(fpkm_matix) <- row.names(sample_sheet)
	pd <- new("AnnotatedDataFrame", data = sample_sheet)
	fd <- new("AnnotatedDataFrame", data = feature_sheet)
	my_data <- newCellDataSet(as.matrix(fpkm_matrix), phenoData = pd, featureData=fd)
	my_data <- detectGenes(my_data, min_expr = 1)
	expressed_genes <- row.names(subset(fData(my_data), num_cells_expressed >= 50))
	valid_cells <- row.names(subset(pData(my_data), total_mass>=100000)
	my_data <- my_data[,valid_cells]
	L <- log(exprs(my_data[expressed_genes,]))
	melted_dens_df <- melt(t(scale(t(L))))
	qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab("Standardized log(FPKM)") + ylab("Density")
	marker_genes <- row.names(subset(fData(my_data), gene_short_name %in% c('Scnn1a', 'Pdpn', 'Hopx', 'Aqp5', 'Ager')))
	diff_test_res_cond <- differentialGeneTest(my_data[marker_genes,], fullModelFormulaStr="expression~condition")
	diff_test_res_day <- differentialGeneTest(my_data[marker_genes,], fullModelFormulaStr="expression~day")
	sig_genes_cond <- subset(diff_test_res_cond, qval < 0.1)
	sig_genes_day <- subset(diff_test_res_day, qval < 0.1)
	sig_genes_cond <- merge(fData(my_data), sig_genes_cond, by="row.names")
	sig_genes_day <- merge(fData(my_data), sig_genes_day, by="row.names")
	sig_genes_cond[,c("gene_short_name", "pval", "qval")]
	sig_genes_day[,c("gene_short_name", "pval", "qval")]
	type1 <- my_data[row.names(subset(fData(my_data),gene_short_name %in% c("Ager", "Aqp5","Pdpn"))),]
	plot_genes_jitter(type1, grouping="condition", ncol=3)
	}