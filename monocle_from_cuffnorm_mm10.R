make_monocle_cuffnorm <- function(fpkm_path, sample_sheet_path){
  library(ggplot2)
  library(monocle)
  #paths to files provided
  #fpkm_path = "hisat2_cuffnorm_gene_fpkm_PR8_x31_Saline_all.txt"
  #sample_sheet_path = "hisat2_cuffnorm_all_cellcroups.txt"
  fpkm_matrix <- read.table(fpkm_path, check.names = F, header = TRUE)
  row.names(fpkm_matrix) <- fpkm_matrix$tracking_id
  fpkm_matrix <- fpkm_matrix[,2:length(colnames(fpkm_matrix))]

  sample_sheet = read.delim(sample_sheet_path)
  row.names(sample_sheet) <- sample_sheet$SampleID
  fpkm2 <- fpkm_matrix[,match(row.names(sample_sheet), colnames(fpkm_matrix))]
  fpkm2 <- fpkm2[,colSums(fpkm2) > 100]
  fpkm2 <- fpkm2[rowSums(fpkm2) > 10,]

  #filters out Gm and Rik and other unannotated genes
  fpkm2$GeneID <- rownames(fpkm2)
  fpkm_filter <- filter(fpkm2, !grepl("^Gm[0-9]{2,5}$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("Rik$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("Rik[0-9]$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("[0-9]{6}$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("[0-9]{6}.[0-9]$", GeneID))
  rownames(fpkm_filter) <- fpkm_filter$GeneID
  fpkm2 <- fpkm_filter[,1:length(colnames(fpkm_filter))-1]
  write.table(fpkm2, file = "fpkm_filtered_to_remove_gm_rik_genes.txt", sep ="\t", col.names = NA)

  #make feature sheet from genes and make monocle cell data set
  feature_sheet <- as.data.frame(row.names(fpkm2))
  row.names(feature_sheet) <- feature_sheet$"row.names(fpkm2)"
  feature_sheet$GeneID <- row.names(feature_sheet)
  fd <- new("AnnotatedDataFrame", data = feature_sheet)
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  HSMM <- newCellDataSet(as.matrix(fpkm2), phenoData = pd, featureData=fd, expressionFamily = negbinomial())
  rpc_matrix <- relative2abs(HSMM)
  # Now, make a new CellDataSet using the RNA counts
  HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit=1,
  expressionFamily=negbinomial.size())

  #pre-calculate size factors and dispersions
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)

  #calulate the number of cells that each gene is expressed in
  HSMM <- detectGenes(HSMM , min_expr = 0.1)
  #create new feature data for the total sum of gene expression across all cells
  fData(HSMM)$GeneSum <- Matrix::rowSums(exprs(HSMM))

  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 6))
  pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
  upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))

  qplot(Total_mRNAs, data=pData(HSMM), color=GroupID, geom="density") + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=upper_bound)
  ggsave("mRNA_density_plot_byGroupID_prefilter.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 9, units = c("in"), dpi = 300)
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & pData(HSMM)$Total_mRNAs < upper_bound]
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  qplot(Total_mRNAs, data=pData(HSMM), color=GroupID, geom="density")
  ggsave("mRNA_density_plot_byGroupID_after_filter.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 9, units = c("in"), dpi = 300)
  relative_expr_matrix <- exprs(HSMM)
  relative_expr_matrix <- apply(relative_expr_matrix, 2, function(x) x / sum(x) * 10^6)
  write.table(relative_expr_matrix, file = paste("cpm_censusnorm_monocle_",fpkm_path), sep ="\t", col.names = NA)

  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 6 & GeneSum >3))
  diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr="~GroupID", cores=3)
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  best_gene_df <- subset(fData(HSMM[ordering_genes,]))

  #save significant genes
  write.table(best_gene_df[order(-best_gene_df[,"GeneSum"]),], file="significant_genes_by_GroupID_ordered_by_expression.txt", sep="\t", col.names=NA)
  #create dataframe and order by highest expression
  best_gene_df <- best_gene_df[order(-best_gene_df[,"GeneSum"]),]
  diff_test_df <- diff_test_res[,c("GeneID", "pval", "qval")]
  df_qval <- diff_test_df[match(row.names(best_gene_df), row.names(diff_test_df)),]
  df_qval$GeneSum <- best_gene_df$GeneSum
  df_qval$num_cells_expressed <- best_gene_df$num_cells_expressed
  top_30_expression_to_plot <- row.names(df_qval[1:min(dim(df_qval)[1],30),])

  #plot top 30 significant genes by expression
  jit <- plot_genes_jitter(HSMM[top_30_expression_to_plot,], grouping="GroupID", color_by="GroupID", nrow=6, ncol=5, plot_trend = TRUE)
  jit + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("top30_byexpression_sig.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 14, height = 11, units = c("in"), dpi = 300)

  #sort to rank by qval
  df_qval <- df_qval[order(df_qval[,"qval"]),]
  top_30_qval_to_plot <- row.names(df_qval[1:min(dim(df_qval)[1],30),])
  #plot top genes by qval
  jit <- plot_genes_jitter(HSMM[top_30_qval_to_plot,], grouping="GroupID", color_by="GroupID", nrow=6, ncol=5, plot_trend = TRUE)
  jit + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("top30_byqval_sig.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 14, height = 11, units = c("in"), dpi = 300)

  return(HSMM)

}
