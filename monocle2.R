make_monocle <- function(fpkm_path, sample_sheet_path, group_name=CellGroup){
  source("/Users/idriver/RockLab-files/RNA-sequence-tools/R_scripts/rename_ensembl_and_remove_dups.R")

  fpkm_matrix <- read.table(fpkm_path, check.names = F, header = TRUE)
  row.names(fpkm_matrix) <- fpkm_matrix$gene_id
  fpkm_matrix <- fpkm_matrix[,2:length(colnames(fpkm_matrix))]
  fpkm_matrix <- rename_and_remove_dups(fpkm_matrix, keep_highest=TRUE, species="mouse")
  sample_sheet <- read.delim(sample_sheet_path)
  row.names(sample_sheet) <- sample_sheet$SampleID
  sample_sheet$CellGroup <- sample_sheet$CellClusters
  fpkm2 <- fpkm_matrix[,match(row.names(sample_sheet), colnames(fpkm_matrix))]
  feature_sheet <- as.data.frame(row.names(fpkm2))
  head(feature_sheet)
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
  #HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & pData(HSMM)$Total_mRNAs < upper_bound]
  #HSMM <- detectGenes(HSMM, min_expr = 0.1)
  qplot(Total_mRNAs, data=pData(HSMM), color=CellGroup, geom="density") + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=upper_bound)
  relative_expr_matrix <- exprs(HSMM)
  relative_expr_matrix <- apply(relative_expr_matrix, 2, function(x) x / sum(x) * 10^6)
  write.table(relative_expr_matrix, file = "cpm_cfms_pnxd7_monocle_norm_filter_from_count.txt", sep ="\t", col.names = NA)
  return(HSMM)

}
