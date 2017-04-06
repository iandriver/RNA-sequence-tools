make_monocle <- function(fpkm_path){
  library(monocle)
  library(dplyr)

  fpkm_matrix <- read.table(fpkm_path, check.names = F, header = TRUE, sep=',')
  row.names(fpkm_matrix) <- fpkm_matrix[,1]
  fpkm_matrix <- fpkm_matrix[,2:length(colnames(fpkm_matrix))]

  #Remove cells with very low counts (empty)
  fpkm2 <- fpkm_matrix[,colSums(fpkm_matrix) > 200]

  #Filter out Gm and Rik and other unannoted genes
  fpkm2$GeneID <- rownames(fpkm2)
  fpkm_filter <- filter(fpkm2, !grepl("^Gm[0-9]{2,5}$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("Rik$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("Rik[0-9]$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("[0-9]{6}$", GeneID))
  fpkm_filter <- filter(fpkm_filter, !grepl("[0-9]{6}.[0-9]$", GeneID))
  rownames(fpkm_filter) <- fpkm_filter$GeneID
  fpkm2 <- fpkm_filter[,1:length(colnames(fpkm_filter))-1]
  write.table(fpkm2, file = "unnormalized_filtered_to_remove_gm_rik_genes.txt", sep ="\t", col.names = NA)

  #Create basic groups from Origin of the cell
  group_term_list = c('Saline', 'x31', 'PR8')
  group_name_list = c('Saline', 'x31', 'PR8')
  group <- c()

  group_map <- mapply(c, group_term_list, group_name_list, SIMPLIFY=FALSE)

  for (x in colnames(fpkm2)){
    for (group_set in group_map){
      if(grepl(group_set[1],x)){ group <- append(group, group_set[2])} } }
  sample_sheet <- data.frame(row.names=colnames(fpkm2), CellGroup=as.factor(c(group)))

  #create a feature sheet with the gene names (no other info intitally)
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
  #HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & pData(HSMM)$Total_mRNAs < upper_bound]
  #HSMM <- detectGenes(HSMM, min_expr = 0.1)
  qplot(Total_mRNAs, data=pData(HSMM), color=CellGroup, geom="density") + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=upper_bound)
  relative_expr_matrix <- exprs(HSMM)
  relative_expr_matrix <- apply(relative_expr_matrix, 2, function(x) x / sum(x) * 10^6)
  write.table(relative_expr_matrix, file = "cpm_censusnorm_from_monocle.txt", sep ="\t", col.names = NA)
  
  return(HSMM)

}
