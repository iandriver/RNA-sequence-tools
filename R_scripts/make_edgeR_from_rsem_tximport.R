
make_cpm <- function(sample_file_path, save_name, group_term_list, group_name_list, convert_ensembl, species){
	library(edgeR)
	library(tximport)
	library(readr)

	samples <- read.csv(file = sample_file_path, header=TRUE, sep="\t", row.names=1)
	files <- file.path(samples$folder, paste0(samples$run, ".genes.results"))
	names(files) <- samples$run
	txi.rsem <- tximport(files, type = "rsem", reader = read_tsv, countsFromAbundance="lengthScaledTPM")
	saveRDS(txi.rsem,file="txi_rsem.RDS")
	cts1 <- txi.rsem$counts
	print('original count diminsion:')
	print(dim(cts1))
	normMat1 <- txi.rsem$length
	cts2 <- cts1[,colSums(cts1)>100000]
	cts3 <- cts2[rowSums(cts2)>10,]
	print('cells removed for low expression:')
	print(colnames(cts1[,colSums(cts1)<100000]))
	print("diminsions after column filtering:")
	print(dim(cts2))
	cts<-cts3[apply(cts3[,-1], 1, function(x) !all(x==0)),]
	print("diminsions after row filtering:")
	print(dim(cts))
	normMat2 <- normMat1[,colSums(cts1)>100000]
	normMat3 <- normMat2[rowSums(cts2)>10,]
	normMat <- normMat3[apply(cts3[,-1], 1, function(x) !all(x==0)),]
	ratio <- cts/normMat
	row.has.na <- apply(ratio, 1, function(x){any(is.na(x))})
	ratio.filtered <- ratio[!row.has.na,]
	cts.filtered <- cts[!row.has.na,]
	print("diminsions after na filtering:")
	print(dim(cts.filtered))
	normMat.filtered <- normMat[!row.has.na,]
	#o <- log(calcNormFactors(ratio.filtered)) + log(colSums(ratio.filtered))

	group <- c()

	group_map <- mapply(c, group_term_list, group_name_list, SIMPLIFY=FALSE)

	for (x in colnames(cts.filtered)){
	  for (group_set in group_map){
	    if(grepl(group_set[1],x)){ group <- append(group, group_set[2])} } }
	print(length(group))
	group <- factor(group)
	y <- DGEList(counts = cts.filtered, group=group)
	#y$offset <- t(t(log(normMat.filtered)) + o)
	y <- calcNormFactors(y)
	design <- model.matrix(~group)
	y <- estimateDisp(y,design)
	fit <- glmFit(y,design)
	saveRDS(y, file=paste(save_name,"edgeR_DGE_object_y_final.rds",sep='_'))
	saveRDS(fit, file=paste(save_name,"edgeR_glm_fit.rds",sep='_'))
	saveRDS(group, file=paste(save_name,"edgeR_group_factor.rds",sep='_'))
	cds <- calcNormFactors(y, method="TMM")
	cps <- cpm(cds, normalized.lib.sizes=TRUE)

	if(convert_ensembl==TRUE){
		source("/Users/idriver/RockLab-files/RNA-sequence-tools/R_scripts/rename_ensembl_and_remove_dups.R")
		renamed_cps <- rename_and_remove_dups(df=cps, TRUE, species)
		write.table(renamed_cps, file=paste(save_name,"normalized_edgeR_counts_TMM_all.txt",sep='_'), sep="\t", col.names=NA)
		}
	else{
		write.table(cps, file=paste(save_name,"normalized_edgeR_counts_TMM_all.txt",sep='_'), sep="\t", col.names=NA)
	}
}
