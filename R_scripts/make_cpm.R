library(edgeR)
make_cpm <- function(path_to_file,file_base_name){
	raw.data <- read.csv(file = file_path, header=TRUE, sep="\t", row.names=1)
	counts <- raw.data[ , -c(1,ncol(raw.data)) ]
	rownames( counts ) <- row.names(raw.data)
	cds <- DGEList( counts , group = colnames(counts) )
	cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
	cds <- calcNormFactors(cds, method="TMM")
	cps <- cpm(cds, normalized.lib.sizes=TRUE)
	write.table(cps, file=paste(file_base_name,"normalized_cpm_all.txt",sep='_'), sep="\t", col.names=NA)
}
