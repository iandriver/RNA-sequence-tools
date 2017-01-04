run_deseq2 <- function(file_path, save_name, group_term_list, group_name_list){
library(DESeq2)
raw.data <- read.csv(file = file_path, header=TRUE, sep="\t", row.names=1)
counts <- raw.data
rownames(counts) <- row.names(raw.data)
group <- c()

group_map <- mapply(c, group_term_list, group_name_list, SIMPLIFY=FALSE)

for (x in colnames( counts )){
  for (group_set in group_map){
    if(grepl(group_set[1],x)){ group <- append(group, group_set[2])} } }
print(length(colnames(counts)))
print(length(group))
samples <- data.frame(row.names=colnames(counts), condition=as.factor(c(group)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- nbinomWaldTest(dds)
matrix <- counts(dds,normalized=TRUE)
vsd <- varianceStabilizingTransformation(dds)
vsd_matrix <- assay(vsd)
write.table(matrix, file=paste("DESeq_", save_name,"matrix_norm.txt", sep = "_"), sep="\t", col.names=NA)
write.table(vsd_matrix, file=paste("DESeq2_", save_name,"matrix_norm_vsd.txt", sep = "_"), sep="\t", col.names=NA)
saveRDS(dds, paste(save_name,".rds", sep = ""))
plotDispEsts(dds)
return(dds)}
