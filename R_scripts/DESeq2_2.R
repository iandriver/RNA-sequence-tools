run_deseq2 <- function(file_path, save_name){
library(DESeq2)
raw.data <- read.csv(file = file_path, header=TRUE, sep="\t", row.names=1)
counts <- raw.data[ , -c(1,ncol(raw.data)) ]
rownames(counts) <- row.names(raw.data)
group <- c()
for (x in colnames( counts )){if(grepl('B_',x)){ group <- append(group, 'Basal')} else if(grepl('L_',x)){ group <- append(group,'Luminal')} }
print(colnames(counts))
print(group)
samples <- data.frame(row.names=colnames(counts), condition=as.factor(c(group)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
matrix <- counts(dds,normalized=TRUE)
write.table(matrix, file=paste("DESeq_", save_name,"matrix_norm.txt", sep = "_"), sep="\t")
saveRDS(dds, paste(save_name,".rds"))
plotDispEsts(dds)
return(dds)}
