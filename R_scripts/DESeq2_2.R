run_deseq2 <- function(file_path, save_name){
library(DESeq2)
raw.data <- read.table(file = file_path, header=TRUE)
counts <- raw.data[ , -c(1,ncol(raw.data)) ]
rownames(counts) <- row.names(raw.data)
group <- c()
for (x in colnames( counts )){if(grepl('ctrl',x)){ group <- append(group, 'ctrl')} else if(grepl('pnx',x)){ group <- append(group,'pnx')} }
print(colnames(counts))
print(group)
samples <- data.frame(row.names=colnames(counts), condition=as.factor(c(group)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
dds <- DESeq(dds)
saveRDS(dds, paste(save_name,".rds"))
plotDispEsts(dds)
return(dds)}