run_deseq2 <- function(file_path, save_name){
library(DESeq2)
raw.data <- read.csv(file = file_path, header=TRUE, sep="\t", row.names=1)
counts <- raw.data[ , -c(1,ncol(raw.data)) ]
rownames(counts) <- row.names(raw.data)
group <- c()
for (x in colnames( counts )){if(grepl('BU3',x)){ group <- append(group, 'BU3')} else if(grepl('ips17',x)){ group <- append(group,'ips17')} }
print(colnames(counts))
print(group)
samples <- data.frame(row.names=colnames(counts), condition=as.factor(c(group)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
dds <- DESeq(dds)
saveRDS(dds, paste(save_name,".rds"))
plotDispEsts(dds)
return(dds)}
