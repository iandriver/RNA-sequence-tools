#R code for monocle2 ips17 BU3 single cell data Hawkins, et al.
#also refer to the monocle2 vignette http://www.bioconductor.org/packages/release/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf

#read in gene cell matrix file
fpkm_matrix <- read.delim("ips17_BU3_normalized_edgeR_counts_TMM_all_outlier_filtered_nojunk_cutoff.txt")

#assign gene symbol names as row.names
row.names(fpkm_matrix) <- fpkm_matrix$X

#rename the matrix as final
final_df <- fpkm_matrix
#filter out any blanks
final_df <- final_df[,!(names(final_df) %in% c('X'))]
fpkm_matrix <- final_df

#read in the Cell names and group assignments
sample_sheet <- read.delim("Cell_classifications_new.txt")
#make the rownames the SampleIDs
row.names(sample_sheet) <- sample_sheet$SampleID
#create the pheno-data DataFrame for monocle2
pd <- new("AnnotatedDataFrame", data = sample_sheet)

#match the SampleIDs in the sample sheet to the matrix
fpkm2 <- fpkm_matrix[,match(row.names(sample_sheet), colnames(fpkm_matrix))]


#create the (gene) feature sheet by taking all of the genes in the matrix and making a dataframe (no extra data needed)
feature_sheet <- as.data.frame(row.names(fpkm2))
row.names(feature_sheet) <- feature_sheet$"row.names(fpkm2)"
feature_sheet$GeneID <- row.names(feature_sheet)
fd <- new("AnnotatedDataFrame", data = feature_sheet)

#create the monocle2 CellDataSet object
ips17_BU3_data <- newCellDataSet(as.matrix(fpkm2), phenoData = pd, featureData=fd)

# As per monocle2 vignette, use it to estimate RNA counts
rpc_matrix <- relative2abs(ips17_BU3_data)
# Now, make a new CellDataSet using the RNA counts
ips17_BU3_data <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
phenoData = pd,
featureData = fd,
lowerDetectionLimit=1,
expressionFamily=negbinomial.size())

#pre-calculate size factors and dispersions
ips17_BU3_data <- estimateSizeFactors(ips17_BU3_data)
ips17_BU3_data <- estimateDispersions(ips17_BU3_data)

#calulate the number of cells that each gene is expressed in
ips17_BU3_data <- detectGenes(ips17_BU3_data , min_expr = 0.1)
#use that calculation to create expressed_genes
expressed_genes <- row.names(subset(fData(ips17_BU3_data), num_cells_expressed >= 6))

#filter cells with way to few or way to many genes expressed
pData(ips17_BU3_data)$Total_mRNAs <- Matrix::colSums(exprs(ips17_BU3_data))
upper_bound <- 10^(mean(log10(pData(ips17_BU3_data)$Total_mRNAs)) + 2*sd(log10(pData(ips17_BU3_data)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(ips17_BU3_data)$Total_mRNAs)) - 2*sd(log10(pData(ips17_BU3_data)$Total_mRNAs)))

#plot mRNA expression by Nkx_group or other
qplot(Total_mRNAs, data=pData(ips17_BU3_data), color=Nkx_group, geom="density") + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=upper_bound)

ips17_BU3_data <- ips17_BU3_data[,pData(ips17_BU3_data)$Total_mRNAs > lower_bound & pData(ips17_BU3_data)$Total_mRNAs < upper_bound]
ips17_BU3_data <- detectGenes(ips17_BU3_data, min_expr = 0.1)

#create gene ids for classification
NKX_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "NKX2-1"))
APOA2_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "APOA2"))
CD47_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "CD47"))
MSX1_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "MSX1"))
SFTA3_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "SFTA3"))
COL19A1_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "COL19A1"))
#create empty newCellTypeHierarchy classifier
cth <- newCellTypeHierarchy()

#define "Lung" classifier expression rules
#expresses NKX2-1 and SFTA3 or
#expresses CD47 and SFTA3 or
#expresses NKX2-1 and CD47

cth <- addCellType(cth, "Lung", classify_func=function(x) {
              (x[NKX_id,] >= 1 & x[SFTA3_id,] >=1) |
              (x[CD47_id,] >= 1 & x[SFTA3_id,]>=1) |
              (x[NKX_id,] >= 1 & x[CD47_id,] >=1)
              })

#define "Liver" classifier expression rules
#doesnt express NKX2-1 and does express APOA2 or
#doesnt express NKX2-1 CD47 and does express MSX1 or does express COL19A1 or
#doesnt express NKX2-1 or CD47
cth <- addCellType(cth, "Liver", classify_func=function(x) {
              (x[NKX_id,] < 1 & x[APOA2_id,] > 1) |
              (x[NKX_id,] < 1 & x[MSX1_id,] > 1 | x[COL19A1_id,]>1) |
              (x[NKX_id,] < 1 & x[CD47_id,] < 1)
              })

ips17_BU3_data <- classifyCells(ips17_BU3_data, cth, 0.1)

pie <- ggplot(pData(ips17_BU3_data), aes(x = factor(1), fill = factor(CellType))) +
geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x=element_blank(), axis.title.y=element_blank())

fData(ips17_BU3_data)$GeneSum <- Matrix::rowSums(exprs(ips17_BU3_data))

expressed_genes <- row.names(subset(fData(ips17_BU3_data), num_cells_expressed >= 6 & GeneSum >=80))

diff_test_res <- differentialGeneTest(ips17_BU3_data[expressed_genes,], fullModelFormulaStr="~CellType", cores=3)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cell_cycle_genes <- read.delim('cell_cycle_genes.txt')

usup_filtered3 <- setdiff(ordering_genes, cell_cycle_genes$GeneID)

best_gene_list <- read.delim('kmeans_2_Best_Gene_list_from_scicast.txt')

best_genes <- best_gene_list$GeneID

best_ordering <- as.factor(c(as.character(best_genes), as.character(usup_filtered3)))

best_ord <- best_ordering[!duplicated(best_ordering)]

ips17_BU3_data <- setOrderingFilter(ips17_BU3_data, best_ord)
ips17_BU3_data <- reduceDimension(ips17_BU3_data, max_components=2)
ips17_BU3_data <- orderCells(ips17_BU3_data, reverse=TRUE)
