run_scLVM_norm <- function(filepath, base_name, species){
	library(scLVM)
	library(genefilter)
	library(statmod)
	require(ggplot2)
	library(gplots)
	require(DESeq2)
	library(EBImage)
	library(rhdf5)
	library(org.Mm.eg.db)
	library(hom.Hs.inp.db)
	data_frame_raw <- read.csv(file = filepath, header=TRUE, sep="\t", row.names=1)
	data_col_filtered <- data_frame_raw[,colSums(data_frame_raw) >= 10000]
	dat_row_col_filtered <- data_col_filtered[!!rowSums(abs(data_col_filtered[-c(1:2)])),]
	i_ERCC <- grep('^*ERCC-', substr( rownames(dat_row_col_filtered), 1, 5 ))
	countsERCC <- dat_row_col_filtered[i_ERCC,]
	lengthERCC <- dat_row_col_filtered[i_ERCC,1]
	i_Gene <- grep('^[^ERCC-]', substr( rownames(dat_row_col_filtered), 1, 5 ))
	countsGene <- dat_row_col_filtered[i_Gene,]
	lengthGene <- dat_row_col_filtered[i_Gene,1]
	
	sfERCC <- estimateSizeFactorsForMatrix(countsERCC)
	sfGene <- sfERCC
	
	nCountsERCC <- t( t(countsERCC) / sfERCC )
	nCountsGene <- t( t(countsGene) / sfGene )
	
	write.table(nCountsERCC, file=paste(base_name,"normalized_ERCC.txt",sep='_'), sep="\t", col.names=NA)
	write.table(nCountsGene, file=paste(base_name,"normalized_Genes_by_ERCC.txt",sep='_'), sep="\t", col.names=NA)
	
	meansGene <- rowMeans( nCountsGene )
	varsGene <- rowVars( nCountsGene )
	cv2Gene <- varsGene / meansGene^2
	
	meansERCC <- rowMeans( nCountsERCC )
	varsERCC <- rowVars( nCountsERCC )
	cv2ERCC <- varsERCC / meansERCC^2
	
	minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
	useForFitA <- meansERCC >= minMeanForFitA
	fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),cv2ERCC[useForFitA] )
	
	plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA)
	xg <- 10^seq( -3, 5, length.out=100 )
	lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
	segments( meansERCC[useForFitA], cv2ERCC[useForFitA], meansERCC[useForFitA], fitA$fitted.values, col="gray" )
	
	minBiolDisp <- .5^2
	xi <- mean( 1 / sfERCC )
	m <- ncol(countsGene)
	psia1thetaA <- mean( 1 / sfERCC ) + ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfGene )
	cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
	testDenomA <- ( meansGene * psia1thetaA + meansGene^2 * cv2thA ) / ( 1 + cv2thA/m )
	pA <- 1 - pchisq( varsGene * (m-1) / testDenomA, m-1 )
	padjA <- p.adjust( pA, "BH" )
	table( padjA < .1 )
	
	
	eps=1
	LogNcountsGene=log10(nCountsGene+eps)
	dLogNcountsGene=1/(meansGene+eps)
	var_techGene=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansGene)*meansGene^2
	LogVar_techGene=(dLogNcountsGene*sqrt(var_techGene))^2
	
	gene_names=rownames(nCountsGene)
	gene_names_het=gene_names[which(padjA<0.1)]
	
	#all Cycle base genes homologs (top 600 genes)
	hu2SYMAll <- intraIDMapper(dataCB[1:600,3],'HOMSA',srcIDType='ENSEMBL',destIDType='SYMBOL')
	ccCBall_gene_indices=na.omit(match(unlist(hu2SYMAll),rownames(nCountsGene)))
	
	xxGO <- as.list(org.Hs.egGO2EG)
	cell_cycleEG <-unlist(xxGO['GO:0007049'])
	x <- org.Hs.egENSEMBL
	mapped_genes <- mappedkeys(x)
	xxE <- as.list(x[mapped_genes])
	ens_ids_cc<-unlist(xxE[cell_cycleEG])
	ens_to_sym_cc <- intraIDMapper(ens_ids_cc,'HOMSA',srcIDType='ENSEMBL',destIDType='SYMBOL')
	cc_gene_indices <- na.omit(match(unlist(ens_to_sym_cc), rownames(countsGene)))
	
	cellcyclegenes <- ens_ids_cc
	cellcyclegenes_filter <- cc_gene_indices
	cell_names <- colnames(nCountsGene)
	Y <- nCountsGene
	genes_heterogen <- (padjA<0.1)*1
	countsERCC_mat=as.matrix(countsERCC * 1)
	countsGene_mat = as.matrix(countsGene * 1)
	
	#plot mean/cv2 relationship and 
	plot( meansGene, cv2Gene, log="xy", col=1+(padjA<0.1),ylim=c(0.1,95), xlab='Mean Counts', ylab='CV2 Counts')
	xg <- 10^seq( -3, 5, length.out=100 )
	lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg,lwd=2,col='blue' )
	points(meansERCC, cv2ERCC,col='blue',pch=15,cex=1.1)
	points(meansGene[cc_gene_indices], cv2Gene[cc_gene_indices],col=rgb(0,255,0,100,maxColorValue=255),pch=2,cex=0.75)
	points(meansGene[ccCBall_gene_indices],cv2Gene[ccCBall_gene_indices],col=rgb(0,255,0,20,maxColorValue=255),pch=2,cex=0.8)
	legend('bottomleft',c('cells (padj >= 0.1)','cells (padj<0.1)','ERCC','Cell Cycle genes'),pch=c(1,1,15),col=c('black','red','blue','green'),cex=0.7)
	
	h5save(ccCBall_gene_indices, gene_names, cellcyclegenes_filter, cellcyclegenes, cell_names, nCountsGene, genes_heterogen, LogVar_techGene, LogNcountsGene, countsGene_mat, sfERCC,countsERCC_mat, file=paste('data', base_name, '_normCounts.h5f',sep='_'))
}