# source("../../compbio_tools/H_DESeq_analysis.R")
# source("H_DESeq_analysis.R")
library("DESeq")
library("RColorBrewer")
library("gplots")
library('corpcor')
#library('nFactors')
#library("ggplot2")

# Get arguments from script
#raw_args = commandArgs()
#args = raw_args[6:length(raw_args)]

dataPrep <- function() {
	# Create source data
	print("Creating source data")
	ensemblTable <<- read.table("r_table_allalign_quick.txt",header=TRUE,row.names=1,sep="\t")
	condition <<-  factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','control','differentiated','hES','teratoma','control','organoid','teratoma','control','hES','organoid','organoid','control','differentiated','organoid','organoid','organoid','organoid','background','teratoma','control','control','differentiated'))
#	condition <<-  factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','control','contaminated','hES','teratoma','control','organoid','teratoma','control','hES','organoid','contaminated','control','differentiated','organoid','organoid','organoid','organoid','background','teratoma','control','control','differentiated'))
#	condition <<-  factor( c( "organoid-DH8","differentiated-DH5","organoid-DH12","hES-DH6","organoid-DH22",'hES-DH7','background-DH25','teratoma-DH9','control-HC1','differentiated-DH2','hES-DH3','teratoma-DH13','control-HC5','organoid-DH20','teratoma-DH14','control-HC4','hES-DH26','organoid-DH16','organoid-DH1','control-HC6','differentiated-DH23','organoid-DH15','organoid-DH18','organoid-DH10','organoid-DH4','background-DH24','teratoma-DH11','control-HC2','control-HC3','differentiated-DH21'))

	# Prepare data
	print("Preparing data")
	ensembl_cds <<- newCountDataSet(ensemblTable,condition)
	ensembl_cds_size <<- estimateSizeFactors(ensembl_cds)
	ensembl_cds_disp <<- estimateDispersions(ensembl_cds_size)

	# Filter data
#	print("Filtering data")
#	rs = rowSums ( counts ( ensembl_cds_disp ))
#	data_filter = (rs > quantile(rs, probs=0.3))
#	ensembl_filt <<- ensembl_cds_disp[ data_filter, ]
	#hist(ensembl_filt_binom$pval, breaks=100, col="skyblue", border="slateblue", main="")

	# Blind data
	print("Blind data prep")
	ensembl_cds_blind <<- estimateDispersions( ensembl_cds_size, method="blind" )
	ensembl_vsd <<- varianceStabilizingTransformation( ensembl_cds_blind )
}

examineFilter <- function() {
	# Compare p-values of filtered data
	print("p-values of filtered data")
	pdf('r_histogram_0.3.pdf')
	h1 = hist(ensembl_cds_binom$pval[!data_filter], breaks=50, plot=FALSE)
	h2 = hist(ensembl_cds_binom$pval[data_filter], breaks=50, plot=FALSE)
	colori = c('do not pass'="khaki",'pass'="powderblue")
	barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")
	text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
	legend("topright", fill=rev(colori), legend=rev(names(colori)))
}

outputTables <- function() {
	# Export size factors
	print("Writing sizes to file")
	write.table(sizeFactors(ensembl_cds_disp), file = "r_ensembl_size.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
	print("Done writing table")

	# Export expression results of processing
	print("Writing vsd expressions to file")
	write.table(exprs(ensembl_vsd), file = "r_ensembl_vsd_no21.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
	print("Done writing table")
}

binomDist <- function() {
	# Graph binomial distribution
	print("Binomial distributions")
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "control","organoid" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="control_notDH1DH2_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "background","organoid" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="background_notDH1DH2_all_genes.csv" )

#	hist(ensembl_cds_binom$pval, breaks=100, col="skyblue", border="slateblue", main="")
#	plotMA(ensembl_cds_binom)

	# List significant genes
#	print("Significantly different genes")
#	ensembl_sig = ensembl_cds_binom[ensembl_cds_binom$padj< 0.1,]
#	write.csv( ensembl_sig[ order(ensembl_sig$pval), ], file="organoid_control_sig_genes_0.1.csv" )
}

allCombos <- function() {
	# Get all combinations of data comparisons
	print("All data combinations")
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","teratoma" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_teratoma_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","organoid" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_organoid_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","differentiated" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_differentiated_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","organoid" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_organoid_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","differentiated" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_differentiated_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","differentiated" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_differentiated_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "differentiated","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="differentiated_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "differentiated","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="differentiated_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "control","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="control_background_all_genes.csv" )
}

heatmapPlot <- function() {
	# Display heatmap
	print("Heatmap")
	select = order(rowMeans(counts(ensembl_cds_size)), decreasing=TRUE)[1:5000]
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
	pdf('r_heatmap_quick.pdf')
	heatmap.2(exprs(ensembl_vsd)[select,], dendrogram = c("column"),col = hmcol, trace="none", labRow = "", margin=c(10, 6))
	dev.off()
#	heatmap.2((exprs(ensembl_vsd))[select,], dendrogram = c("column"),col = hmcol, rowVal = "", trace="none", margin=c(10, 6))
}

heatmapGenes <- function() {
	# Display heatmap
	print("Heatmap")
	# Prep gene list by replacing newlines with \\b","\\b
	geneList = paste(c("\\bAXIN2\\b","\\bDUOX2\\b","\\bKLF5\\b","\\bKRT20\\b","\\bLGR5\\b","\\bMEP1A\\b","\\bMGAM\\b","\\bMUC13\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bSLC15A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bVIL1\\b","\\bVILL\\b","\\bREG4\\b","\\bMUC17\\b","\\bSLC15A1\\b","\\bEPCAM\\b","\\bGUCA2A\\b","\\bCLCA4\\b","\\bTFF2\\b"),collapse="|")
	rows = grep(geneList,featureNames(ensembl_vsd),value=FALSE)
	cols_cat <- c(19,25,1,22,24,18,3,23,14,5,10,2,30,21) #categorized grouping of org and diff samples
	cols_pair <- c(19,10,25,2,14,30,5,21) #pairwise grouping of org and diff samples
	selected_vsd = exprs(ensembl_vsd)[rows,cols_pair]
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#	norm_vsd = 100*selected_vsd/rowSums(selected_vsd)

	pdf('r_heatmap_onlypairs.pdf',width=10, height=10)
	heatmap.2(selected_vsd, dendrogram = c("none"),col = hmcol, trace="none", margin=c(15, 6), Colv = F)
	dev.off()
#	dev.new(height=10)
#	heatmap.2(selected_vsd, dendrogram = c("column"),col = hmcol, trace="none", margin=c(10, 6))
}

similarityPlot <- function() {
	# Display similarity heatmap
	print("Similarity matrix")
	dists = dist( t( exprs(ensembl_vsd) ) )
	mat = as.matrix( dists )
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
	rownames(mat) = colnames(mat) = with(pData(ensembl_vsd), paste(colnames(ensemblTable), sep=" : "))
	pdf('r_similarity_quick.pdf')
	heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
	dev.off()
#	heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
}

deseqPCA <- function() {
	# Principle Component Analysis plot
	print("DESeq PCA")
	pdf('r_PCA_quick.pdf')
	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=120000))
	dev.off()
#	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=120000)) # 42859 genes, 120000 transcripts
}

rFactorAnalysis <- function() {
	# R Factor Analysis
	fit <- factanal(exprs(ensembl_vsd), 5, rotation="varimax")
	print(fit, digits=2, cutoff=.3, sort=TRUE)
	# plot factor 1 by factor 2
	load <- fit$loadings[,1:2]
	cond_color =  factor( c( "orange","yellow","orange","red","orange",'red','black','pink','blue','yellow','red','pink','blue','orange','pink','blue','red','orange','orange','blue','yellow','orange','orange','orange','orange','black','pink','blue','blue'))
	pdf('r_FA_3factors_genes.pdf')
	plot(load,) # set up plot
	text(load,labels = colnames(exprs(ensembl_vsd)),col=as.character(cond_color),cex=.7) # add variable names
	dev.off()
	plot(load,type="n") # set up plot
	text(load,labels = colnames(exprs(ensembl_vsd)),col=as.character(cond_color),cex=.7) # add variable names
}

prcomp.shrink =	function( X ){

  n = dim(X)[1];
  cent = sweep(X,2,colMeans(X),"-");
  covar = cov.shrink( cent );
  M = eigen( covar );

  to_ret = list( sdev = sqrt( M$values[1:n] ), rotation = M$vectors[,1:n] );
  to_ret$x = cent %*% to_ret$rotation;

  return( to_ret );

}

rPCA <- function() {
	# R PCA
	print("R PCA w/out VSD")
	ncounts <- t( t(counts(ensembl_cds_blind)) / sizeFactors(ensembl_cds_blind) ) # does not include dispersions!
	blah <<- cov.shrink(ncounts)
	#print blah

	#pca <- prcomp.shrink(t(ncounts))
	#print(summary(pca)) # print variance accounted for
	#print(loadings(pca)) # pc loadings
#	plot(pca,type="lines") # scree plot
	#print(pca$scores) # the principal components
#	biplot(pca)

#	pdf('r_PCA_withnames.pdf')
	#plot(pca$x[,1],pca$x[,2],type="p",col=as.character('red'))
	#text(pca$x[,1],pca$x[,2],labels=rownames(pca$x))
#	dev.off()

	# Export component and loading information
#	print("Save PCA components")
#	write.table(pca$loadings, file = "r_PCA_loadings_genes.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
#	write.table(pca$scores, file = "r_PCA_scores_genes.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
}

calcDist <- function() {
	# Calculate distances based on PCA
	print("Calculating MDS")
	d = dist(sweep(fit$x,2,fit$sdev,"/")[,1:2])
	mds <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
	mds <- isoMDS(d, k=2) # k is the number of dim

	# Calculate distances based on factor analysis
	print("Calculating MDS")
	d = dist(sweep(fit$x,2,fit$sdev,"/")[,1:2])
	mds <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
	mds <- isoMDS(d, k=2) # k is the number of dim

	# plot distances
	x <- mds$points[,1]
	y <- mds$points[,2]
	ggplot(scopx, aes(x=x, y=y))+ geom_point(aes(colour=condition))+labs(colour = "Samples")
	cond_color =  factor( c( "orange","yellow","orange","red","orange",'red','black','pink','blue','yellow','red','pink','blue','orange','pink','blue','red','orange','orange','blue','yellow','orange','orange','orange','orange','black','pink','blue','blue'))
	pdf('r_isomds_5dim_plot.pdf')
	plot(x, y)
	text(x, y,  col=as.character(cond_color), labels = colnames(exprs(ensembl_vsd)), cex=.7)
	dev.off()
	plot(x, y)
	text(x, y,  col=as.character(cond_color), cex=.7)

	# Export distances
	print("Save distances")
	write.table(fit$loadings, file = "r_PCA_loadings_no21.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
}

calcProb <- function() {
	# probability of overlap between lists
	total = 42859
#	total = 2000
	sam1 = 585
	sam2 = 1364
	overlap = 381

	prob_smaller = phyper(overlap,sam1,total-sam1,sam2)
	prob_larger = 1- phyper(overlap-1,sam1,total-sam1,sam2)
	print(prob_smaller)
	print(prob_larger)
	print(format(prob_smaller, scientific = TRUE, digits = 10))
}

#dataPrep()

#binomDist()
heatmapPlot()
#heatmapGenes()
similarityPlot()
deseqPCA()
#allCombos()
#rPCA()
#calcProb()
