# source("../../compbio_tools/H_DESeq_analysis.R")
# source("H_DESeq_analysis.R")
library("DESeq")
library("RColorBrewer")
library("gplots")
library('corpcor')
library('nFactors')
library('psych')
library('psy')
#library("ggplot2")

# Get arguments from script
#raw_args = commandArgs()
#args = raw_args[6:length(raw_args)]

dataPrep <- function() {
	# Create source data
	print("Creating source data")
	ensemblTable <<- read.table("r_table_noa_newdata_genenames.txt",header=TRUE,row.names=1,sep="\t")

	condition <<- factor (c("background","organoid","organoid","organoid","teratoma","organoid","control","teratoma","differentiated","control","hES","hES","teratoma","organoid","organoid","organoid","new_mixed","organoid","organoid","teratoma","differentiated","background","hES","new_differentiated","organoid","differentiated","control","control","differentiated","hES","new_differentiated","control","new_mixed","new_organoid","new_organoid","control")) #_noa_newdata_genenames new samples
#	condition <<- factor (c("background","organoid","organoid","organoid","teratoma","organoid","control","teratoma","differentiated","control","hES","hES","teratoma","organoid","organoid","organoid","Nterm","organoid","organoid","teratoma","differentiated","background","hES","Nterm","organoid","differentiated","control","control","differentiated","hES","Cterm","control","Cterm","Nterm","Cterm","control")) #_noa_newdata_genenames new samples
#	condition <<- factor (c("background","organoid","organoid","organoid","teratoma","organoid","control","teratoma","differentiated","control","hES","hES","teratoma","organoid","organoid","organoid","mixed","organoid","organoid","teratoma","differentiated","background","hES","differentiated","organoid","differentiated","control","control","differentiated","hES","differentiated","control","mixed","organoid","organoid","control")) #_noa_newdata_genenames
#	condition <<- factor (c("background","organoid","organoid","organoid","teratoma","organoid","rectum","teratoma","differentiated","duodenum","hES","hES","teratoma","organoid","organoid","organoid","mixed","organoid","organoid","teratoma","differentiated","background","hES","differentiated","organoid","differentiated","ileum","ileum","differentiated","hES","differentiated","rectum","mixed","organoid","organoid","duodenum")) #_noa_newdata_genenames w/specific controls
#	condition <<- factor( c( "background","organoid","organoid","organoid","teratoma","organoid","control","teratoma","differentiated","control","hES","hES","teratoma","organoid","organoid","liver","organoid","organoid","organoid","teratoma","differentiated","background","hES","organoid","differentiated","control","breast","control","differentiated","hES","control","colon","control")) #full _noa w/IBM
#	condition <<- factor(c("organoid","organoid","differentiated","control","background","organoid","teratoma","hES","organoid","differentiated","hES","organoid","control","differentiated","organoid","teratoma","control","organoid","organoid","teratoma","control","organoid","control","teratoma","organoid","differentiated","background","control","hES","hES")) #_quick
#	condition <<- factor( c( "organoid","teratoma","background","organoid","hES","hES","teratoma","organoid","organoid","differentiated","control","differentiated","differentiated","organoid","organoid","differentiated","teratoma","control","hES","teratoma","control","organoid","control","organoid","control","background","organoid","organoid","hES")) #_a without HC2
#	condition <<- factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','control','differentiated','hES','teratoma','control','organoid','teratoma','control','hES','organoid','organoid','control','differentiated','organoid','organoid','organoid','organoid','background','teratoma','control','control','differentiated')) #_newfused
#	condition <<- factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','duodenum','differentiated','hES','teratoma','rectum','organoid','teratoma','ileum','hES','organoid','organoid','rectum','differentiated','organoid','organoid','organoid','organoid','background','teratoma','duodenum','ileum','differentiated')) #_newfused w/specific controls
#	condition <<- factor( c( "organoid-DH8","differentiated-DH5","organoid-DH12","hES-DH6","organoid-DH22",'hES-DH7','background-DH25','teratoma-DH9','control-HC1','differentiated-DH2','hES-DH3','teratoma-DH13','control-HC5','organoid-DH20','teratoma-DH14','control-HC4','hES-DH26','organoid-DH16','organoid-DH1','control-HC6','differentiated-DH23','organoid-DH15','organoid-DH18','organoid-DH10','organoid-DH4','background-DH24','teratoma-DH11','control-HC2','control-HC3','differentiated-DH21'))

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
	write.table(sizeFactors(ensembl_cds_disp), file = "r_noa_ensembl_size.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
	print("Done writing table")

	# Export expression results of processing
	print("Writing vsd expressions to file")
	write.table(exprs(ensembl_vsd), file = "r_noa_ensembl_vsd.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
	print("Done writing table")
}

binomDist <- function() {
	# Graph binomial distribution
	print("Binomial distributions")
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "new_organoid","new_mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_organoid_mixed_new_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "new_organoid","new_differentiated" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_organoid_differentiated_new_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "new_differentiated","new_mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_differentiated_mixed_new_genes.csv" )

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
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_hES_teratoma_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","organoid" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_hES_organoid_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","differentiated" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_hES_differentiated_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_hES_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_hES_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_hES_mixed_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","organoid" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_teratoma_organoid_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","differentiated" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_teratoma_differentiated_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_teratoma_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_teratoma_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_teratoma_mixed_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","differentiated" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_organoid_differentiated_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_organoid_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_organoid_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_organoid_mixed_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "differentiated","control" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_differentiated_control_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "differentiated","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_differentiated_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "differentiated","mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_differentiated_mixed_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "control","background" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_control_background_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "control","mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_control_mixed_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "background","mixed" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="decemberdata_background_mixed_all_genes.csv" )

}

heatmapPlot <- function() {
	# Display heatmap
	print("Heatmap")
	select = order(rowMeans(counts(ensembl_cds_size)), decreasing=TRUE)[1:5000]
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
	pdf('decemberdata_r_heatmap_genes.pdf')
	heatmap.2(exprs(ensembl_vsd)[select,], dendrogram = c("column"),col = hmcol, trace="none", labRow = "", margin=c(10, 6))
	dev.off()
	heatmap.2(exprs(ensembl_vsd)[select,], dendrogram = c("column"),col = hmcol, trace="none", labRow = "", margin=c(10, 6))
}

heatmapGenes <- function(listname, rawgenelist) {
	# Display heatmap. TABLE MUST USE GENE NAMES, NOT TRANSCRIPT NAMES
	print("Heatmap")
	# Prep gene list by replacing newlines with \\b","\\b
	geneList = paste(rawgenelist,collapse="|")
	rows = grep(geneList,featureNames(ensembl_vsd),value=FALSE)
	cols_nodiff = c(1,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29)
	cols_orgdiff = c(1,2,3,5,10,14,18,19,21,22,23,24,25,30)
	cols_paired = c(19,10,25,2,14,30,5,21)
#	selected_vsd = exprs(ensembl_vsd)[rows,cols]
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#	norm_vsd = 100*selected_vsd/rowSums(selected_vsd)

	pdf(paste('r_heatmap_noa_',listname,'.pdf'),width=10, height=10)
	heatmap.2(exprs(ensembl_vsd)[rows,],col = hmcol, trace="none", margin=c(15, 6),dendrogram = c("both"),)
	dev.off()
	dev.new(height=10)
	heatmap.2(exprs(ensembl_vsd)[rows,],col = hmcol, trace="none", margin=c(15, 6),dendrogram = c("both"),)
#	heatmap.2(exprs(ensembl_vsd)[rows,cols_paired],col = hmcol, trace="none", margin=c(15, 6),dendrogram = c("row"), Colv=FALSE) # no column reordering

}

similarityPlot <- function() {
	# Display similarity heatmap
	print("Similarity matrix")
	dists = dist( t( exprs(ensembl_vsd) ) )
	mat = as.matrix( dists )
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
	rownames(mat) = colnames(mat) = with(pData(ensembl_vsd), paste(colnames(ensemblTable), sep=" : "))
	pdf('r_similarity_a.pdf')
	heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
	dev.off()
#	heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
}

screePlot <- function() {
	# Determine Number of Factors to Extract


	#ev <- eigen(cor(exprs(ensembl_vsd))) # get eigenvalues
	#ap <- parallel(subject=nrow(exprs(ensembl_vsd)),var=ncol(exprs(ensembl_vsd)),
	#  rep=100,cent=.05)
	#nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
	#plotnScree(nS)
}

deseqPCA <- function() {
	# Principle Component Analysis plot
	print("DESeq PCA")

	pdf('decemberdata_r_PCA.pdf')
	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=42859))
	dev.off()
	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=42859)) # 42859 genes, 120000 transcripts
}

rFactorAnalysis <- function() {
	# R Factor Analysis
	fit <- factanal(exprs(ensembl_vsd), 3, rotation="varimax")
    scree.plot(fit$correlations)
	print(fit, digits=2, cutoff=.3, sort=TRUE)
	# plot factor 1 by factor 2
	load <- fit$loadings[,1:2]

#	cond_color = factor(c("orange","orange","yellow","black","pink","orange","red","blue","orange","yellow","blue","orange","black","yellow","orange","red","black","orange","orange","red","black","orange","black","red","orange","yellow","pink","black","blue","blue")) #_quick
	cond_color = factor( c( "orange","red","pink","orange","blue","blue","red","orange","orange","yellow","black","yellow","yellow","orange","orange","yellow","red","black","blue","red","black","orange","black","orange","black","pink","orange","orange","blue")) #_a without HC2

#	pdf('r_FA_3factors_quick.pdf')
#	plot(load,) # set up plot
#	text(load,labels = colnames(exprs(ensembl_vsd)),col=as.character(cond_color),cex=.7) # add variable names
#	dev.off()
#	plot(load,type="n") # set up plot
#	text(load,labels = colnames(exprs(ensembl_vsd)),col=as.character(cond_color),cex=.7) # add variable names
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

getVariances <- function() {
	# Get variances per condition
	organoid_cols = which( condition == 'organoid' )
	differentiated_cols = which( condition == 'differentiated' )
	control_cols = which( condition == 'control' )
	hES_cols = which( condition == 'hES' )
	background_cols = which( condition == 'background' )
	teratoma_cols = which( condition == 'teratoma' )

	# Check data is correct
	print(head(counts( ensembl_cds_size, normalized=TRUE )["OPA3",organoid_cols]))
	print(head(genefilter::rowVars( counts( ensembl_cds_size, normalized=TRUE )["OPA3",organoid_cols])))
	print(head(rowMeans(counts(ensembl_cds_size, normalized=TRUE)["OPA3",organoid_cols])))

	print("Wrinting variances to file")
	organoid_var = genefilter::rowVars( counts( ensembl_cds_size, normalized=TRUE )[,organoid_cols])
	differentiated_var = genefilter::rowVars( counts( ensembl_cds_size, normalized=TRUE )[,differentiated_cols])
	control_var = genefilter::rowVars( counts( ensembl_cds_size, normalized=TRUE )[,control_cols])
	hES_var = genefilter::rowVars( counts( ensembl_cds_size, normalized=TRUE )[,hES_cols])
	background_var = genefilter::rowVars( counts( ensembl_cds_size, normalized=TRUE )[,background_cols])
	teratoma_var = genefilter::rowVars( counts( ensembl_cds_size, normalized=TRUE )[,teratoma_cols])

	var_table <- t(rbind(organoid_var,differentiated_var,control_var,hES_var,background_var,teratoma_var))
	write.table(var_table, file = "variances.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
	print("Done writing variances")

	print("Writing means to file")
	organoid_m = rowMeans(counts(ensembl_cds_size, normalized=TRUE)[,organoid_cols])
	differentiated_m = rowMeans(counts(ensembl_cds_size, normalized=TRUE)[,differentiated_cols])
	control_m = rowMeans(counts(ensembl_cds_size, normalized=TRUE)[,control_cols])
	hES_m = rowMeans(counts(ensembl_cds_size, normalized=TRUE)[,hES_cols])
	background_m = rowMeans(counts(ensembl_cds_size, normalized=TRUE)[,background_cols])
	teratoma_m = rowMeans(counts(ensembl_cds_size, normalized=TRUE)[,teratoma_cols])

	var_table <- t(rbind(organoid_m,differentiated_m,control_m,hES_m,background_m,teratoma_m))
	write.table(var_table, file = "means.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
	print("Done writing means")
}

dataPrep()

binomDist()
#heatmapPlot()
#similarityPlot()
#deseqPCA()
#allCombos()
#screePlot()
#rPCA()
#rFactorAnalysis()
#calcProb()
#outputTables()
#getVariances()

#heatmapGenes("aquaporins",c("\\bAQP1\\b","\\bAQP2\\b","\\bAQP3\\b","\\bAQP4\\b","\\bAQP5\\b","\\bAQP6\\b","\\bAQP7\\b","\\bAQP8\\b","\\bAQP9\\b","\\bAQP10\\b","\\bAQP11\\b"))
#heatmapGenes("maturity",c("\\bAXIN2\\b","\\bCDH1\\b","\\bCDH17\\b","\\bCDX2\\b","\\bDEFA5\\b","\\bDEF6\\b","\\bDUOX2\\b","\\bELF3\\b","\\bEPCAM\\b","\\bKLF4\\b","\\bKLF5\\b","\\bLGR5\\b","\\bLYZ\\b","\\bMUC13\\b","\\bMUC2\\b","\\bSOX9\\b","\\bVIL1\\b","\\bSI\\b","\\bMMP7\\b","\\bPLA2G2A\\b","\\bWNT3A\\b","\\bATOH1\\b","\\bTFF3\\b","\\bTERT\\b"))
#heatmapGenes("adultstem",c("\\bLGR5\\b","\\bSOX9\\b","\\bKLF4\\b","\\bKLF5\\b","\\bOLFM4\\b","\\bTERT\\b"))
#heatmapGenes("pattern",c("\\bBARX1\\b","\\bSFRP1\\b","\\bHHEX\\b","\\bSOX2\\b","\\bDKK1\\b","\\bFABP2\\b","\\bCDX1\\b","\\bCDX2\\b","\\bHOXC5\\b"))
#heatmapGenes('differentiated',c("\\bANPEP\\b","\\bASCL2\\b","\\bCA1\\b","\\bLGR5\\b","\\bMKI67\\b","\\bMUC2\\b","\\bTFF3\\b","\\bAXIN2\\b","\\bCLCA4\\b","\\bDUOX2\\b","\\bEPCAM\\b","\\bGUCA2A\\b","\\bKLF5\\b","\\bLGR5\\b","\\bMEP1A\\b","\\bMGAM\\b","\\bMGAM\\b","\\bMUC13\\b","\\bMUC17\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bREG4\\b","\\bSLC15A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF2\\b","\\bVIL1\\b","\\bVILL\\b"))
#heatmapGenes('secretory_new',c("\\bAXIN2\\b","\\bCDH1\\b","\\bCDH17\\b","\\bCDX2\\b","\\bDEF6\\b","\\bDUOX2\\b","\\bELF3\\b","\\bEPCAM\\b","\\bKLF4\\b","\\bKLF5\\b","\\bLGR5\\b","\\bLYZ\\b","\\bMUC13\\b","\\bMUC2\\b","\\bSOX9\\b","\\bVIL1\\b","\\bSI\\b","\\bMMP7\\b","\\bPLA2G2A\\b","\\bWNT3A\\b","\\bTFF3 \\b","\\bTERT\\b","\\bCA2\\b","\\bSCNN1A\\b","\\bSLC5A1\\b"))
#heatmapGenes('adultstem_new',c("\\bLGR5\\b","\\bSOX9\\b","\\bKLF4\\b","\\bKLF5\\b","\\bOLFM4\\b","\\bTERT\\b"))
#heatmapGenes('orgonly1',c("\\bAXIN2\\b","\\bCLCA4\\b","\\bDUOX2\\b","\\bEPCAM\\b","\\bGUCA2A\\b","\\bKLF5\\b","\\bKRT20\\b","\\bKRT20\\b","\\bLGR5\\b","\\bMEP1A\\b","\\bMGAM\\b","\\bMGAM\\b","\\bMUC13\\b","\\bMUC17\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bREG4\\b","\\bSLC15A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF2\\b","\\bVIL1\\b","\\bCA2\\b","\\bSCNN1A\\b"))
#heatmapGenes('orgonly2',c("\\bANPEP\\b","\\bASCL2\\b","\\bCA1\\b","\\bLGR5\\b","\\bMKI67\\b","\\bMUC2\\b","\\bTFF3\\b","\\bAXIN2\\b","\\bCLCA4\\b","\\bDUOX2\\b","\\bEPCAM\\b","\\bGUCA2A\\b","\\bKLF5\\b","\\bLGR5\\b","\\bMEP1A\\b","\\bMGAM\\b","\\bMGAM\\b","\\bMUC13\\b","\\bMUC17\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bREG4\\b","\\bSLC15A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF2\\b","\\bVIL1\\b","\\bVILL\\b"))
#heatmapGenes('orgonlyfinal',c("\\bAXIN2\\b","\\bCA2\\b","\\bCLCA4\\b","\\bDUOX2\\b","\\bEPCAM\\b","\\bGUCA2A\\b","\\bKLF5\\b","\\bKRT20\\b","\\bKRT20\\b","\\bLGR5\\b","\\bMEP1A\\b","\\bMGAM\\b","\\bMGAM\\b","\\bMUC13\\b","\\bMUC17\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bREG4\\b","\\bSCNN1A\\b","\\bSLC15A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF2\\b","\\bVIL1\\b","\\bANPEP\\b","\\bASCL2\\b","\\bMKI67\\b","\\bTFF3\\b","\\bVILL\\b"))
#heatmapGenes('column1',c("\\bAXIN2\\b","\\bCA2\\b","\\bCDH1\\b","\\bCDH17\\b","\\bCDX2\\b","\\bDEF6\\b","\\bDUOX2\\b","\\bELF3\\b","\\bEPCAM\\b","\\bKLF4\\b","\\bKLF5\\b","\\bKRT20\\b","\\bLGR5\\b","\\bLYZ\\b","\\bMMP7\\b","\\bMUC13\\b","\\bMUC2\\b","\\bPLA2G2A\\b","\\bSCNN1A\\b","\\bSI\\b","\\bSLC5A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF3 \\b","\\bVIL1\\b","\\bWNT3A\\b"))
#heatmapGenes('column2',c("\\bAXIN2\\b","\\bCA2\\b","\\bCDH1\\b","\\bCDH17\\b","\\bCDX2\\b","\\bCFTR\\b","\\bDEF6\\b","\\bDUOX2\\b","\\bELF3\\b","\\bEPCAM\\b","\\bKLF4\\b","\\bKLF5\\b","\\bKRT20\\b","\\bLGR5\\b","\\bLYZ\\b","\\bMMP7\\b","\\bMUC13\\b","\\bMUC2\\b","\\bPLA2G2A\\b","\\bSCNN1A\\b","\\bSI\\b","\\bSLC5A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF3 \\b","\\bVIL1\\b","\\bWNT3A\\b"))
#heatmapGenes('column3',c("\\bAXIN2\\b","\\bCA2\\b","\\bCDH1\\b","\\bCDH17\\b","\\bCDX2\\b","\\bCFTR\\b","\\bDEF6\\b","\\bDUOX2\\b","\\bELF3\\b","\\bEPCAM\\b","\\bKLF4\\b","\\bKLF5\\b","\\bKRT20\\b","\\bLGR5\\b","\\bLYZ\\b","\\bMMP7\\b","\\bMUC13\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bPLA2G2A\\b","\\bSCNN1A\\b","\\bSI\\b","\\bSLC5A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF3 \\b","\\bVIL1\\b","\\bWNT3A\\b"))
#heatmapGenes('column4',c("\\bAXIN2\\b","\\bCA2\\b","\\bCDH1\\b","\\bCDH17\\b","\\bCDX2\\b","\\bDEF6\\b","\\bDUOX2\\b","\\bELF3\\b","\\bEPCAM\\b","\\bKLF4\\b","\\bKLF5\\b","\\bKRT20\\b","\\bLGR5\\b","\\bLYZ\\b","\\bMMP7\\b","\\bMUC13\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bPDX1\\b","\\bPLA2G2A\\b","\\bSCNN1A\\b","\\bSI\\b","\\bSOX9\\b","\\bTERT\\b","\\bTFF3 \\b","\\bVIL1\\b","\\bWNT3A\\b"))
