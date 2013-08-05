library("DESeq")
library("RColorBrewer")
library("gplots")
#library('nFactors')
#library("ggplot2")

# Get arguments from script
#raw_args = commandArgs()
#args = raw_args[6:length(raw_args)]

# Create source data
#print("Creating source data")
ensemblTable = read.table("r_table_newfused.txt",header=TRUE,row.names=1,sep="\t")
condition =  factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','control','differentiated-DH2','hES','teratoma','control','organoid','teratoma','control','hES','organoid','organoid-DH1','control','differentiated','organoid','organoid','organoid','organoid','background','teratoma','control','control','differentiated-DH21'))
#condition =  factor( c( "organoid-DH8","differentiated-DH5","organoid-DH12","hES-DH6","organoidDH22",'hESDH7','background-DH25','teratoma-DH9','control-HC1','differentiated-DH2','hES-DH3','teratoma-DH13','control-HC5','organoid-DH20','teratoma-DH14','control-HC4','hES-DH26','organoid-DH16','organoid-DH1','control-HC6','differentiated-DH23','organoid-DH15','organoid-DH18','organoid-DH10','organoid-DH4','background-DH24','teratoma-DH11','control-HC2','control-HC3','differentiated-DH21'))

# Prepare data
print("Preparing data")
ensembl_cds = newCountDataSet(ensemblTable,condition)
ensembl_cds_size = estimateSizeFactors(ensembl_cds)
ensembl_cds_disp = estimateDispersions(ensembl_cds_size)

# Filter data
print("Filtering data")
rs = rowSums ( counts ( ensembl_cds_disp ))
data_filter = (rs > quantile(rs, probs=0.3))
ensembl_filt = ensembl_cds_disp[ data_filter, ]
#ensembl_filt_binom = nbinomTest( ensembl_filt, "control","organoid" )
#hist(ensembl_filt_binom$pval, breaks=100, col="skyblue", border="slateblue", main="")

# Get all combinations of data comparisons
print("All data combinations")
ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","organoid" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_organoid_all_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "hES","differentiated" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_differentiated_all_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","organoid" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_organoid_all_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "teratoma","differentiated" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_differentiated_all_genes.csv" )
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

ensembl_cds_binom = nbinomTest( ensembl_filt, "hES","organoid" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_organoid_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "hES","differentiated" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="hES_differentiated_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "teratoma","organoid" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_organoid_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "teratoma","differentiated" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="teratoma_differentiated_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "organoid","differentiated" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_differentiated_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "organoid","control" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_control_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "organoid","background" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_background_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "differentiated","control" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="differentiated_control_filtered_genes.csv" )
ensembl_cds_binom = nbinomTest( ensembl_filt, "differentiated","background" )
write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="differentiated_background_filtered_genes.csv" )

# Export expression results of processing
#print("Writing sizes to file")
#write.table(sizeFactors(ensembl_cds_disp), file = "r_ensembl_size.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
#print("Done writing table")

# Graph binomial distribution
#print("Binomial distributions")
#ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "control","organoid" )
#hist(ensembl_cds_binom$pval, breaks=100, col="skyblue", border="slateblue", main="")
#plotMA(ensembl_cds_binom)

# List significant genes
#print("Significantly different genes")
#ensembl_sig = ensembl_cds_binom[ensembl_cds_binom$padj< 0.1,]
#write.csv( ensembl_sig[ order(ensembl_sig$pval), ], file="organoid_control_sig_genes_0.1.csv" )

# Compare p-values of filtered data
#print("p-values of filtered data")
#pdf('r_histogram_0.3.pdf')
#h1 = hist(ensembl_cds_binom$pval[!data_filter], breaks=50, plot=FALSE)
#h2 = hist(ensembl_cds_binom$pval[data_filter], breaks=50, plot=FALSE)
#colori = c('do not pass'="khaki",'pass'="powderblue")
#barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")
#text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
#legend("topright", fill=rev(colori), legend=rev(names(colori)))

# Display heatmap
#print("Heatmap")
#ensembl_cds_blind = estimateDispersions( ensembl_cds_size, method="blind" )
#ensembl_vsd = varianceStabilizingTransformation( ensembl_cds_blind )
#jpeg('r_heatmap_no21.jpg')
#select = order(rowMeans(counts(ensembl_cds_size)), decreasing=TRUE)[1:5000]
#hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#heatmap.2(exprs(ensembl_vsd)[select,], dendrogram = c("column"),col = hmcol, trace="none", labRow = "", margin=c(10, 6))
#dev.off()

# Export expression results of processing
#print("Writing vsd expressions to file")
#write.table(exprs(ensembl_vsd), file = "r_ensembl_vsd_no21.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
#print("Done writing table")

# Export expression results of processing
#print("Writing vsd to file")
#write.table(ensembl_vsd, file = "r_ensembl_vsd.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
#print("Done writing table")

# Display similarity heatmap
#print("Similarity matrix")
#jpeg('r_simmatrix_no21.jpg')
#dists = dist( t( exprs(ensembl_vsd) ) )
#mat = as.matrix( dists )
#rownames(mat) = colnames(mat) = with(pData(ensembl_vsd), paste(colnames(ensemblTable), sep=" : "))
#heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
#dev.off()

# Principle Component Analysis plot
#print("DESeq PCA")
#print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=100000))

# R Factor Analysis
#fit <- factanal(exprs(ensembl_vsd), 5, rotation="varimax")
#print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2
#load <- fit$loadings[,1:2]
#cond_color =  factor( c( "orange","yellow","orange","red","orange",'red','black','pink','blue','yellow','red','pink','blue','orange','pink','blue','red','orange','orange','blue','yellow','orange','orange','orange','orange','black','pink','blue','blue'))
#jpeg('r_FA_3factors.jpg')
#plot(load,) # set up plot
#text(load,labels = colnames(exprs(ensembl_vsd)),col=as.character(cond_color),cex=.7) # add variable names
#dev.off()
#plot(load,type="n") # set up plot
#text(load,labels = colnames(exprs(ensembl_vsd)),col=as.character(cond_color),cex=.7) # add variable names

# Find best number of factors
#ev <- eigen(cor(exprs(ensembl_vsd))) # get eigenvalues
#ap <- parallel(subject=nrow(exprs(ensembl_vsd)),var=ncol(exprs(ensembl_vsd)), rep=100,cent=.05)
#nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
#plotnScree(nS)

# R PCA
#print("R PCA")
#fit <- prcomp(t(exprs(ensembl_vsd)))
#print(summary(fit)) # print variance accounted for
#print(loadings(fit)) # pc loadings
#plot(fit,type="lines") # scree plot
#print(fit$scores) # the principal components
#biplot(fit)

# Export component and loading information
#print("Save PCA components")
#write.table(fit$loadings, file = "r_PCA_loadings_no21.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
#write.table(fit$scores, file = "r_PCA_scores_no21.txt", row.names=TRUE, col.names=TRUE, sep = "\t")

# Calculate distances based on PCA
#print("Calculating MDS")
#d = dist(sweep(fit$x,2,fit$sdev,"/")[,1:2])
#mds <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#mds <- isoMDS(d, k=2) # k is the number of dim

# Calculate distances based on ML FA
#print("Calculating MDS")
#d = dist(sweep(fit$x,2,fit$sdev,"/")[,1:2])
#mds <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#mds <- isoMDS(d, k=2) # k is the number of dim

# plot distances
#x <- mds$points[,1]
#y <- mds$points[,2]
#ggplot(scopx, aes(x=x, y=y))+ geom_point(aes(colour=condition))+labs(colour = "Samples")
#cond_color =  factor( c( "orange","yellow","orange","red","orange",'red','black','pink','blue','yellow','red','pink','blue','orange','pink','blue','red','orange','orange','blue','yellow','orange','orange','orange','orange','black','pink','blue','blue'))
#jpeg('r_isomds_5dim_plot.jpg')
#plot(x, y)
#text(x, y,  col=as.character(cond_color), labels = colnames(exprs(ensembl_vsd)), cex=.7)
#dev.off()
#plot(x, y)
#text(x, y,  col=as.character(cond_color), cex=.7)

# Export distances
#print("Save distances")
#write.table(fit$loadings, file = "r_PCA_loadings_no21.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
