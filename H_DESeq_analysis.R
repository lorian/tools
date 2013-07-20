library("DESeq")
library("RColorBrewer")
library("gplots")

# Get arguments from script
#raw_args = commandArgs()
#args = raw_args[6:length(raw_args)]

# Create source data
#ensemblTable = read.table("r_table_ensembl.txt",header=TRUE,row.names=1,sep="\t")
#condition =  factor( c( "organoid","organoid","organoid-DH12","organoid","hES","hES","control","teratoma","control","teratoma","organoid-DH1","differentiated","teratoma","control","hES","control","control","teratoma","control","organoid","differentiated") )

# Prepare data
#ensembl_cds = newCountDataSet(ensemblTable,condition)
#ensembl_cds_size = estimateSizeFactors(ensembl_cds)
#ensembl_cds_disp = estimateDispersions(ensembl_cds_size)

# Filter data
#rs = rowSums ( counts ( ensembl_cds_disp ))
#data_filter = (rs > quantile(rs, probs=0.3))
#ensembl_filt = ensembl_cds_disp[ data_filter, ]
#ensembl_filt_binom = nbinomTest( ensembl_filt, "control","organoid" )
#hist(ensembl_filt_binom$pval, breaks=100, col="skyblue", border="slateblue", main="")

# Graph binomial distribution
#ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "control","organoid" )
#hist(ensembl_cds_binom$pval, breaks=100, col="skyblue", border="slateblue", main="")
#plotMA(ensembl_cds_binom)

# Compare p-values of filtered data
#pdf('r_histogram_0.3.pdf')
#h1 = hist(ensembl_cds_binom$pval[!data_filter], breaks=50, plot=FALSE)
#h2 = hist(ensembl_cds_binom$pval[data_filter], breaks=50, plot=FALSE)
#colori = c('do not pass'="khaki",'pass'="powderblue")
#barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")
#text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
#legend("topright", fill=rev(colori), legend=rev(names(colori)))
#dev.off()

# Display heatmap
#ensembl_cds_blind = estimateDispersions( ensembl_cds_size, method="blind" )
#ensembl_vsd = varianceStabilizingTransformation( ensembl_cds_blind )
select = order(rowMeans(counts(ensembl_cds_size)), decreasing=TRUE)[1:1000]
#hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#heatmap.2(exprs(ensembl_vsd)[select,], dendrogram = c("column"),col = hmcol, trace="none", labRow = "", margin=c(10, 6))

# Display similarity heatmap
#dists = dist( t( exprs(ensembl_vsd) ) )
#mat = as.matrix( dists )
#rownames(mat) = colnames(mat) = with(pData(ensembl_filt), paste(colnames(ensemblTable), sep=" : "))
#heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# Principle Component Analysis plot
#print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=100000))

# R PCA
fit <- princomp(exprs(ensembl_vsd))
print(summary(fit)) # print variance accounted for
#print(loadings(fit)) # pc loadings
#plot(fit,type="lines") # scree plot
#print(fit$scores) # the principal components
biplot(fit)
# Export component and loading information
#write.table(fit$loadings, file = "r_PCA_loadings.txt", row.names=TRUE, col.names=TRUE, sep = "\t")
#write.table(fit$scores, file = "r_PCA_scores.txt", row.names=TRUE, col.names=TRUE, sep = "\t")

