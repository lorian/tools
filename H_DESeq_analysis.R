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
	ensemblTable <<- read.table("r_table_newfused_genenames.txt",header=TRUE,row.names=1,sep="\t")

#	condition <<- factor(c("organoid","organoid","differentiated","control","background","organoid","teratoma","hES","organoid","differentiated","hES","organoid","control","differentiated","organoid","teratoma","control","organoid","organoid","teratoma","control","organoid","control","teratoma","organoid","differentiated","background","control","hES","hES")) #_quick
#	condition <<-  factor( c( "organoid","teratoma","background","organoid","hES","hES","teratoma","organoid","organoid","differentiated","control","differentiated","differentiated","organoid","organoid","differentiated","teratoma","control","hES","teratoma","control","organoid","control","organoid","control","background","organoid","organoid","hES")) #_a without HC2
<<<<<<< HEAD
#	condition <<-  factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','control','differentiated','hES','teratoma','control','organoid','teratoma','control','hES','organoid','organoid','control','differentiated','organoid','organoid','organoid','organoid','background','teratoma','control','control','differentiated')) #_newfused
	condition <<-  factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','duodenum','differentiated','hES','teratoma','rectum','organoid','teratoma','ileum','hES','organoid','organoid','rectum','differentiated','organoid','organoid','organoid','organoid','background','teratoma','duodenum','ileum','differentiated')) #_newfused w/specific controls
=======
	condition <<-  factor( c( "organoid","differentiated","organoid","hES","organoid",'hES','background','teratoma','control','contaminated','hES','teratoma','control','organoid','teratoma','control','hES','organoid','contaminated','control','differentiated','organoid','organoid','organoid','organoid','background','teratoma','control','control','differentiated')) #_genes
>>>>>>> origin/master
#	condition <<-  factor( c( "organoid-DH8","differentiated-DH5","organoid-DH12","hES-DH6","organoid-DH22",'hES-DH7','background-DH25','teratoma-DH9','control-HC1','differentiated-DH2','hES-DH3','teratoma-DH13','control-HC5','organoid-DH20','teratoma-DH14','control-HC4','hES-DH26','organoid-DH16','organoid-DH1','control-HC6','differentiated-DH23','organoid-DH15','organoid-DH18','organoid-DH10','organoid-DH4','background-DH24','teratoma-DH11','control-HC2','control-HC3','differentiated-DH21'))
#	condition <<-  factor( c( "background","organoid","organoid","organoid","teratoma","organoid","control","teratoma","differentiated","control","hES","hES","teratoma","organoid","organoid","liver","organoid","organoid","organoid","teratoma","differentiated","background","hES","organoid","differentiated","control","breast","control","differentiated","hES","control","colon","control")) #full _noa w/IBM


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
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","ileum" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_ileum_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","rectum" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_rectum_all_genes.csv" )
	ensembl_cds_binom = nbinomTest( ensembl_cds_disp, "organoid","duodenum" )
	write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file="organoid_duodenum_all_genes.csv" )

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
	pdf('r_heatmap_a_genes.pdf')
	heatmap.2(exprs(ensembl_vsd)[select,], dendrogram = c("column"),col = hmcol, trace="none", labRow = "", margin=c(10, 6))
	dev.off()
	heatmap.2(exprs(ensembl_vsd)[select,], dendrogram = c("column"),col = hmcol, trace="none", labRow = "", margin=c(10, 6))
}

heatmapGenes <- function(listname, genelist) {
	# Display heatmap. TABLE MUST USE GENE NAMES, NOT TRANSCRIPT NAMES
	print("Heatmap")
	# Prep gene list by replacing newlines with \\b","\\b
<<<<<<< HEAD
	geneList = paste(c(genelist),collapse="|")
#	geneList = paste(c("\\bCLDN6\\b","\\bEOMES\\b","\\bFABP1\\b","\\bFABP2\\b","\\bFOXA1\\b","\\bFOXA2\\b","\\bGATA4\\b","\\bGSC\\b","\\bHNF1B\\b","\\bKRT19\\b","\\bSOX17\\b","\\bSOX7\\b","\\bKIT\\b","\\bCXCR4\\b"),collapse="|") #definitive endoderm2
#	geneList = paste(c("\\bAFP\\b","\\bCTNNB1\\b","\\bGATA4\\b","\\bGATA6\\b","\\bGDF1\\b","\\bGDF3\\b","\\bHNF4A\\b","\\bMIXL1\\b","\\bSALL4\\b","\\bSOX17\\b","\\bSOX7\\b"),collapse="|") #primative endoderm2
#	geneList = paste(c("\\bCD34 \\b","\\bPROM1\\b","\\bSOX9\\b","\\bLYZ\\b","\\bMMP7\\b","\\bMUC2\\b","\\bPLA2G2A\\b","\\bLGR5\\b"),collapse="|") #adult2
#	geneList = paste(c("\\bCDX1\\b","\\bCDX2\\b","\\bPDX1\\b","\\bHHEX\\b","\\bSOX2\\b","\\bBARX1\\b","\\bSFRP1\\b","\\bCD34\\b","\\bCDX4\\b","\\bHOXC5\\b","\\bFABP2\\b","\\bGATA4\\b","\\bNKX2-1\\b","\\bDKK1\\b"),collapse="|") #patterning2
#	geneList = paste(c("\\bFKBP9\\b","\\bGLUD1\\b","\\bKRT19\\b","\\bKRT8\\b","\\bIGFBP3\\b","\\bFKBP9\\b","\\bPERP\\b","\\bCDX1\\b","\\bCDX2\\b","\\bPDX1\\b","\\bGATA4\\b","\\bSOX17\\b","\\bHHEX\\b","\\bSOX2\\b","\\bNES\\b","\\bNEUROG3\\b","\\bPAX6\\b","\\bSOX1\\b","\\bSOX10\\b","\\bLEFTY1\\b","\\bESRRG\\b","\\bBMP1\\b","\\bCXCR4\\b","\\bDKK1\\b","\\bFOXC1\\b","\\bGATA3\\b","\\bKIT\\b","\\bLYZ\\b","\\bATOH1\\b","\\bMIXL1\\b","\\bMMP7\\b","\\bMUC2\\b","\\bPLA2G2A\\b","\\bSOX9\\b","\\bT\\b","\\bWNT3A\\b","\\bHNF4A\\b","\\bONECUT1\\b","\\bONECUT2\\b","\\bHOXA1\\b","\\bHOXA2\\b","\\bHOXC5\\b","\\bLRIG1\\b","\\bCER1\\b","\\bGREM2\\b","\\bCXCR4\\b","\\bEPCAM\\b","\\bFOXA2\\b","\\bSERPINA1\\b","\\bAFP\\b","\\bALB\\b","\\bASCL2\\b","\\bATOH1\\b","\\bAXIN2\\b","\\bCD44\\b"),collapse="|") #list of all
=======
	geneList = paste(c("\\bFKBP9\\b","\\bGLUD1\\b","\\bKRT19\\b","\\bKRT8\\b","\\bIGFBP3\\b","\\bFKBP9\\b","\\bPERP\\b","\\bCDX1\\b","\\bCDX2\\b","\\bPDX1\\b","\\bGATA4\\b","\\bSOX17\\b","\\bHHEX\\b","\\bSOX2\\b","\\bNES\\b","\\bNEUROG3\\b","\\bPAX6\\b","\\bSOX1\\b","\\bSOX10\\b","\\bLEFTY1\\b","\\bESRRG\\b","\\bBMP1\\b","\\bCXCR4\\b","\\bDKK1\\b","\\bFOXC1\\b","\\bGATA3\\b","\\bKIT\\b","\\bLYZ\\b","\\bATOH1\\b","\\bMIXL1\\b","\\bMMP7\\b","\\bMUC2\\b","\\bPLA2G2A\\b","\\bSOX9\\b","\\bT\\b","\\bWNT3A\\b","\\bHNF4A\\b","\\bONECUT1\\b","\\bONECUT2\\b","\\bHOXA1\\b","\\bHOXA2\\b","\\bHOXC5\\b","\\bLRIG1\\b","\\bCER1\\b","\\bGREM2\\b","\\bCXCR4\\b","\\bEPCAM\\b","\\bFOXA2\\b","\\bSERPINA1\\b","\\bAFP\\b","\\bALB\\b","\\bASCL2\\b","\\bATOH1\\b","\\bAXIN2\\b","\\bCD44\\b"),collapse="|") #list of all
>>>>>>> origin/master
#	geneList = paste(c("\\bFKBP9\\b","\\bGLUD1\\b","\\bKRT19\\b","\\bKRT8\\b","\\bIGFBP3\\b","\\bFKBP9\\b","\\bPERP\\b"),collapse="|") #shortlist dark
#	geneList = paste(c("\\bCDX1\\b","\\bCDX2\\b","\\bPDX1\\b","\\bGATA4\\b","\\bSOX17\\b","\\bHHEX\\b","\\bSOX2\\b"),collapse="|") #shortlist pattern
#	geneList = paste(c("\\bNES\\b","\\bNEUROG3\\b","\\bPAX6\\b","\\bSOX1\\b","\\bSOX10\\b"),collapse="|") #shortlist neuronal
#	geneList = paste(c("\\bLEFTY1\\b","\\bESRRG\\b","\\bBMP1\\b","\\bCXCR4\\b","\\bDKK1\\b","\\bFOXC1\\b","\\bGATA3\\b","\\bKIT\\b"),collapse="|") #shortlist other
#	geneList = paste(c("\\bLYZ\\b","\\bATOH1\\b","\\bMIXL1\\b","\\bMMP7\\b","\\bMUC2\\b","\\bPLA2G2A\\b","\\bSOX9\\b","\\bT\\b","\\bWNT3A\\b"),collapse="|") #shortlist4
#	geneList = paste(c("\\bHNF4A\\b","\\bONECUT1\\b","\\bONECUT2\\b","\\bHOXA1\\b","\\bHOXA2\\b","\\bHOXC5\\b","\\bLRIG1\\b"),collapse="|") #shortlist3
#	geneList = paste(c("\\bCER1\\b","\\bGREM2\\b","\\bCXCR4\\b","\\bEPCAM\\b","\\bFOXA2\\b"),collapse="|") #shortlist2
#	geneList = paste(c("\\bSERPINA1\\b","\\bAFP\\b","\\bALB\\b","\\bASCL2\\b","\\bATOH1\\b","\\bAXIN2\\b","\\bCD44\\b"),collapse="|") #shortlist1
#	geneList = paste(c("\\bTDX1\\b","\\bSox2\\b","\\bCD44\\b","\\bPDX1\\b","\\bCDX1\\b","\\bCDX2\\b","\\bHOXC5\\b","\\bT\\b","\\bMixl1\\b","\\bCXCR4\\b","\\bHHEX\\b","\\bFOXA2\\b","\\bCERB\\b","\\bSOX17\\b","\\bHOXA1\\b","\\bGATA4\\b","\\bHNF4a\\b","\\bEpcam\\b","\\bHOXA2\\b","\\bAFP\\b","\\bALB\\b","\\bA1AT\\b","\\bHNF4ALPHA\\b","\\bCK18\\b","\\bHNF6\\b","\\bHLXB9\\b","\\bNGN3\\b"),collapse="|") # Patterning
#	geneList = paste(c("\\bADAM19\\b","\\bAPOA1\\b","\\bAPOE\\b","\\bBAMBI\\b","\\bBMP2\\b","\\bBMP7\\b","\\bCDKN1C\\b","\\bCER1\\b","\\bCOL5A2\\b","\\bCOL4A5\\b","\\bCOL9A2\\b","\\bCRIP1\\b","\\bCXCR4\\b","\\bDAB2\\b","\\bDKK1\\b","\\bDKK3\\b","\\bEOMES\\b","\\bEPHA2\\b","\\bESRRG\\b","\\bEYA1\\b","\\bFKBP9\\b","\\bFOXA1\\b","\\bFOXA2\\b","\\bFOXC1\\b","\\bFZD4\\b","\\bGAD1\\b","\\bGATA3\\b","\\bGATA4\\b","\\bGATA6\\b","\\bGLIS3\\b","\\bGLUD1\\b","\\bGLUD2\\b","\\bGSC\\b","\\bH2AFY2\\b","\\bHHEX\\b","\\bHSZFP36\\b","\\bID1\\b","\\bID3\\b","\\bIGF2\\b","\\bIGFBP3\\b","\\bKIT\\b","\\bKRT19\\b","\\bKRT8\\b","\\bLEFTY1\\b","\\bLEFTY2\\b","\\bLHX1\\b","\\bMANEA\\b","\\bMANEAL\\b","\\bMIXL1\\b","\\bMSX2\\b","\\bNID2\\b","\\bNODAL\\b","\\bNRP1\\b","\\bOTX2\\b","\\bPDZK1\\b","\\bPERP\\b","\\bPPOX\\b","\\bSOX17\\b","\\bSYTL5\\b","\\bTBC1D9\\b","\\bTNNC1\\b","\\bTYRO3\\b"),collapse="|") # Spence Wells ST1b
#	geneList = paste(c("\\bAXIN2\\b","\\bDUOX2\\b","\\bKLF5\\b","\\bKRT20\\b","\\bLGR5\\b","\\bMEP1A\\b","\\bMGAM\\b","\\bMUC13\\b","\\bMUC2\\b","\\bOLFM4\\b","\\bSLC15A1\\b","\\bSOX9\\b","\\bTERT\\b","\\bVIL1\\b","\\bVILL\\b","\\bREG4\\b","\\bMUC17\\b","\\bSLC15A1\\b","\\bEPCAM\\b","\\bGUCA2A\\b","\\bCLCA4\\b","\\bTFF2\\b"),collapse="|")

	rows = grep(geneList,featureNames(ensembl_vsd),value=FALSE)
	cols_cat <- c(19,25,1,22,24,18,3,23,14,5,10,2,30,21) #categorized grouping of org and diff samples ()
	cols_pair <- c(19,10,25,2,14,30,5,21) #pairwise grouping of org and diff samples ()
#	selected_vsd = exprs(ensembl_vsd)[rows,cols]
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#	norm_vsd = 100*selected_vsd/rowSums(selected_vsd)

<<<<<<< HEAD
	pdf('r_heatmap_noa_definitiveendoderm2.pdf',width=10, height=10)
=======
	pdf('r_heatmap_noa_alllist.pdf',width=10, height=10)
>>>>>>> origin/master
	heatmap.2(exprs(ensembl_vsd)[rows,],col = hmcol, trace="none", margin=c(15, 6))
#	heatmap.2(selected_vsd, dendrogram = c("none"),col = hmcol, trace="none", margin=c(15, 6), Colv = F)
	dev.off()
	dev.new(height=10)
	heatmap.2(exprs(ensembl_vsd)[rows,],col = hmcol, trace="none", margin=c(15, 6))
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
<<<<<<< HEAD
	pdf('r_PCA_noa_withcontrols.pdf')
	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=42859))
	dev.off()
	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=42859)) # 42859 genes, 120000 transcripts
=======
#	pdf('r_PCA_genenames.pdf')
	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=42859))
#	dev.off()
#	print(plotPCA(ensembl_vsd, intgroup=c('condition'), ntop=120000)) # 42859 genes, 120000 transcripts
>>>>>>> origin/master
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

#dataPrep()

#binomDist()
#heatmapPlot()
#similarityPlot()
#deseqPCA()
#allCombos()
#screePlot()
#rPCA()
#rFactorAnalysis()
#calcProb()

heatmapGenes('aquaporins',["AQP1","AQP2","AQP3","AQP4","AQP5","AQP6","AQP7","AQP8","AQP9","AQP10","AQP11"])
heatmapGenes('maturity',['AXIN2','CDH1','CDH17','CDX2','DEFA5','DEF6','DUOX2','ELF3','EPCAM','KLF4','KLF5','LGR5','LYZ','MUC13','MUC2','SOX9','VIL1','SI','MMP7','PLA2G2A','WNT3A','ATOH1','TFF3','TERT'])
heatmapGenes('adultstem',['LGR5','SOX9','KLF4','KLF5','OLFM4','TERT'])
heatmapGenes('pattern',['BARX1','SFRP1','HHEX','SOX2','DKK1','FABP2','CDX1','CDX2','HOXC5'])

