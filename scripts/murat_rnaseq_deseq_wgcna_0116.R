library(WGCNA);
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/home/atimms/murat_rnaseq/all_tissues";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('murat_etoh_0116.clustering.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('murat_etoh_0116.clustering.metadata.txt', header=T, row.names=1)

head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ tissue + individual)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="murat_etoh_0116.deseq.rlog_counts.csv")

##compare log transformation with rlog - not really needed
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)

##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$individual, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='murat_etoh_0116.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
plotPCA(rld, intgroup = c("tissue", "individual"))
ggsave('murat_etoh_0116.sample_pca.pdf') 

##multidimensional scaling plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=tissue,shape=individual)) + geom_point(size=3)
ggsave('murat_etoh_0116.sample_mds.pdf') 

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","individual")])
pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='murat_etoh_0116.deseq.gene_clustering.pdf', width = 7, height = 5)

##differeential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
summary(res)

##to get a specific test i.e. pk vs gr
res2 <- results(dds, contrast=c("tissue", "gr", "pk"))
##get summary
summary(res2) #lots of significant genes

##subset genes to get most significant and sort by fold change
resSig <- subset(res2, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

##plot all points for individual genes
##using inbuilt function
plotCounts(dds, gene="CALB1", intgroup=c("tissue"))
##better with ggplot
data <- plotCounts(dds, gene="GRIA1", intgroup=c("tissue","individual"), returnData=TRUE)
ggplot(data, aes(x=tissue, y=count,color=tissue,shape=individual)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=3)
ggsave('dwm0116.GRIA1.tissue.pdf') 

##save differentially expression results
##sort results by adjusted p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)
##save results as dataframe and take top 10k results, then write csv file
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF, file="results.etoh_vs_ctl.csv")


##wgcna from rlog transformed values
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in the data set
dataset_1 = read.csv("dwm0116.deseq.rlog_counts.csv", header=T);

#transpose
dataset_2 = as.data.frame(t(dataset_1[, -c(1)]));
# Set the column headings
names(dataset_2) = dataset_1$X;
rownames(dataset_2) = names(dataset_1)[-c(1)];


##cluster samples to look for outliers
sampleTree = hclust(dist(dataset_2), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.copy2pdf(file="dwm0116.wgcna.sample_clustering.pdf", width = 7, height = 5)

##read in trait data
traitData <- read.table('dwm0116.clustering.metadata.txt', header=T);
dim(traitData)
names(traitData)
head(traitData)
# Form a data frame analogous to expression data that will hold the clinical traits.


# Re-cluster samples
sampleTree2 = hclust(dist(dataset_2), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

