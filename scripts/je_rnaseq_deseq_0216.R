library(WGCNA);
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/home/atimms/je_rnaseq";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('je_rnaseq_0216.clustering.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('je_rnaseq_0216.clustering.metadata.txt', header=T, row.names=1)

head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ experiment + mutation)
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
write.csv(assay(rld), file="je_rnaseq_0216.deseq.rlog_counts.csv")

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
rownames(sampleDistMatrix) <- paste( rld$experiment, rld$mutation, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='je_rnaseq_0216.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
plotPCA(rld, intgroup = c("experiment", "mutation"))
ggsave('je_rnaseq_0216.sample_pca.pdf') 

##multidimensional scaling plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=experiment,shape=mutation)) + geom_point(size=3)
ggsave('je_rnaseq_0216.sample_mds.pdf') 

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("experiment","mutation")])
pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='je_rnaseq_0216.deseq.gene_clustering.pdf', width = 7, height = 5)

##differential expression -- for both together
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
summary(res)

##subset genes to get most significant and sort by fold change
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
#head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

##save differentially expression results
##sort results by adjusted p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)
##save results as dataframe and take top x results, then write csv file
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF, file="je_rnaseq_0216.results.combined_mut.csv")

##look at data individually
##jessica
##read in count and metadata
countData1 <- read.table('je_rnaseq_0216.jess.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('je_rnaseq_0216.jess_metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~mutation)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##rlog transform
rld <- rlogTransformation(dds)
##check
head(assay(rld), 3)
head(assay(dds),3)
##get rlog fold changes
res <- data.frame(
  assay(rld), 
  avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
  rLogFC_a = assay(rld)[,2] - assay(rld)[,1],
  rLogFC_b = assay(rld)[,1] - assay(rld)[,2])
over_res <- ( res[ order(res$rLogFC_a),] )
under_res <- ( res[ order(res$rLogFC_b),] )

##save differentially expression results
##save results as dataframe and take top x results, then write csv file
over_resDF <- as.data.frame(over_res)[1:100,]
under_resDF <- as.data.frame(under_res)[1:100,]
write.csv(over_resDF, file="je_rnaseq_0216.jessica.results.overexpressed.csv")
write.csv(under_resDF, file="je_rnaseq_0216.jessica.results.underexpressed.csv")
##save all genes
over_resDF <- as.data.frame(over_res)[1:17108,]
write.csv(over_resDF, file="je_rnaseq_0216.jessica.results.all_genes.csv")

##esra
##read in count and metadata
countData1 <- read.table('je_rnaseq_0216.esra.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('je_rnaseq_0216.esra_metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~mutation)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##rlog transform
rld <- rlogTransformation(dds)
##check
head(assay(rld), 3)
head(assay(dds),3)

##compare log transformation with rlog - not really needed
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)


##get rlog fold changes
res <- data.frame(
  assay(rld), 
  avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
  rLogFC_a = assay(rld)[,2] - assay(rld)[,1],
  rLogFC_b = assay(rld)[,1] - assay(rld)[,2])
head( res[ order(res$rLogFC_a),] )
head( res[ order(res$rLogFC_b),] )
over_res <- ( res[ order(res$rLogFC_a),] )
under_res <- ( res[ order(res$rLogFC_b),] )

##save differentially expression results
##save results as dataframe and take top x results, then write csv file
over_resDF <- as.data.frame(over_res)[1:100,]
under_resDF <- as.data.frame(under_res)[1:100,]
write.csv(over_resDF, file="je_rnaseq_0216.esra.results.overexpressed.all_genes.csv")
write.csv(under_resDF, file="je_rnaseq_0216.esra.results.underexpressed.csv")
##save all genes
over_resDF <- as.data.frame(over_res)[1:18000,]
write.csv(over_resDF, file="je_rnaseq_0216.esra.results.all_genes.csv")
##esra-- using jessica wt as control
##read in count and metadata
countData1 <- read.table('je_rnaseq_0216.esra_alt.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('je_rnaseq_0216.esra_alt.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~mutation)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##rlog transform
rld <- rlogTransformation(dds)
##check
head(assay(rld), 3)
head(assay(dds),3)
##get rlog fold changes
res <- data.frame(
  assay(rld), 
  avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
  rLogFC_a = assay(rld)[,2] - assay(rld)[,1],
  rLogFC_b = assay(rld)[,1] - assay(rld)[,2])
head( res[ order(res$rLogFC_a),] )
head( res[ order(res$rLogFC_b),] )
over_res <- ( res[ order(res$rLogFC_a),] )
under_res <- ( res[ order(res$rLogFC_b),] )

##save differentially expression results
##save results as dataframe and take top x results, then write csv file
over_resDF <- as.data.frame(over_res)[1:100,]
under_resDF <- as.data.frame(under_res)[1:100,]
write.csv(over_resDF, file="je_rnaseq_0216.esra_alt.results.overexpressed.csv")
write.csv(under_resDF, file="je_rnaseq_0216.esra_alt.results.underexpressed.csv")
##save all genes
over_resDF <- as.data.frame(over_res)[1:18000,]
write.csv(over_resDF, file="je_rnaseq_0216.esra_alt.results.all_genes.csv")