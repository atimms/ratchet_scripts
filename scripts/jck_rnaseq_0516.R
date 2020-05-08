library(WGCNA);
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/home/atimms/jck_rnaseq";
setwd(workingDir);

##all data

##read in count and metadata
#countData1 <- read.table('jck_rnaseq_0516.htseq_counts.txt', header=T, row.names=1)
#colData1 <- read.table('jck_rnaseq_0516.metadata.txt', header=T, row.names=1)
countData1 <- read.table('jck_rnaseq_0516.htseq_counts.0916.txt', header=T, row.names=1)
colData1 <- read.table('jck_rnaseq_0516.metadata.0916.txt', header=T, row.names=1)
head(countData1)
head(colData1)


##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ age + mutant)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
count_data <- counts(dds, normalized=TRUE)
write.csv(count_data, file="jck_rnaseq_0516.deseq.norm_counts.0916.csv")

##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="jck_rnaseq_0516.deseq.rlog_counts.0916.csv")

##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$age, rld$mutant, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='jck.sample_heatmap.0916.pdf', width = 7, height = 5)

##principal components analysis 
plotPCA(rld, intgroup = c("age", "mutant"))
dev.copy2pdf(file='jck.sample_pca.0916.pdf', width = 7, height = 5)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("age","mutant")])
pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='jck.deseq.gene_clustering.0916.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
summary(res)

##to get a specific test i.e. jck vs wt
res2 <- results(dds, contrast=c("mutant", "jck", "wt"))
##get summary
summary(res2) #lots of significant genes

##subset genes to get most significant and sort by fold change
resSig <- subset(res2, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

##plot all points for individual genes
##using inbuilt function
#plotCounts(dds, gene="Uroc1", intgroup=c("mutant"))
#dev.copy2pdf(file='jck.Uroc1.pdf', width = 7, height = 5)
#plotCounts(dds, gene="Tacstd2", intgroup=c("mutant"))
#dev.copy2pdf(file='jck.Tacstd2.pdf', width = 7, height = 5)
#plotCounts(dds, gene="Mal2", intgroup=c("mutant"))
#dev.copy2pdf(file='jck.Mal2.pdf', width = 7, height = 5)
               
##save differentially expression results
##sort results by adjusted p-value
res2Ordered <- res2[order(res2$padj),]
head(res2Ordered)
##save results as dataframe and take top 1k results, then write csv file
resOrderedDF <- as.data.frame(res2Ordered)[1:20000,]
write.csv(resOrderedDF, file="jck.results.jck_vs_wt.0916.csv")

##p5
#workingDir = "/home/atimms/jck_rnaseq/p5";
#setwd(workingDir);
##read in count and metadata
countData1 <- read.table('jck_rnaseq_0516.p5.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('jck_rnaseq_0516.p5.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ mutant)
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
##differential expression
##do the test
dds <- DESeq(dds)
res <- results(dds, contrast=c("mutant", "jck", "wt"))
summary(res)
##save differentially expression results
##sort results by adjusted p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)
##save results as dataframe and take top 1k results, then write csv file
resOrderedDF <- as.data.frame(resOrdered)[1:20000,]
write.csv(resOrderedDF, file="jck.results.jck_vs_wt.p5.0916.csv")

##p10
#workingDir = "/home/atimms/jck_rnaseq/p10";
#setwd(workingDir);
##read in count and metadata
countData1 <- read.table('jck_rnaseq_0516.p10.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('jck_rnaseq_0516.p10.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ mutant)
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
##differential expression
##do the test
dds <- DESeq(dds)
res <- results(dds, contrast=c("mutant", "jck", "wt"))
summary(res)
##save differentially expression results
##sort results by adjusted p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)
##save results as dataframe and take top 1k results, then write csv file
resOrderedDF <- as.data.frame(resOrdered)[1:20000,]
write.csv(resOrderedDF, file="jck.results.jck_vs_wt.p10.0916.csv")

##p15
#workingDir = "/home/atimms/jck_rnaseq/p15";
#setwd(workingDir);
##read in count and metadata
countData1 <- read.table('jck_rnaseq_0516.p15.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('jck_rnaseq_0516.p15.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ mutant)
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
##differential expression
##do the test
dds <- DESeq(dds)
res <- results(dds, contrast=c("mutant", "jck", "wt"))
summary(res)
##save differentially expression results
##sort results by adjusted p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)
##save results as dataframe and take top 1k results, then write csv file
resOrderedDF <- as.data.frame(resOrdered)[1:20000,]
write.csv(resOrderedDF, file="jck.results.jck_vs_wt.p15.0916.csv")

##p20
#workingDir = "/home/atimms/jck_rnaseq/p20";
#setwd(workingDir);
##read in count and metadata
countData1 <- read.table('jck_rnaseq_0516.p20.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('jck_rnaseq_0516.p20.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ mutant)
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
##differential expression
##do the test
dds <- DESeq(dds)
res <- results(dds, contrast=c("mutant", "jck", "wt"))
summary(res)
##save differentially expression results
##sort results by adjusted p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)
##save results as dataframe and take top 1k results, then write csv file
resOrderedDF <- as.data.frame(resOrdered)[1:20000,]
write.csv(resOrderedDF, file="jck.results.jck_vs_wt.p20.0916.csv")