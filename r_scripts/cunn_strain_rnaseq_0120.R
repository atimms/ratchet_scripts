library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/cunn_strain_rnaseq_0220";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sample_only + strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="strain_rnaseq_0220.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="strain_rnaseq_0220.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$sample_only, rld$strain ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='strain_rnaseq_0220.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('strain_rnaseq_0220.sample_pca.pdf')
plotPCA(rld, intgroup = c("strain"))
ggsave('strain_rnaseq_0220.strain_pca.pdf')
plotPCA(rld, intgroup = c("sex"))
ggsave('strain_rnaseq_0220.sex_pca.pdf')
plotPCA(rld, intgroup = c("gene"))
ggsave('strain_rnaseq_0220.gene_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_only", "sex", "strain", "gene")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=2)
dev.copy2pdf(file='strain_rnaseq_0220.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
dds
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.all_samples.strain0_vs_strain100.csv")

##to get a specific test i.e. 0 vs 100 strain in individual samples
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sample_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##redo DE
dds <- DESeq(dds)
dds
##get specific tests
res2 <- results(dds, contrast=c("sample_strain", "1017_0", "1017_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.1017.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "3035_0", "3035_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.3035.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C1653_0", "C1653_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1653.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C1671_0", "C1671_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1671.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C1856_0", "C1856_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1856.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C1859_0", "C1859_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1859.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C1926_0", "C1926_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1926.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C2038_0", "C2038_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C2038.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C2084_0", "C2084_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C2084.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "3007_0", "3007_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.3007.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C3066_0", "C3066_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C3066.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "STL19_0", "STL19_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.STL19.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "1032_0", "1032_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.1032.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "4025_0", "4025_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.4025.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "4032_0", "4032_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.4032.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "AUT25_0", "AUT25_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.AUT25.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C1625_0", "C1625_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1625.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "C3049_0", "C3049_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C3049.strain0_vs_strain100.csv")
res2 <- results(dds, contrast=c("sample_strain", "OST85_0", "OST85_100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23179,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.OST85.strain0_vs_strain100.csv")





#split the data into individual samples and go from there......
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_1017.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_1017.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18583,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.1017_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_3035.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_3035.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18724,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.3035_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C1653.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C1653.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19266,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1653_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C1671.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C1671.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19355,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1671_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C1856.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C1856.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18946,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1856_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C1859.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C1859.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18871,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1859_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C1926.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C1926.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19096,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1926_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C1957.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C1957.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18479,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1957_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C2038.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C2038.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18800,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C2038_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C2084.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C2084.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18622,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C2084_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_3007.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_3007.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18388,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.3007_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C3066.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C3066.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19012,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C3066_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_STL19.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_STL19.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18893,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.STL19_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_1032.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_1032.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19958,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.1032_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_4025.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_4025.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19797,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.4025_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_4032.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_4032.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19342,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.4032_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_AUT25.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_AUT25.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20189,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.AUT25_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C1625.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C1625.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19406,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C1625_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_C3049.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_C3049.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19766,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.C3049_ind.strain0_vs_strain100.csv")
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_OST85.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_OST85.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19675,]
write.csv(resOrdered2DF, file="strain_rnaseq_0220.OST85_ind.strain0_vs_strain100.csv")

