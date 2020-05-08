library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/kim_hipscs_0318";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('kim_hipacs_0318_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_hipacs_0318_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~Sample + Source + Genotype)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
colData(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="kim_hipacs_0318.norm_counts.csv")

##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_hipacs_0318.deseq.rlog_counts.csv")

##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Sample, rld$Source, rld$Genotype, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_hipacs_0318.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
#plotPCA(rld, intgroup = c("tissue", "sex", "rnaaccess_batch"))
plotPCA(rld, intgroup = "Sample")
ggsave('kim_hipacs_0318.sample_pca.Sample.pdf') 
plotPCA(rld, intgroup = "Source")
ggsave('kim_hipacs_0318.sample_pca.Source.pdf') 
plotPCA(rld, intgroup = "Genotype")
ggsave('kim_hipacs_0318.sample_pca.Genotype.pdf') 

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("Sample","Source","Genotype")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='kim_hipacs_0318.deseq.gene_clustering.pdf', width = 7, height = 5)

##differeential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))


##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77", "78"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77_vs_78.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77", "78"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77_vs_78.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77", "77null"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77_vs_77null.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77", "SAH01px"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77_vs_SAH01px.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77", "SAH02"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77_vs_SAH02.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77null", "SAH02"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77null_vs_SAH02.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77null", "SAH01px"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77null_vs_SAH01px.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "77null", "78"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.77null_vs_78.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "78", "SAH01px"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.78_vs_SAH01px.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "78", "SAH02"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.78_vs_SAH02.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("Sample", "SAH01px", "SAH02"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20208,]
write.csv(resOrdered2DF, file="kim_hipacs_0318.de.SAH01px_vs_SAH02.csv")


####now the comparisons withe geo dataset and rnaseq 0218
##read in count and metadata
countData1 <- read.table('kim_hipacs_0318_comparison.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_hipacs_0318_comparison.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
colData(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="kim_hipacs_0318.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
#head(assay(rld), 3)
#head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_hipacs_0318.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$id, rld$tissue, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_hipacs_0318_comparison.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
#plotPCA(rld, intgroup = c("tissue", "sex", "rnaaccess_batch"))
plotPCA(rld, intgroup = "id")
ggsave('kim_hipacs_0318_comparison.sample_pca.id.pdf') 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_hipacs_0318_comparison.sample_pca.tissue.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("id","tissue")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='kim_hipacs_0318_comparison.deseq.gene_clustering.pdf', width = 7, height = 5)


##divide data to get rlog counts from encode data
##read in count and metadata
countData1 <- read.table('GSE78474_0318.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('GSE78474_0318.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~1 )
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
colData(dds)
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="GSE78474_encode_pcl_0318.deseq.rlog_counts.csv")

