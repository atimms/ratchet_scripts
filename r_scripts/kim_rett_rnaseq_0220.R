library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kim_rett_rnaseq_0119";
setwd(workingDir);

###all samples
##read in count and metadata
countData1 <- read.table('kim_rett_rnaseq_0119.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rett_rnaseq_0119_bias.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
#dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~region + dx + region_dx)
#dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~bias + region_dx)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~bias_3 + region + dx)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="kim_rett_rnaseq_0119_bias.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_rett_rnaseq_0119.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$region, rld$dx, rld$bias_3 ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rett_rnaseq_0119_bias.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("region"))
ggsave('kim_rett_rnaseq_0119_bias.region_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("dx"))
ggsave('kim_rett_rnaseq_0119_bias.dx_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("bias_3"))
ggsave('kim_rett_rnaseq_0119_bias.bias_pca.pdf', width=6, height = 6)
##get just data
pcaData_bias <- plotPCA(rld, intgroup = c("bias_3"), returnData = TRUE)
pcaData_bias
write.csv(pcaData_bias,file='kim_rett_rnaseq_0119_bias.bias_pca_data.csv')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('kim_rett_rnaseq_0119_bias.sample_pca.pdf', width=6, height = 6)
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("region", "dx", "bias_3")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='kim_rett_rnaseq_0119_bias.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression --- just rett vs control
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("dx", "Rett", "control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19623,]
write.csv(resOrdered2DF, file="kim_rett_rnaseq_0119_bias.Rett_vs_control.csv")

#redo, so splitting the data into region and dx
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~bias_3 + region_dx)
dds
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)
##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("region_dx", "temporal_Rett", "temporal_control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19623,]
write.csv(resOrdered2DF, file="kim_rett_rnaseq_0119_bias.temporal.Rett_vs_control.csv")
res2 <- results(dds, contrast=c("region_dx", "cingulate_Rett", "cingulate_control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19623,]
write.csv(resOrdered2DF, file="kim_rett_rnaseq_0119_bias.cingulate.Rett_vs_control.csv")

###for a slightly different gtf which differtiates 
##read in count and metadata
countData1 <- read.table('kim_rett_rnaseq_mecp1_0319.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rett_rnaseq_mecp1_0319_bias.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ bias_3 + region + dx)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="kim_rett_rnaseq_mecp1_0319_bias.norm_counts.csv")
##differential expression --- just rett vs control
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("dx", "Rett", "control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20413,]
write.csv(resOrdered2DF, file="kim_rett_rnaseq_mecp1_0319_bias.Rett_vs_control.csv")

