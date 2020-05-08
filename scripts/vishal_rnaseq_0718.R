library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/vishal_rnaseq_0618";
setwd(workingDir);

###striatum
##read in count and metadata
countData1 <- read.table('vishal_0618_striatum.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_0618_striatum.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~animal + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="visal_0618_striatum.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="visal_0618_striatum.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$animal, rld$treatment ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='visal_0618_striatum.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("treatment"))
ggsave('visal_0618_striatum.treatment_pca.pdf') 
plotPCA(rld, intgroup = c("animal"))
ggsave('visal_0618_striatum.animal_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("animal", "treatment")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='visal_0618_striatum.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "sham", "CPB_DHCA"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_striatum.sham_vs_CPB_DHCA.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "sham", "CPB_DHCA_G_CSF"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_striatum.sham_vs_CPB_DHCA_G_CSF.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "CPB_DHCA", "CPB_DHCA_G_CSF"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_striatum.CPB_DHCA_vs_CPB_DHCA_G_CSF.csv")



###hippocampus
##read in count and metadata
countData1 <- read.table('vishal_0618_hippocampus.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_0618_hippocampus.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~animal + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="visal_0618_hippocampus.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="visal_0618_hippocampus.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$animal, rld$treatment ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='visal_0618_hippocampus.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("treatment"))
ggsave('visal_0618_hippocampus.treatment_pca.pdf') 
plotPCA(rld, intgroup = c("animal"))
ggsave('visal_0618_hippocampus.animal_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("animal", "treatment")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='visal_0618_hippocampus.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "sham", "CPB_DHCA"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_hippocampus.sham_vs_CPB_DHCA.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "sham", "CPB_DHCAG_CSF"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_hippocampus.sham_vsCPB_DHCAG_CSF.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "CPB_DHCA", "CPB_DHCAG_CSF"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_hippocampus.CPB_DHCA_vsCPB_DHCAG_CSF.csv")



###14samples
##read in count and metadata
countData1 <- read.table('vishal_0618_14samples.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_0618_14samples.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~animal + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="visal_0618_14samples.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="visal_0618_14samples.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$animal, rld$treatment ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='visal_0618_14samples.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("treatment"))
ggsave('visal_0618_14samples.treatment_pca.pdf') 
plotPCA(rld, intgroup = c("animal"))
ggsave('visal_0618_14samples.animal_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("animal", "treatment")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='visal_0618_14samples.25_var.gene_clustering.pdf', width = 7, height = 5)



###striatum_reduced removed one outlier
##read in count and metadata
countData1 <- read.table('vishal_0618_striatum_reduced.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_0618_striatum_reduced.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~animal + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="visal_0618_striatum_reduced.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="visal_0618_striatum_reduced.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$animal, rld$treatment ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='visal_0618_striatum_reduced.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("treatment"))
ggsave('visal_0618_striatum_reduced.treatment_pca.pdf') 
plotPCA(rld, intgroup = c("animal"))
ggsave('visal_0618_striatum_reduced.animal_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("animal", "treatment")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='visal_0618_striatum_reduced.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "sham", "CPB_DHCA"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_striatum_reduced.sham_vs_CPB_DHCA.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "sham", "CPB_DHCA_G_CSF"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_striatum_reduced.sham_vs_CPB_DHCA_G_CSF.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("treatment", "CPB_DHCA", "CPB_DHCA_G_CSF"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21223,]
write.csv(resOrdered2DF, file="visal_0618_striatum_reduced.CPB_DHCA_vs_CPB_DHCA_G_CSF.csv")


###striatum_nogcsf remove gcsf treatments
##read in count and metadata
countData1 <- read.table('vishal_0618_striatum_nogcsf.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_0618_striatum_nogcsf.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~animal + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="vishal_0618_striatum_nogcsf.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="vishal_0618_striatum_nogcsf.deseq.rlog_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("treatment")) + theme_bw()
ggsave('vishal_0618_striatum_nogcsf.treatment_pca.png')
plotPCA(rld, intgroup = c("animal")) + theme_bw()
ggsave('vishal_0618_striatum_nogcsf.animal_pca.png') 


