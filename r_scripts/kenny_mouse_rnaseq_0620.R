##just use common library, so check which one to use and then set that parameter
.libPaths()
.libPaths( .libPaths()[2] )
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kenny_mouse_0620";
setwd(workingDir);

###eucomm_0620
##read in count and metadata
countData1 <- read.table('eucomm_0620.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('eucomm_0620.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="eucomm_0620.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="eucomm_0620.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$genotype, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='eucomm_0620.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('eucomm_0620.genotype_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sample_name"))
ggsave('eucomm_0620.sample_name_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='eucomm_0620.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression --- just ko vs wt
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "ko", "wt"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:17776,]
write.csv(resOrdered2DF, file="eucomm_0620.ko_vs_wt.csv")

##plots for genes of interest
plotCounts(dds, gene="Sox8", intgroup="genotype")
dev.copy2pdf(file='eucomm_0620.Sox8_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Tsr3", intgroup="genotype")
dev.copy2pdf(file='eucomm_0620.Tsr3_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Ddah2", intgroup="genotype")
dev.copy2pdf(file='eucomm_0620.Ddah2_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Tubb5", intgroup="genotype")
dev.copy2pdf(file='eucomm_0620.Tubb5_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Zfp760", intgroup="genotype")
dev.copy2pdf(file='eucomm_0620.Zfp760_counts.pdf', width = 7, height = 5)



###wegner_0620
##read in count and metadata
countData1 <- read.table('wegner_0620.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('wegner_0620.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="wegner_0620.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="wegner_0620.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$genotype, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='wegner_0620.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('wegner_0620.genotype_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sample_name"))
ggsave('wegner_0620.sample_name_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='wegner_0620.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression --- just ko vs wt
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "ko", "wt"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:17071,]
write.csv(resOrdered2DF, file="wegner_0620.ko_vs_wt.csv")

##plots for genes of interest
plotCounts(dds, gene="Sox8", intgroup="genotype")
dev.copy2pdf(file='wegner_0620.Sox8_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Tsr3", intgroup="genotype")
dev.copy2pdf(file='wegner_0620.Tsr3_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Ddah2", intgroup="genotype")
dev.copy2pdf(file='wegner_0620.Ddah2_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Tubb5", intgroup="genotype")
dev.copy2pdf(file='wegner_0620.Tubb5_counts.pdf', width = 7, height = 5)
plotCounts(dds, gene="Zfp760", intgroup="genotype")
dev.copy2pdf(file='wegner_0620.Zfp760_counts.pdf', width = 7, height = 5)
