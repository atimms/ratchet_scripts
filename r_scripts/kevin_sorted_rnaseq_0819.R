library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kevin_rnaseq_0819";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('sorted_cell_0819_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sorted_cell_0819_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##using all 'timepoints'
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ clone + reporter)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="sorted_cell_0819_all.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sorted_cell_0819_all.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$clone, rld$reporter ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='sorted_cell_0819_all.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("clone"))
ggsave('sorted_cell_0819_all.clone_pca.pdf')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('sorted_cell_0819_all.sample_pca.pdf')
plotPCA(rld, intgroup = c("reporter"))
ggsave('sorted_cell_0819_all.reporter_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("clone", "reporter")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='sorted_cell_0819_all.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("reporter", "CRALBP-GFP", "rhoK-tdT"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:24422,]
write.csv(resOrdered2DF, file="sorted_cell_0819_all.CRALBP-GFP_vs_rhoK-tdT.csv")


##just look at individual clones
####### c4 #####
##read in count and metadata
countData1 <- read.table('sorted_cell_0819_c4.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sorted_cell_0819_c4.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ reporter)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##principal components analysis 
plotPCA(rld, intgroup = c("reporter"))
ggsave('sorted_cell_0819_c4.reporter_pca.pdf')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('sorted_cell_0819_c4.sample_pca.pdf')
##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("reporter", "CRALBP-GFP", "rhoK-tdT"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22789,]
write.csv(resOrdered2DF, file="sorted_cell_0819_c4.CRALBP-GFP_vs_rhoK-tdT.csv")

####### c5 #####
##read in count and metadata
countData1 <- read.table('sorted_cell_0819_c5.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sorted_cell_0819_c5.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ reporter)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##principal components analysis 
plotPCA(rld, intgroup = c("reporter"))
ggsave('sorted_cell_0819_c5.reporter_pca.pdf')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('sorted_cell_0819_c5.sample_pca.pdf')
##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("reporter", "CRALBP-GFP", "rhoK-tdT"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:23546,]
write.csv(resOrdered2DF, file="sorted_cell_0819_c5.CRALBP-GFP_vs_rhoK-tdT.csv")


