library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kevin_rpe_library_0619";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('kevin_rpe_ribo_depleted_0619.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kevin_rpe_ribo_depleted_0619.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##using all samples
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~status + sex + rs1742162)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="kevin_rpe_ribo_depleted_0619.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kevin_rpe_ribo_depleted_0619.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$rs1742162, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kevin_rpe_ribo_depleted_0619.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("status"))
ggsave('kevin_rpe_ribo_depleted_0619.status_pca.pdf')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('kevin_rpe_ribo_depleted_0619.sample_pca.pdf')
plotPCA(rld, intgroup = c("sex"))
ggsave('kevin_rpe_ribo_depleted_0619.sex_pca.pdf')
plotPCA(rld, intgroup = c("rs1742162"))
ggsave('kevin_rpe_ribo_depleted_0619.rs1742162_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("status", "sex", "rs1742162")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
#pheatmap(mat, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='kevin_rpe_ribo_depleted_0619.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("rs1742162", "het", "wt"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:26727,]
write.csv(resOrdered2DF, file="kevin_rpe_ribo_depleted_0619.het_vs_all_wt.csv")

