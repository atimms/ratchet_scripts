library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kim_rnaseq_hb_1018";
setwd(workingDir);

###all samples
##read in count and metadata
countData1 <- read.table('kim_rnaseq_hb_1018.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_hb_1018.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~brain+ sample_type)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="visal_0618_striatum.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_hb_1018.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$brain, rld$sample_type ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_hb_1018.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("brain"))
ggsave('kim_rnaseq_hb_1018.brain_pca.pdf') 
plotPCA(rld, intgroup = c("sample_type"))
ggsave('kim_rnaseq_hb_1018.sample_type_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("brain", "sample_type")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='kim_rnaseq_hb_1018.25_var_gene_clustering.pdf', width = 7, height = 5)








