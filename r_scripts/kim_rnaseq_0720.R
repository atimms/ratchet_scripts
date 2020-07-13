##just use common library, so check which one to use and then set that parameter
.libPaths()
.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/data/atimms/kim_rnaseq_0720";
setwd(workingDir);

###kim_rnaseq_0720_all
##read in count and metadata
countData1 <- read.table('kim_rnaseq_0720_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0720_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Age_GW <- as.factor(colData1$Age_GW)
colData1$Batch <- as.factor(colData1$Batch)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Batch + Region + Group)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="kim_rnaseq_0720_all.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_rnaseq_0720_all.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Region, rld$Group, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0720_all.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("Specimen"))
ggsave('kim_rnaseq_0720_all.Specimen_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sample_name"))
ggsave('kim_rnaseq_0720_all.sample_name_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Source"))
ggsave('kim_rnaseq_0720_all.Source_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Group"))
ggsave('kim_rnaseq_0720_all.Group_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('kim_rnaseq_0720_all.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Age_GW"))
ggsave('kim_rnaseq_0720_all.Age_GW_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sequencing"))
ggsave('kim_rnaseq_0720_all.Sequencing_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Batch"))
ggsave('kim_rnaseq_0720_all.Batch_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Region"))
ggsave('kim_rnaseq_0720_all.Region_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Region", "Group", "Specimen", "Source", "Sex", "Age_GW", "Sequencing", "Batch")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='kim_rnaseq_0720_all.25_var_gene_clustering.pdf', width = 10, height = 10)
