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
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~Sample)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
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



