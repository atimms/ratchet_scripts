library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/kim_rnaseq_0517";
setwd(workingDir);


###4 tissues#########
##read in count and metadata
countData1 <- read.table('kim_rnaseq_1016.star_fc.4_tissues.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_1016.star_fc.4_tissues.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~sex + rnaaccess_batch + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="kim_rnaseq_1016.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_rnaseq_1016.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$individual, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0517.pcl_egl_bulk_rl.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0517.pcl_egl_bulk_rl.sample_pca.tissue_only.pdf') 



###3 tissues#########
##read in count and metadata
countData1 <- read.table('kim_rnaseq_1016.star_fc.3_tissues.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_1016.star_fc.3_tissues.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~sex + rnaaccess_batch + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="kim_rnaseq_1016.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_rnaseq_1016.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$individual, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0517.pcl_egl_bulk.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0517.pcl_egl_bulk.sample_pca.tissue_only.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "individual")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='kim_rnaseq_0517.pcl_egl_bulk.25vargene_clustering.pdf', width = 7, height = 5)
##take genes from kim's list
####from gene list
#genes <- read.csv('ASD_162gene_list.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- read.csv('kim_gl_061417.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("tissue", "individual")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
#dev.copy2pdf(file='kim_rnaseq_1016.deseq.gene_clustering.asd162.pdf', width = 7, height = 10)
dev.copy2pdf(file='kim_rnaseq_0517.pcl_egl_bulk.kim_gl_061417.pdf', width = 7, height = 10)

###EGL PCL#########
##read in count and metadata
countData1 <- read.table('kim_rnaseq_1016.star_fc.egl_pcl.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_1016.star_fc.egl_pcl.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~sex + rnaaccess_batch + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="kim_rnaseq_1016.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_rnaseq_1016.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$individual, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0517.pcl_egl.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0517.pcl_egl.sample_pca.tissue_only.pdf') 
##differeential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "EGL", "PCL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18515,]
write.csv(resOrdered2DF, file="kim_rnaseq_0517.pcl_egl.de.csv")



