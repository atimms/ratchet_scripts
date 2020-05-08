library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/greb1l_rnaseq_0817";
setwd(workingDir);


##read in count and metadata
countData1 <- read.table('greb1l_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('greb1l_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="greb1l_rnaseq_0817.norm_counts.csv")

##rlog transform
rld<- rlogTransformation(dds, blind=TRUE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="greb1l_rnaseq_0817.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='greb1l_rnaseq_0817.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "genotype")
ggsave('greb1l_rnaseq_0817.sample_pca.pdf') 
#gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)["genotype"])
pheatmap(mat, annotation_col=df, fontsize_row=10, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_0817.25vargene_clustering.pdf', width = 7, height = 5)
####from gene list
genes <- read.csv('smrt_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_0817.gene_clustering.smrt.pdf', width = 7, height = 5)
genes <- read.csv('ra_reg_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_0817.gene_clustering.ra_reg.pdf', width = 7, height = 5)
genes <- read.csv('up_in_3_geo.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=0.1, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_0817.gene_clustering.in_3_geo.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "WT", "KO"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:16848,]
write.csv(resOrdered2DF, file="greb1l_rnaseq_0817.wt_ko.de.csv")
plotCounts(dds, gene = "Greb1l", intgroup=c("genotype"))
dev.copy2pdf(file='greb1l_rnaseq_0817.Greb1l_expression.pdf', width = 7, height = 5)


##comapre against GSE75616
##read in count and metadata
countData1 <- read.table('greb1l_gse75616.counts.txt', header=T, row.names=1)
colData1 <- read.table('greb1l_gse75616.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##rlog transform
rld <- rlog(dds, blind=FALSE)
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='greb1l_rnaseq_0817_gse75616.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "genotype")
ggsave('greb1l_rnaseq_0817_gse75616.sample_pca.pdf') 
#gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)["genotype"])
pheatmap(mat, annotation_col=df, fontsize_row=10, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_0817_gse75616.25vargene_clustering.pdf', width = 7, height = 5)

