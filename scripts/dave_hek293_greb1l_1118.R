library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/dave_hek293_greb1l_1118";
setwd(workingDir);

###all samples
##read in count and metadata
countData1 <- read.table('hek293_greb1l_1118.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('hek293_greb1l_1118.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~transfected+ genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)

##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="hek293_greb1l_rnaseq_1118.all_samples.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$transfected, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='hek293_greb1l_rnaseq_1118.all_samples.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("transfected"))
ggsave('hek293_greb1l_rnaseq_1118.all_samples.transfected_pca.pdf') 
plotPCA(rld, intgroup = c("genotype"))
ggsave('hek293_greb1l_rnaseq_1118.all_samples.genotype_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("transfected", "genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='hek293_greb1l_rnaseq_1118.all_samples.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "WT", "KO"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19407,]
write.csv(resOrdered2DF, file="hek293_greb1l_rnaseq_1118.WT_vs_KO.csv")


##receode the genotype so can differetiate WT/KO
###all samples
##read in count and metadata
countData1 <- read.table('hek293_greb1l_1118.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('hek293_greb1l_1118.star_fc.metadata_trans.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "WT_T", "KO_T"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19407,]
write.csv(resOrdered2DF, file="hek293_greb1l_rnaseq_1118.not_transfected.WT_vs_KO.csv")


##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "KO", "KO_T"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19407,]
write.csv(resOrdered2DF, file="hek293_greb1l_rnaseq_1118.KO_vs_KOtransfected.csv")
