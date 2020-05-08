ibrary("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/dave_greb1l_trna_1118";
setwd(workingDir);


##all timepoints
##read in count and metadata
countData1 <- read.table('greb1l_all_1118.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('greb1l_all_1118.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert day to a factor
colData1$day <- as.factor(colData1$day)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~genotype + day)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="greb1l_rnaseq_1017.norm_counts.csv")
##rlog transform
rld<- rlogTransformation(dds, blind=TRUE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="greb1l_rnaseq_1017.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- c(rld$genotype, rld$day)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='greb1l_trna_1118.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype", "day"))
ggsave('greb1l_rnaseq_1017.sample_pca.pdf') 
#gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(mat, annotation_col=df, fontsize_row=10, fontsize_col=5)
dev.copy2pdf(file='greb1l_trna_1118.25vargene_clustering.pdf', width = 7, height = 5)
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
resOrdered2DF <- as.data.frame(resOrdered2)[1:201,]
write.csv(resOrdered2DF, file="greb1l_trna_1118.wt_ko.de.csv")
