library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/cherry_rnaseq_1019";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('cherry_rnaseq_1019_ret_rpe.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cherry_rnaseq_1019_ret_rpe.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sex + ethnicity + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cherry_rnaseq_1019.ret_rpe.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cherry_rnaseq_1019.ret_rpe.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$name, rld$tissue ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cherry_rnaseq_1019.ret_rpe.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("name"))
ggsave('cherry_rnaseq_1019.ret_rpe.name_pca.pdf')
plotPCA(rld, intgroup = c("tissue"))
ggsave('cherry_rnaseq_1019.ret_rpe.tissue_pca.pdf')


##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[,c("tissue","name","sex")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cherry_rnaseq_1019.ret_rpe.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "retina", "rpe"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:27910,]
write.csv(resOrdered2DF, file="cherry_rnaseq_1019.ret_rpe.retina_vs_rpe.csv")

##plot individual genes
plotCounts(dds, gene="PHGDH", intgroup=c("tissue"))
dev.copy2pdf(file='cherry_rnaseq_1019.ret_rpe.PHGDH.pdf')
plotCounts(dds, gene="RPE65", intgroup=c("tissue"))
dev.copy2pdf(file='cherry_rnaseq_1019.ret_rpe.RPE65.pdf')
plotCounts(dds, gene="RHO", intgroup=c("tissue"))
dev.copy2pdf(file='cherry_rnaseq_1019.ret_rpe.RHO.pdf')
plotCounts(dds, gene="NR2E3", intgroup=c("tissue"))
dev.copy2pdf(file='cherry_rnaseq_1019.ret_rpe.NR2E3.pdf')
plotCounts(dds, gene="ABCA4", intgroup=c("tissue"))
dev.copy2pdf(file='cherry_rnaseq_1019.ret_rpe.ABCA4.pdf')

