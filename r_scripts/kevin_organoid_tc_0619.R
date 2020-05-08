library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kevin_organoid_tc_0619";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('kevin_organiod_tc_0619.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kevin_organiod_tc_0619.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##using all 'timepoints'
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~experiment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="kevin_organiod_tc_0619.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kevin_organiod_tc_0619.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$experiment, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kevin_organiod_tc_0619.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("experiment"))
ggsave('kevin_organiod_tc_0619.experiment_pca.pdf')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('kevin_organiod_tc_0619.sample_pca.pdf')
plotPCA(rld, intgroup = c("timepoint"))
ggsave('kevin_organiod_tc_0619.timepoint_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("experiment", "timepoint")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='kevin_organiod_tc_0619.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "IPSC", "RC5"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.IPSC_vs_RC5.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "IPSC", "RC20"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.IPSC_vs_RC20.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "IPSC", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.IPSC_vs_RC28.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "RC5", "RC20"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.RC5_vs_RC20.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "RC5", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.RC5_vs_RC28.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "RC20", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.RC20_vs_RC28.csv")



##change to using early and late 
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~timepoint)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("timepoint", "IPSC", "early"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.IPSC_vs_early.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("timepoint", "IPSC", "late"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.IPSC_vs_late.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("timepoint", "early", "late"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:25783,]
write.csv(resOrdered2DF, file="kevin_organiod_tc_0619.early_vs_late.csv")
