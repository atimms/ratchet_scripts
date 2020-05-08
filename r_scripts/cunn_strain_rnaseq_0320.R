library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/cunn_strain_rnaseq_0320";
setwd(workingDir);

##get differential expression per sample comparing strains
##read in count and metadata
countData1 <- read.table('strain_rnaseq_1017.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_1017.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18583,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1017.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18583,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.1017.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18583,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.1017.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_1032.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_1032.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19958,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1032.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19958,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.1032.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19958,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.1032.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_3007.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_3007.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18388,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3007.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18388,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.3007.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18388,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.3007.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_3035.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_3035.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18724,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3035.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18724,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.3035.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18724,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.3035.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_4025.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_4025.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19797,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4025.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19797,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.4025.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19797,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.4025.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_4032.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_4032.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19342,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4032.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19342,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.4032.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19342,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.4032.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_AUT25.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_AUT25.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:20189,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.AUT25.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20189,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.AUT25.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:20189,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.AUT25.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1625.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1625.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19406,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1625.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19406,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C1625.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19406,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C1625.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1653.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1653.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19266,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1653.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19266,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C1653.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19266,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C1653.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1671.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1671.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19355,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1671.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19355,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C1671.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19355,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C1671.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1856.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1856.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18946,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1856.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18946,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C1856.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18946,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C1856.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1859.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1859.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18871,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1859.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18871,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C1859.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18871,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C1859.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1926.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1926.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19096,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1926.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19096,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C1926.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19096,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C1926.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1957.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1957.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18479,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1957.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18479,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C1957.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18479,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C1957.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C2038.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C2038.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18800,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2038.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18800,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C2038.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18800,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C2038.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C2084.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C2084.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18622,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2084.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18622,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C2084.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18622,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C2084.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C3049.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C3049.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19766,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3049.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19766,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C3049.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19766,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C3049.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C3066.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C3066.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19012,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3066.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19012,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.C3066.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19012,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.C3066.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_OST85.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_OST85.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:19675,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.OST85.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19675,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.OST85.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19675,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.OST85.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_STL19.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_STL19.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:18893,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.STL19.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18893,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.STL19.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:18893,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.STL19.strain15_vs_strain100.csv")





###each case vs controls
##read in count and metadata
countData1 <- read.table('strain_rnaseq_1017_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_1017_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21167,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1017_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21167,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1017_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21167,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1017_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1653_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1653_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21328,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1653_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21328,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1653_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21328,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1653_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1856_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1856_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21281,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1856_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21281,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1856_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21281,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1856_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_3007_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_3007_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21204,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3007_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21204,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3007_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21204,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3007_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C2084_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C2084_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21265,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2084_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21265,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2084_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21265,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2084_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1625_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1625_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21565,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1625_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21565,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1625_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21565,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1625_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C3049_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C3049_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21522,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3049_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21522,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3049_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21522,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3049_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C3066_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C3066_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21339,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3066_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21339,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3066_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21339,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C3066_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_1032_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_1032_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21565,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1032_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21565,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1032_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21565,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.1032_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_4032_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_4032_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21522,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4032_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21522,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4032_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21522,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4032_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_3035_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_3035_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21202,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3035_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21202,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3035_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21202,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.3035_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C2038_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C2038_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21241,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2038_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21241,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2038_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21241,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C2038_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_4025_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_4025_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21567,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4025_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21567,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4025_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21567,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.4025_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1671_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1671_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21346,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1671_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21346,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1671_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21346,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1671_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1859_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1859_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21250,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1859_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21250,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1859_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21250,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1859_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1926_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1926_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21271,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1926_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21271,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1926_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21271,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1926_vs_control.strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_C1957_ctrls.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_C1957_ctrls.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ cc_strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. case vs control in 0 strain
res1 <- results(dds, contrast=c("cc_strain", "case_0", "control_0"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21216,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1957_vs_control.strain0.csv")
##to get a specific test i.e. case vs control in 15 strain
res1 <- results(dds, contrast=c("cc_strain", "case_15", "control_15"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21216,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1957_vs_control.strain15.csv")
##to get a specific test i.e. case vs control in 100 strain
res1 <- results(dds, contrast=c("cc_strain", "case_100", "control_100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21216,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.C1957_vs_control.strain100.csv")

##get differential expression per gene group comparing strains
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_AXL.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_AXL.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21427,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.AXL.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21427,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.AXL.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:21427,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.AXL.strain15_vs_strain100.csv")

##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_CTRL.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_CTRL.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:22201,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.CTRL.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22201,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.CTRL.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:22201,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.CTRL.strain15_vs_strain100.csv")

##get differential expression per gene group comparing strains
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_FLNA.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_FLNA.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21666,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.FLNA.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21666,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.FLNA.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:21666,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.FLNA.strain15_vs_strain100.csv")

##get differential expression per gene group comparing strains
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_FLNB.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_FLNB.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21728,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.FLNB.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21728,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.FLNB.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:21728,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.FLNB.strain15_vs_strain100.csv")

##get differential expression per gene group comparing strains
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_FLNC.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_FLNC.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21337,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.FLNC.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21337,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.FLNC.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:21337,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.FLNC.strain15_vs_strain100.csv")


##get differential expression per gene group comparing strains
##read in count and metadata
countData1 <- read.table('strain_rnaseq_0220_PIEZO1.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_PIEZO1.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ strain)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##differential expression
dds <- DESeq(dds)
##to get a specific test i.e. 0 vs 100 strain in all 
res1 <- results(dds, contrast=c("strain", "s0", "s100"))
##get summary
summary(res1) #lots of significant genes
##sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
##save results as dataframe and take top 20k results, then write csv file
resOrdered1DF <- as.data.frame(resOrdered1)[1:21523,]
write.csv(resOrdered1DF, file="strain_rnaseq_0320.PIEZO1.strain0_vs_strain100.csv")
##to get a specific test i.e. 0 vs 15 strain in all 
res2 <- results(dds, contrast=c("strain", "s0", "s15"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21523,]
write.csv(resOrdered2DF, file="strain_rnaseq_0320.PIEZO1.strain0_vs_strain15.csv")
##to get a specific test i.e. 15 vs 100 strain in all 
res3 <- results(dds, contrast=c("strain", "s15", "s100"))
##get summary
summary(res3) #lots of significant genes
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:21523,]
write.csv(resOrdered3DF, file="strain_rnaseq_0320.PIEZO1.strain15_vs_strain100.csv")
