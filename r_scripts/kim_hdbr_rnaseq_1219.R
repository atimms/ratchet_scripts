library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kim_hdbr_rnaseq_1219";
setwd(workingDir);


##original data
##read in count and metadata
countData1 <- read.table('kim_hdbr_rnaseq_all_1219.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_hdbr_rnaseq_all_1219.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ developmental_stage + organism_part)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 20, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="kim_hdbr_rnaseq_1219.norm_counts.csv")
