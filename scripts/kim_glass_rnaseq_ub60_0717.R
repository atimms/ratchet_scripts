library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/glass_rnaseq_0717";
setwd(workingDir);


##read in count and metadata
countData1 <- read.table('glass_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('glass_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~sample.1 + library)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="glass_rnaseq_0717.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="glass_rnaseq_071.rlog_counts.csv")

##histogram from featue counts data
de_data = read.table('glass_feature_counts_long.txt', header=T)
ggplot(data=de_data, aes(x=sample, y=percentage, fill=type)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('glass_rnaseq_0717.featue_counts.histogram.pdf') 

##graph normalized counts
norm_data = read.table('brain1_norm_counts.txt', header=T)
ggplot(norm_data, aes(x=rlog_counts, colour=library_prep)) + geom_density()
ggsave('glass_rnaseq_0717.brain1.density_plot.pdf') 
ggplot(norm_data, aes(x=rlog_counts, colour=library_prep)) + geom_histogram(binwidth=1, alpha=.1,position="dodge")
ggsave('glass_rnaseq_0717.brain1.histogram.pdf') 
