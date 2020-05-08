library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/milena_rnaseq_0818";
setwd(workingDir);

###all samples
##read in count and metadata
countData1 <- read.table('milena_rnaseq_0818_all_samples.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('milena_rnaseq_0818_all_samples.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~timepoint+ genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)

##principal components analysis 
pcaData <- plotPCA(rld, intgroup=c("genotype", "timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=timepoint)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  scale_colour_manual(values=c("blue", "red", "green"))
ggsave('milena_rnaseq_0818.all_samples.genotype_timepoint.pca.pdf', width = 7, height = 5) 



##gene clustering
##take 25 most variable gene for all samples
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
#topVarGenes_group <- head(order(rowVars(group_means),decreasing=TRUE),25)
##used python to get means of different groups, and manually made new metadata
##get into R
group_means <- read.csv('milena_rnaseq_0818_all_samples.deseq.vst_counts_by_group.csv', header=T, row.names=1)
group_mean_meta <- read.table('milena_rnaseq_0818.group_metadata.txt', header=T, row.names=1) 
##get new data
mat <- group_means[ topVarGenes, ]
#mat <- group_means[ topVarGenes_group, ]
##make the heatmap
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(group_mean_meta[c("timepoint", "genotype")])
col.pal <- brewer.pal(9,"Blues")
annotation_colors = list(
  genotype = c(WT="green", LIG4_Mutant="red", LIG4_Corrected='blue'),
  timepoint = c(Day_4="pink", Day_6="orange", Day_22='yellow'))
pheatmap(mat, annotation_col=newdf, annotation_colors = annotation_colors, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='milena_rnaseq_0818.groups_combined.25_var_genes.heatmap.pdf', width = 14, height = 10)
##take genes of intrest
gene_order <- c('APEX1', 'ATM', 'ATR', 'BRCA1', 'BRCA2', 'CETN2', 'ERCC1', 'XRCC6', 'XRCC5', 'LIG3', 'LIG4', 'MGMT', 'MPG', 'MRE11A', 'MSH2', 'MSH6', 'NBN', 'PCNA', 'RAD50', 'RAD51', 'RAD52', 'TP53BP1', 'XPC', 'XRCC1', 'XRCC4', 'AIFM1', 'BAD', 'BAK1', 'BAX', 'BID', 'CAD', 'CASP3', 'CASP9', 'DIABLO', 'BBC3', 'TGFB1', 'TNF', 'BCL2', 'BCL2L1', 'MCL1', 'BIRC5', 'ASPM', 'CCND1', 'CDK1', 'MKI67', 'NOTCH1', 'BTG1', 'BTG2', 'BTG3', 'BTG4', 'GADD45A', 'ARHGEF2', 'DOCK7', 'ETV5', 'GOLGA2', 'INSC', 'PARD3', 'RAB10', 'SOX5', 'TEAD3', 'ZBTB16', 'NES', 'DCX', 'NEUROD1', 'PAX6', 'SOX1', 'SOX2', 'TUBB3', 'POU3F2', 'BCL11B', 'MAP2', 'RBFOX3', 'SATB2', 'TBR1', 'EOMES')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
##used python to get means of different groups, and manually made new metadata
##get into R
group_means <- read.csv('milena_rnaseq_0818_all_samples.deseq.vst_counts_by_group.csv', header=T, row.names=1)
group_mean_meta <- read.table('milena_rnaseq_0818.group_metadata.txt', header=T, row.names=1) 
##get new data
mat <- group_means[ gene_order, ]
newmat <- assay(rld)[gene_order,] #subset to target genes
#mat <- group_means[ topVarGenes_group, ]
##make the heatmap
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(group_mean_meta[c("timepoint", "genotype")])
annotation_colors = list(
  genotype = c(WT="green", LIG4_Mutant="red", LIG4_Corrected='blue'),
  timepoint = c(Day_4="pink", Day_6="orange", Day_22='yellow'))
pheatmap(mat, annotation_col=newdf, annotation_colors = annotation_colors, fontsize_row=8, fontsize_col=8, cluster_rows=F)
dev.copy2pdf(file='milena_rnaseq_0818.groups_combined.genelist.heatmap.pdf', width = 14, height = 10)


##gene clustering on just days 4 and 6
##take genes of intrest
gene_order <- c('APEX1', 'ATM', 'ATR', 'BRCA1', 'BRCA2', 'CETN2', 'ERCC1', 'XRCC6', 'XRCC5', 'LIG3', 'LIG4', 'MGMT', 'MPG', 'MRE11A', 'MSH2', 'MSH6', 'NBN', 'PCNA', 'RAD50', 'RAD51', 'RAD52', 'TP53BP1', 'XPC', 'XRCC1', 'XRCC4', 'AIFM1', 'BAD', 'BAK1', 'BAX', 'BID', 'CAD', 'CASP3', 'CASP9', 'DIABLO', 'BBC3', 'TGFB1', 'TNF', 'BCL2', 'BCL2L1', 'MCL1', 'BIRC5', 'ASPM', 'CCND1', 'CDK1', 'MKI67', 'NOTCH1', 'BTG1', 'BTG2', 'BTG3', 'BTG4', 'GADD45A', 'ARHGEF2', 'DOCK7', 'ETV5', 'GOLGA2', 'INSC', 'PARD3', 'RAB10', 'SOX5', 'TEAD3', 'ZBTB16', 'NES', 'DCX', 'NEUROD1', 'PAX6', 'SOX1', 'SOX2', 'TUBB3', 'POU3F2', 'BCL11B', 'MAP2', 'RBFOX3', 'SATB2', 'TBR1', 'EOMES')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
##used python to get means of different groups, and manually made new metadata
##get into R
group_means <- read.csv('milena_rnaseq_0818_day4_6.deseq.vst_counts_by_group.csv', header=T, row.names=1)
group_mean_meta <- read.table('milena_rnaseq_0818.day4_6.group_metadata.txt', header=T, row.names=1) 
##get new data
mat <- group_means[ gene_order, ]
newmat <- assay(rld)[gene_order,] #subset to target genes
#mat <- group_means[ topVarGenes_group, ]
##make the heatmap
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(group_mean_meta[c("timepoint", "genotype")])
annotation_colors = list(
  genotype = c(WT="green", LIG4_Mutant="red", LIG4_Corrected='blue'),
  timepoint = c(Day_4="pink", Day_6="orange"))
pheatmap(mat, annotation_col=newdf, annotation_colors = annotation_colors, fontsize_row=8, fontsize_col=8, cluster_rows=F)
dev.copy2pdf(file='milena_rnaseq_0818.day4_6.collapsed.genelist.heatmap.pdf', width = 14, height = 10)

##and get for top 25 variable genes
group_means <- read.csv('milena_rnaseq_0818_day4_6.deseq.vst_counts_by_group.csv', header=T, row.names=1)
group_mean_meta <- read.table('milena_rnaseq_0818.day4_6.group_metadata.txt', header=T, row.names=1) 
##take 25 most variable gene for all samples
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- group_means[ topVarGenes, ]
#mat <- group_means[ topVarGenes_group, ]
##make the heatmap
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(group_mean_meta[c("timepoint", "genotype")])
col.pal <- brewer.pal(9,"Blues")
annotation_colors = list(
  genotype = c(WT="green", LIG4_Mutant="red", LIG4_Corrected='blue'),
  timepoint = c(Day_4="pink", Day_6="orange"))
pheatmap(mat, annotation_col=newdf, annotation_colors = annotation_colors, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='milena_rnaseq_0818.day4_6.collapsed.25_var_genes.heatmap.pdf', width = 14, height = 10)


