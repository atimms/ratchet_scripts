library(WGCNA);
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/home/atimms/ipex_rnaseq";
setwd(workingDir);

##read in count and metadata
#countData1 <- read.table('ipex_gene_edit_0716.htseq_counts.txt', header=T, row.names=1)
#colData1 <- read.table('ipex_gene_edit_0716.metadata.txt', header=T, row.names=1)
#colData1 <- read.table('ipex_gene_edit_0716_with_numbers.metadata.txt', header=T, row.names=1)
##just ctl data
countData1 <- read.table('ipex_gene_edit_0716.ctls.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('ipex_gene_edit_0716.ctls.metadata.with_numbers.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ group)
dds

##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)

##write normalized counts
dds <- estimateSizeFactors(dds)
count_data <- counts(dds, normalized=TRUE)
write.csv(count_data, file="ipex.all_samples.norm_counts.csv")

##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)

##and write to csv file
#write.csv(assay(rld), file="lm.deseq.rlog_counts.csv")

##compare log transformation with rlog - not really needed
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)



##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- (rld$group )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#dev.copy2pdf(file='ipex.sample_heatmap.pdf', width = 7, height = 5)
dev.copy2pdf(file='ipex_ctl.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
plotPCA(rld, intgroup = "group")
#ggsave('ipex.sample_pca.pdf') 
ggsave('ipex_ctl.sample_pca.pdf') 

##multidimensional scaling plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=group)) + geom_point(size=3)
#ggsave('ipex.sample_mds.pdf') 
ggsave('ipex_ctl.sample_mds.pdf') 


##gene clustering
##take 25 most variable genes
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("group", "sample_number")])
pheatmap(mat, annotation_col=df, fontsize_row=10, fontsize_col=10)
#dev.copy2pdf(file='ipex.deseq.gene_clustering.25.pdf', width = 7, height = 5)
dev.copy2pdf(file='ipex_ctl.deseq.gene_clustering.25.pdf', width = 7, height = 5)
##or 50 most variable genes
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("group", "sample_number")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#dev.copy2pdf(file='ipex.deseq.gene_clustering.50.pdf', width = 7, height = 5)
dev.copy2pdf(file='ipex_ctl.deseq.gene_clustering.50.pdf', width = 7, height = 5)

####from gene list
genes <- read.csv('ipex_genelist_1116.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize

foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("group", "sample_number")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=10, fontsize_col=10)
dev.copy2pdf(file='ipex_1016.deseq.gene_clustering.gl_1116.pdf', width = 7, height = 10)





##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
summary(res)

##to get a specific test i.e. pk vs gr
res2 <- results(dds, contrast=c("group", "ctl_edit", "ctl_mock_edit"))
##get summary
summary(res2) #lots of significant genes

##subset genes to get most significant and sort by fold change
resSig <- subset(res2, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

##save differentially expression results
##sort results by adjusted p-value
res2Ordered <- res2[order(res2$padj),]
head(res2Ordered)
##save results as dataframe and take top 2k results, then write csv file
resOrderedDF <- as.data.frame(res2Ordered)[1:2000,]
write.csv(resOrderedDF, file="ipex.de_results.ctl_edit_vs_ctl_mock_edit.csv")


##to get a specific test i.e. pk vs gr
res3 <- results(dds, contrast=c("group", "ctl_edit", "ctl_activated_Treg"))
##get summary
summary(res3) #lots of significant genes

##subset genes to get most significant and sort by fold change
resSig <- subset(res3, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

##save differentially expression results
##sort results by adjusted p-value
res3Ordered <- res3[order(res3$padj),]
head(res3Ordered)
##save results as dataframe and take top 1k results, then write csv file
resOrderedDF <- as.data.frame(res3Ordered)[1:2500,]
write.csv(resOrderedDF, file="ipex.de_results.ctl_edit_vs_ctl_activated_Treg.csv")

resOrderedDF <- as.data.frame(res3Ordered)[1:20000,]
write.csv(resOrderedDF, file="ipex.de_results.ctl_edit_vs_ctl_activated_Treg.all_genes.csv")

##to get a specific test i.e. pk vs gr
res4 <- results(dds, contrast=c("group", "ctl_mock_edit", "ctl_activated_Treg"))
##get summary
summary(res4) #lots of significant genes

##subset genes to get most significant and sort by fold change
resSig <- subset(res4, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

##save differentially expression results
##sort results by adjusted p-value
res4Ordered <- res4[order(res4$padj),]
head(res4Ordered)
##save results as dataframe and take top 1k results, then write csv file
resOrderedDF <- as.data.frame(res4Ordered)[1:2500,]
write.csv(resOrderedDF, file="ipex.de_results.ctl_mock_edit_vs_ctl_activated_Treg.csv")

resOrderedDF <- as.data.frame(res4Ordered)[1:20000,]
write.csv(resOrderedDF, file="ipex.de_results.ctl_mock_edit_vs_ctl_activated_Treg.all_genes.csv")







##plot all points for individual genes
##better with ggplot
gene_list <- c("IL2RA", "IL7", "IL7R", "CTLA4", "IL2", "IL4", "CD40LG", "CCR7", "CCR4", "LAG3", "IFNG", "IL10", "FOXP3", "TGFB1", "IKZF2")
for (i in gene_list) {
  data <- plotCounts(dds, gene=i, intgroup="group", returnData=TRUE)
  ggplot(data, aes(x=group, y=count,color=group)) +
    scale_y_log10() + 
    geom_point(position=position_jitter(width=.1,height=0), size=3)
  pdf_file = paste(i, 'all_data.pdf',sep=".")
  ggsave(pdf_file) 
}



##different samples 11/16

##read in count and metadata
countData1 <- read.table('ipex_gene_edit_0716.tt1.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('ipex_gene_edit_0716.tt1.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ group)
dds

##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)

##gene clustering
####from gene list
genes <- read.csv('ipex_genelist_1116.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize

foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("group", "sample_number")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=10, fontsize_col=10)
dev.copy2pdf(file='ipex_1016.deseq.gene_clustering.with_ipex_patients.gl_1116.pdf', width = 7, height = 10)

