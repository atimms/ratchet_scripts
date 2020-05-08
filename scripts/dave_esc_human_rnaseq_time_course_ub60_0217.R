library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/dave_esc_0117/human_esc_bulk_rnaseq_0217";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('human_esc_rnaseq_time_course_0217.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('human_esc_rnaseq_time_course_0217.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ time_point)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="human_esc_rnaseq_time_course_0217.norm_counts.csv")

##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="esc_rnaseq_0117.deseq.rlog_counts.csv")

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
rownames(sampleDistMatrix) <- paste( rld$time_point,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='human_esc_rnaseq_time_point_0217.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
plotPCA(rld, intgroup = c("time_point"))
ggsave('human_esc_rnaseq_time_point_0217.sample_pca.pdf') 

##multidimensional scaling plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=time_point)) + geom_point(size=3)
ggsave('human_esc_rnaseq_time_point_0217.sample_mds.pdf') 

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)["time_point"])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='human_esc_rnaseq_time_point_0217.25_var.gene_clustering.pdf', width = 7, height = 5)


##take genes from het list
####from gene list
genes <- read.csv('top50_human_gene_names.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize

foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["time_point"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=7, fontsize_col=7)
dev.copy2pdf(file='human_esc_rnaseq_time_point_0217.top_50_het_genes.pdf', width = 7, height = 10)




##differeential expression

##just use early/late classification
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds2 <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ timing)
dds2

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds2)
dds2 <- dds2[ rowSums(counts(dds2)) > 1, ]
nrow(dds2)


##do the test
dds2 <- DESeq(dds2)
##get results and summary
##this just gets the last test
(res2 <- results(dds2))


##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:17722,]
write.csv(resOrdered2DF, file="human_esc_rnaseq_time_point_0217.early_vs_late.csv")
