library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/daryl_human_kidney_rnaseq_0719";
setwd(workingDir);

###all samples
##read in count and metadata
countData1 <- read.table('GSE126805_all.gene_counts.formatted.txt', header=T, row.names=1)
colData1 <- read.table('GSE126805.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ kidney + timepoint)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="GSE126805_0719.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="GSE126805_0719.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$kidney, rld$timepoint ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='GSE126805_0719.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("kidney"))
ggsave('GSE126805_0719.kidney_pca.pdf')
plotPCA(rld, intgroup = c("timepoint"))
ggsave('GSE126805_0719.timepoint_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("kidney", "timepoint")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=4)
dev.copy2pdf(file='GSE126805_0719.25_var_gene_clustering.pdf', width = 7, height = 5)

##take from genelist
genes <- read.csv('day2_1.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day2_1.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day2_2.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day2_2.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day2_3.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day2_3.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day2_4.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day2_4.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day2_5.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day2_5.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day2_6.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day2_6.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day5_1.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day5_1.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day5_1.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day5_1.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day5_2.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day5_2.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day5_3.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day5_3.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day5_4.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day5_4.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day5_5.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day5_5.heatmap.pdf', width = 7, height = 10)
genes <- read.csv('day5_6.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["timepoint"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=5)
dev.copy2pdf(file='GSE126805_0719.day5_6.heatmap.pdf', width = 7, height = 10)



##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("timepoint",  "baseline_pre", "baseline_post"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:37414,]
write.csv(resOrdered2DF, file="GSE126805_0719.pre_vs_post.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("timepoint",  "baseline_pre", "3months"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:37414,]
write.csv(resOrdered2DF, file="GSE126805_0719.pre_vs_3months.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("timepoint",  "baseline_pre", "1year"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:37414,]
write.csv(resOrdered2DF, file="GSE126805_0719.pre_vs_1year.csv")


##to get a specific test:
res2 <- results(dds, contrast=c("timepoint",  "baseline_post", "3months"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:37414,]
write.csv(resOrdered2DF, file="GSE126805_0719.post_vs_3months.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("timepoint",  "baseline_post", "1year"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:37414,]
write.csv(resOrdered2DF, file="GSE126805_0719.post_vs_1year.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("timepoint",  "3months", "1year"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:37414,]
write.csv(resOrdered2DF, file="GSE126805_0719.3months_vs_1year.csv")




