library(WGCNA);
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/kim_rnaseq_1016";
setwd(workingDir);

##read in count and metadata
#countData1 <- read.table('kim_rnaseq_1016.star_fc.counts2.txt', header=T, row.names=1)
#colData1 <- read.table('kim_rnaseq_1016.star_fc.metadata2.txt', header=T, row.names=1)
#countData1 <- read.table('kim_rnaseq_1016.star_fc.counts3.txt', header=T, row.names=1)
#colData1 <- read.table('kim_rnaseq_1016.star_fc.metadata3.txt', header=T, row.names=1)
countData1 <- read.table('kim_rnaseq_1016.star_fc.counts5.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_1016.star_fc.metadata5.txt', header=T, row.names=1)

head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~sex + rnaaccess_batch + tissue)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
colData(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="kim_rnaseq_1016.norm_counts.csv")

##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_1016.deseq.rlog_counts.csv")

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
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$sex, rld$rnaaccess_batch, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_1016.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
#plotPCA(rld, intgroup = c("tissue", "sex", "rnaaccess_batch"))
plotPCA(rld, intgroup = "tissue")
#ggsave('kim_rnaseq_1016.pk_egl_whole.sample_pca.tissue_only.pdf') 
ggsave('kim_rnaseq_1016.pk_egl_whole_rl.sample_pca.tissue_only.pdf') 

##multidimensional scaling plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=tissue,shape=rnaaccess_batch)) + geom_point(size=3)
ggsave('kim_rnaseq_1016.sample_mds.pdf') 

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "individual")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='kim_rnaseq_1016.deseq.gene_clustering.pdf', width = 7, height = 5)


##take genes from kim's list
####from gene list
#genes <- read.csv('ASD_GeneList176.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- read.csv('genelist_041317.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize

foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "individual")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='kim_rnaseq_1016.deseq.gene_clustering.genes0413.pdf', width = 7, height = 10)


##differeential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))

##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "EGL", "PK"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19390,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.egl_vs_pk.csv")

##to get a specific test:
res3 <- results(dds, contrast=c("tissue", "Whole", "PK"))
##get summary
summary(res3) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19390,]
write.csv(resOrdered3DF, file="kim_rnaseq_1016.de.whole_vs_pk.csv")

##to get a specific test:
res4 <- results(dds, contrast=c("tissue", "Whole", "EGL"))
##get summary
summary(res4) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered4 <- res4[order(res4$padj),]
head(resOrdered4)
##save results as dataframe and take top 20k results, then write csv file
resOrdered4DF <- as.data.frame(resOrdered4)[1:19390,]
write.csv(resOrdered4DF, file="kim_rnaseq_1016.de.whole_vs_egl.csv")

##to get a specific test:
res5 <- results(dds, contrast=c("tissue", "RL", "PK"))
##get summary
summary(res5) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered5 <- res5[order(res5$padj),]
head(resOrdered5)
##save results as dataframe and take top 20k results, then write csv file
resOrdered5DF <- as.data.frame(resOrdered5)[1:19483,]
write.csv(resOrdered5DF, file="kim_rnaseq_1016.de.rl_vs_pk.csv")

##to get a specific test:
res6<- results(dds, contrast=c("tissue", "RL", "Whole"))
##get summary
summary(res6) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered6 <- res6[order(res6$padj),]
head(resOrdered6)
##save results as dataframe and take top 20k results, then write csv file
resOrdered6DF <- as.data.frame(resOrdered6)[1:19483,]
write.csv(resOrdered6DF, file="kim_rnaseq_1016.de.rl_vs_whole.csv")

##to get a specific test:
res7 <- results(dds, contrast=c("tissue", "RL", "EGL"))
##get summary
summary(res7) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered7 <- res7[order(res7$padj),]
head(resOrdered7)
##save results as dataframe and take top 20k results, then write csv file
resOrdered7DF <- as.data.frame(resOrdered7)[1:19483,]
write.csv(resOrdered7DF, file="kim_rnaseq_1016.de.rl_vs_egl.csv")


###just pk and gr
##read in count and metadata
countData1 <- read.table('kim_rnaseq_1016.star_fc.counts4.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_1016.star_fc.metadata4.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~sex + rnaaccess_batch + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
colData(dds)
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##take genes from kim's list
####from gene list
#genes <- read.csv('gsea_118_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
#genes <- read.csv('core_40_pk_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- read.csv('gensat_genelist_0517.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- read.csv('dougherty_genelist_0517.txt', header = FALSE, col.names = 'genes', as.is = TRUE)

genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "individual")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
#dev.copy2pdf(file='kim_rnaseq_1016.deseq.gene_clustering.gsea118_genes.pdf', width = 7, height = 10)
#dev.copy2pdf(file='kim_rnaseq_1016.deseq.gene_clustering.core_40_pk_genes.pdf', width = 7, height = 10)
dev.copy2pdf(file='kim_rnaseq_1016.deseq.gene_clustering.gensat_genes.pdf', width = 7, height = 10)
dev.copy2pdf(file='kim_rnaseq_1016.deseq.gene_clustering.dougherty_genes.pdf', width = 7, height = 10)
