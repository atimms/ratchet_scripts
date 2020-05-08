library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kim_rnaseq_0218";
setwd(workingDir);


##whole egl pcl
##read in count and metadata
countData1 <- read.table('kim_rnaseq_0218_wep.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0218_wep.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData1$rnaaccess_batch <- as.factor(colData1$rnaaccess_batch)
colData1$age_pcw <- as.factor(colData1$age_pcw)
colData1$donor <- as.factor(colData1$donor)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~rnaaccess_batch + sex + age_pcw + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##rlog transform and check
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_0218_wep.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$donor, sep="-" )
#rownames(sampleDistMatrix) <- rownames(colData1)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0218_wep.sample_heatmap.pdf', width = 7, height = 5)
#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_wep.pca.tissue.pdf') 
plotPCA(rld, intgroup = "age_pcw")
ggsave('kim_rnaseq_0218_wep.pca.age.pdf') 
plotPCA(rld, intgroup = "donor")
ggsave('kim_rnaseq_0218_wep.pca.donor.pdf') 
plotPCA(rld, intgroup = "rnaaccess_batch")
ggsave('kim_rnaseq_0218_wep.pca.batch.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
##remove first 2 genes as giving strange signal in one sample
#topVarGenes = topVarGenes[3:25]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "donor")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='kim_rnaseq_0218_wep.deseq.gene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "EGL", "PCL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19617,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.egl_vs_pcl.csv")
##to get a specific test:
res3 <- results(dds, contrast=c("tissue", "EGL", "Bulk"))
##get summary
summary(res3) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:19617,]
write.csv(resOrdered3DF, file="kim_rnaseq_1016.de.egl_vs_bulk.csv")
##to get a specific test:
res4 <- results(dds, contrast=c("tissue", "PCL", "Bulk"))
##get summary
summary(res4) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered4 <- res4[order(res4$padj),]
head(resOrdered4)
##save results as dataframe and take top 20k results, then write csv file
resOrdered4DF <- as.data.frame(resOrdered4)[1:19617,]
write.csv(resOrdered4DF, file="kim_rnaseq_1016.de.pcl_vs_bulk.csv")




##other samples
##read in count and metadata
countData2 <- read.table('kim_rnaseq_0218_other.star_fc.counts.txt', header=T, row.names=1)
colData2 <- read.table('kim_rnaseq_0218_other.star_fc.metadata.txt', header=T, row.names=1)
head(countData2)
head(colData2)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData2$rnaaccess_batch <- as.factor(colData2$rnaaccess_batch)
colData2$age_pcw <- as.factor(colData2$age_pcw)
colData2$donor <- as.factor(colData2$donor)
dds <- DESeqDataSetFromMatrix(countData = countData2, colData = colData2, design = ~sex + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##rlog transform and check
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_0218_other.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$donor, sep="-" )
#rownames(sampleDistMatrix) <- rownames(colData1)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0218_other.sample_heatmap.pdf', width = 7, height = 5)
#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_other.pca.tissue.pdf') 
plotPCA(rld, intgroup = "age_pcw")
ggsave('kim_rnaseq_0218_other.pca.age.pdf') 
plotPCA(rld, intgroup = "donor")
ggsave('kim_rnaseq_0218_other.pca.donor.pdf') 
plotPCA(rld, intgroup = "rnaaccess_batch")
ggsave('kim_rnaseq_0218_other.pca.batch.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "donor")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='kim_rnaseq_0218_other.deseq.gene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "RL_svz", "RL_vz"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19014,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.RL_svz_vs_RL_vz.csv")
##to do a lot more:
res2 <- results(dds, contrast=c("tissue", "VZ_ntzegl", "RL"))
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
resOrdered2DF <- as.data.frame(resOrdered2)[1:19014,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.VZ_ntzegl_vs_RL.csv")
res2 <- results(dds, contrast=c("tissue", "VZ_ntzegl", "VZ_svz"))
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
resOrdered2DF <- as.data.frame(resOrdered2)[1:19014,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.VZ_ntzegl_vs_VZ_svz.csv")
res2 <- results(dds, contrast=c("tissue", "VZ_ntzegl", "VZ_vz"))
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
resOrdered2DF <- as.data.frame(resOrdered2)[1:19014,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.VZ_ntzegl_vs_VZ_vz.csv")
res2 <- results(dds, contrast=c("tissue", "RL", "VZ_svz"))
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
resOrdered2DF <- as.data.frame(resOrdered2)[1:19014,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.RL_vs_VZ_svz.csv")
res2 <- results(dds, contrast=c("tissue", "RL", "VZ_vz"))
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
resOrdered2DF <- as.data.frame(resOrdered2)[1:19014,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.RL_vs_VZ_vz.csv")
res2 <- results(dds, contrast=c("tissue", "VZ_svz", "VZ_vz"))
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
resOrdered2DF <- as.data.frame(resOrdered2)[1:19014,]
write.csv(resOrdered2DF, file="kim_rnaseq_1016.de.VZ_svz_vs_VZ_vz.csv")


##all samples so can compare RL with EGL 
##read in count and metadata
countData1 <- read.table('kim_rnaseq_0218_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0218_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData1$rnaaccess_batch <- as.factor(colData1$rnaaccess_batch)
colData1$age_pcw <- as.factor(colData1$age_pcw)
colData1$donor <- as.factor(colData1$donor)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~sex + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##rlog transform and check
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_0218_all.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$donor, sep="-" )
#rownames(sampleDistMatrix) <- rownames(colData1)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0218_all.sample_heatmap.pdf', width = 7, height = 5)
#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_all.pca.tissue.pdf') 
plotPCA(rld, intgroup = "age_pcw")
ggsave('kim_rnaseq_0218_all.pca.age.pdf') 
plotPCA(rld, intgroup = "donor")
ggsave('kim_rnaseq_0218_all.pca.donor.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "donor")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='kim_rnaseq_0218_all.deseq.gene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "EGL", "RL_vz"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20159,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218.de.egl_vs_rl_vz.csv")
##to get a specific test:
res3 <- results(dds, contrast=c("tissue", "EGL", "RL_svz"))
##get summary
summary(res3) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)
##save results as dataframe and take top 20k results, then write csv file
resOrdered3DF <- as.data.frame(resOrdered3)[1:20159,]
write.csv(resOrdered3DF, file="kim_rnaseq_0218.de.egl_vs_rl_svz.csv")


##all samples plus the homebew library prep version so can compare
##read in count and metadata
countData1 <- read.table('kim_rnaseq_0218_homebrew_vs_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0218_homebrew_vs_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData1$rnaaccess_batch <- as.factor(colData1$rnaaccess_batch)
colData1$age_pcw <- as.factor(colData1$age_pcw)
colData1$donor <- as.factor(colData1$donor)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~donor + tissue + library_prep)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##rlog transform and check
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_0218_homebrew_vs_all.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$donor, sep="-" )
#rownames(sampleDistMatrix) <- rownames(colData1)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0218_homebrew_vs_all.sample_heatmap.pdf', width = 7, height = 5)
#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_homebrew_vs_all.pca.tissue.pdf') 
plotPCA(rld, intgroup = "library_prep")
ggsave('kim_rnaseq_0218_homebrew_vs_all.pca.library_prep.pdf') 
plotPCA(rld, intgroup = "donor")
ggsave('kim_rnaseq_0218_homebrew_vs_all.pca.donor.pdf') 



##compare homebrew with all other samples
##read in count and metadata
countData1 <- read.table('kim_rnaseq_0218_homebrew_vs_all_updated.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0218_homebrew_vs_all_updated.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData1$donor <- as.factor(colData1$donor)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~library_prep + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##rlog transform and check
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_rnaseq_0218_homebrew_vs_all.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$library_prep, sep="-" )
#rownames(sampleDistMatrix) <- rownames(colData1)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0218_homebrew_vs_all.sample_heatmap.pdf', width = 7, height = 5)
#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_homebrew_vs_all.pca.tissue.pdf') 
plotPCA(rld, intgroup = "library_prep")
ggsave('kim_rnaseq_0218_homebrew_vs_all.pca.library_prep.pdf') 
plotPCA(rld, intgroup = c("tissue", "library_prep"))
ggsave('kim_rnaseq_0218_homebrew_vs_all.pca.library_prep_tissue.pdf') 

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","library_prep")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='kim_rnaseq_0218_homebrew_vs_all.deseq.gene_clustering.pdf', width = 7, height = 5)


##rl combined set
##read in count and metadata
countData1 <- read.table('kim_rnaseq_0218_rl_combined.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0218_rl_combined.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData1$rnaaccess_batch <- as.factor(colData1$rnaaccess_batch)
colData1$age_pcw <- as.factor(colData1$age_pcw)
colData1$donor <- as.factor(colData1$donor)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~rnaaccess_batch + sex + age_pcw + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 15, ]
nrow(dds)
##rlog transform and check
#rld <- rlog(dds, blind=FALSE)
rld <- vst(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_0218_rl_combined.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$donor, sep="-" )
#rownames(sampleDistMatrix) <- rownames(colData1)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0218_rl_combined.sample_heatmap.pdf', width = 7, height = 5)
#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_rl_combined.pca.tissue.pdf') 
plotPCA(rld, intgroup = "age_pcw")
ggsave('kim_rnaseq_0218_rl_combined.pca.age.pdf') 
plotPCA(rld, intgroup = "donor")
ggsave('kim_rnaseq_0218_rl_combined.pca.donor.pdf') 
plotPCA(rld, intgroup = "rnaaccess_batch")
ggsave('kim_rnaseq_0218_rl_combined.pca.batch.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
##remove first 2 genes as giving strange signal in one sample
#topVarGenes = topVarGenes[3:25]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "donor")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='kim_rnaseq_0218_rl_combined.deseq.gene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test, need to add filtering to remove rows with low counts
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 10
dds <- dds[filter,]
#dds <- DESeq(dds)
##didn't work so use.. https://support.bioconductor.org/p/65091/
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=2000)

##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "EGL", "PCL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.egl_vs_pcl.csv")
res2 <- results(dds, contrast=c("tissue", "RL", "PCL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.rl_vs_pcl.csv")
res2 <- results(dds, contrast=c("tissue", "RL", "EGL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.rl_vs_egl.csv")

##new analysis 0819
res2 <- results(dds, contrast=c("tissue", "RL", "Bulk"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.rl_vs_bulk.csv")
res2 <- results(dds, contrast=c("tissue", "PCL", "Bulk"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.pcl_vs_bulk.csv")
res2 <- results(dds, contrast=c("tissue", "EGL", "Bulk"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.egl_vs_bulk.csv")


