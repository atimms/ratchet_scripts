library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/lan_mouse_hypoxia_rnaseq_0719";
setwd(workingDir);

###all samples
##read in count and metadata
countData1 <- read.table('lan_mouse_p15_hypoxia_0719_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('lan_mouse_p15_hypoxia_0719_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sex + litter + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="lan_mouse_p15_hypoxia_0719_all.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="lan_mouse_p15_hypoxia_0719_all.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sex,  rld$litter, rld$treatment ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='lan_mouse_p15_hypoxia_0719_all.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("litter"))
ggsave('lan_mouse_p15_hypoxia_0719_all.litter_pca.pdf')
plotPCA(rld, intgroup = c("treatment"))
ggsave('lan_mouse_p15_hypoxia_0719_all.treatment_pca.pdf')
plotPCA(rld, intgroup = c("sex"))
ggsave('lan_mouse_p15_hypoxia_0719_all.sex_pca.pdf')

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sex","litter", "treatment")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='lan_mouse_p15_hypoxia_0719_all.25_var_gene_clustering.pdf', width = 7, height = 5)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("treatment",  "Hypoxia", "Control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18321,]
write.csv(resOrdered2DF, file="lan_mouse_p15_hypoxia_0719_all.Hypoxia_vs_Control.csv")

##read in count and metadata for ind litters
##litter1
countData1 <- read.table('lan_mouse_p15_hypoxia_0719_l1.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('lan_mouse_p15_hypoxia_0719_l1.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sex + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("treatment",  "Hypoxia", "Control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18207,]
write.csv(resOrdered2DF, file="lan_mouse_p15_hypoxia_0719_l1.Hypoxia_vs_Control.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##principal components analysis 
plotPCA(rld, intgroup = c("treatment"))
ggsave('lan_mouse_p15_hypoxia_0719_l1.treatment_pca.pdf')
plotPCA(rld, intgroup = c("sex"))
ggsave('lan_mouse_p15_hypoxia_0719_l1.sex_pca.pdf')


##read in count and metadata for ind litters
##litter1
countData1 <- read.table('lan_mouse_p15_hypoxia_0719_l2.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('lan_mouse_p15_hypoxia_0719_l2.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sex + treatment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("treatment",  "Hypoxia", "Control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18207,]
write.csv(resOrdered2DF, file="lan_mouse_p15_hypoxia_0719_l2.Hypoxia_vs_Control.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##principal components analysis 
plotPCA(rld, intgroup = c("treatment"))
ggsave('lan_mouse_p15_hypoxia_0719_l2.treatment_pca.pdf')
plotPCA(rld, intgroup = c("sex"))
ggsave('lan_mouse_p15_hypoxia_0719_l2.sex_pca.pdf')

