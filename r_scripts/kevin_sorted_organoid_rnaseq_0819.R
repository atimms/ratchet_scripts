library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/kevin_rnaseq_0819";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('sorted_cell_0819_with_organoid.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sorted_cell_0819_with_organoid.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##using all 'timepoints'
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ experiment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="sorted_cell_0819_with_organoid.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sorted_cell_0819_with_organoid.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$sample_name, rld$experiment ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='sorted_cell_0819_with_organoid.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("experiment"))
ggsave('sorted_cell_0819_with_organoid.experiment_pca.pdf')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('sorted_cell_0819_with_organoid.sample_pca.pdf')


##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)["experiment"])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='sorted_cell_0819_with_organoid.25_var_gene_clustering.pdf', width = 7, height = 5)

####from gene list
genes <- read.csv('tim_marker_genes_1019.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize

foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["experiment"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='sorted_cell_0819_with_organoid.marker_genes.clustered.pdf', width = 7, height = 10)
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5, cluster_rows = F)
dev.copy2pdf(file='sorted_cell_0819_with_organoid.marker_genes.not_clustered.pdf', width = 7, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "C4_GFP", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:26437,]
write.csv(resOrdered2DF, file="sorted_cell_0819_with_organoid.C4_GFP_vs_RC28.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "C4_rhoK", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:26437,]
write.csv(resOrdered2DF, file="sorted_cell_0819_with_organoid.C4_rhoK_vs_RC28.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "C5_GFP", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:26437,]
write.csv(resOrdered2DF, file="sorted_cell_0819_with_organoid.C5_GFP_vs_RC28.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "C5_rhoK", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:26437,]
write.csv(resOrdered2DF, file="sorted_cell_0819_with_organoid.C5_rhoK_vs_RC28.csv")


##using just week 28 organoid
##read in count and metadata
countData1 <- read.table('sorted_cell_0819_with_organoid28.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sorted_cell_0819_with_organoid28.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ experiment)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="sorted_cell_0819_with_organoid28.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sorted_cell_0819_with_organoid28.deseq.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("experiment"))
ggsave('sorted_cell_0819_with_organoid28.experiment_pca.pdf')
plotPCA(rld, intgroup = c("sample_name"))
ggsave('sorted_cell_0819_with_organoid28.sample_pca.pdf')

####from gene list
genes <- read.csv('tim_marker_genes_1019.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["experiment"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=8, fontsize_col=10)
dev.copy2pdf(file='sorted_cell_0819_with_organoid28.marker_genes.clustered.pdf', width = 7, height = 10)
##repaeat using genes in order given, CNGT2 not in dataset 
gene_order = c("MAP2","TUBB3","RHO","PDE6A","AIPL1","SAG","CNGB1","NRL","RP1","GNB1","CNGA1","PRPH2","OPN1SW","OPN1MW","OPN1LW","ARR3","GUCA1C","KCNJ10","ASCL1","SOX2","RGR","GLUL","APOE","PMP2","CLRN1","USH1G","USH1C","CDH23","CFI","RLBP1","GFAP","SLC1A3","PRKCA","GRIK1","VSX2","ONECUT1","LHX1","ONECUT2","LNP1","GAD1","TFAP2A","TFAP2B","SLC6A9","POU4F2","GAP43","NEFL","NEFH","NEFM")
newmat <- assay(rld)[gene_order,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["experiment"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=8, fontsize_col=10, cluster_rows = F)
dev.copy2pdf(file='sorted_cell_0819_with_organoid28.marker_genes.not_clustered.pdf', width = 7, height = 10)


##differential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "CRALBP-GFP", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:24875,]
write.csv(resOrdered2DF, file="sorted_cell_0819_with_organoid28.CRALBP-GFP_vs_RC28.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("experiment", "rhoK-tdT", "RC28"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:24875,]
write.csv(resOrdered2DF, file="sorted_cell_0819_with_organoid28.rhoK-tdT_vs_RC28.csv")
