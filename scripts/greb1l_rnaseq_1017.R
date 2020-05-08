library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/dave_greb1l_rnaseq_1017";
setwd(workingDir);


##all timepoints
##read in count and metadata
countData1 <- read.table('greb1l_1017.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('greb1l_1017.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert day to a factor
colData1$day <- as.factor(colData1$day)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~genotype + day)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="greb1l_rnaseq_1017.norm_counts.csv")
##rlog transform
rld<- rlogTransformation(dds, blind=TRUE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="greb1l_rnaseq_1017.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- c(rld$genotype, rld$day)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='greb1l_rnaseq_1017.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype", "day"))
ggsave('greb1l_rnaseq_1017.sample_pca.pdf') 
#gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(mat, annotation_col=df, fontsize_row=10, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.25vargene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "WT", "KO"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:16848,]
write.csv(resOrdered2DF, file="greb1l_rnaseq_1017.wt_ko.de.csv")
plotCounts(dds, gene = "Greb1l", intgroup=c("genotype"))
dev.copy2pdf(file='greb1l_rnaseq_1017.Greb1l_expression.pdf', width = 7, height = 5)
####from gene list
genes <- read.csv('smrt_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.smrt.pdf', width = 7, height = 5)
genes <- read.csv('ra_reg_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_0817.gene_clustering.ra_reg.pdf', width = 7, height = 5)
genes <- read.csv('up_in_3_geo.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=0.1, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.in_3_geo.pdf', width = 7, height = 5)
##using genes from gsea leading edge analysis
genes <- read.csv('gsea_le_day4_neg_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.le_day4_neg.pdf', width = 7, height = 5)
genes <- read.csv('gsea_le_day4_pos_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.le_day4_pos.pdf', width = 7, height = 5)
genes <- read.csv('exatlas_genes_1017.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.exatlas_1017.clustered.pdf', width = 7, height = 5)
sample_order = c("C7_D0_1", "C7_D0_2", "C7_D0_3", "F1_D0_1", "F1_D0_2", "F1_D0_3", "C7_D4_1", "C7_D4_2", "C7_D4_3", "F1_D4_1", "F1_D4_2", "F1_D4_3", "C7_1", "C7_2", "C7_3", "F1_1", "F1_2", "F1_3")
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.exatlas_1017.not_clustered.pdf', width = 7, height = 5)
####from gene list --- germlayer markers
tests = c("ectoderm", "endoderm", "mesendoderm", "mesoderm", "pleuripotent", "all_markers")
for(i in 1:length(tests)){
  csv_file = paste(tests[i],"txt",sep=".")
  genes <- read.csv(csv_file, header = FALSE, col.names = 'genes', as.is = TRUE)
  genes <- genes$genes #vectorize
  foo <- assay(rld) #save assay as matrix
  assaygenes <- row.names(foo) #extract rownames
  idx <- assaygenes %in% genes #find target genes
  newmat <- assay(rld)[idx,] #subset to target genes
  newmat <- newmat - rowMeans(newmat)
  newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
  pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
  pdf_file <- paste("greb1l_rnaseq_1017.gene_clustering",tests[i],"pdf",sep=".")
  dev.copy2pdf(file=pdf_file, width = 7, height = 5)
}
##imprinted genes
genes <- read.csv('dave_imprinted_genes_0118.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.imprinted_genes_0118.clustered.pdf', width = 7, height = 5)
sample_order = c("C7_D0_1", "C7_D0_2", "C7_D0_3", "F1_D0_1", "F1_D0_2", "F1_D0_3", "C7_D4_1", "C7_D4_2", "C7_D4_3", "F1_D4_1", "F1_D4_2", "F1_D4_3", "C7_1", "C7_2", "C7_3", "F1_1", "F1_2", "F1_3")
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.imprinted_genes_0118.not_clustered.pdf', width = 7, height = 5)
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5, cluster_cols=F, border_color=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.imprinted_genes_0118.not_clustered.no_border.pdf', width = 7, height = 5)
##bivalent genes 0118
genes <- read.csv('dave_bivalent_genes_0118.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.bivalent_genes_0118.clustered.pdf', width = 7, height = 5)
sample_order = c("C7_D0_1", "C7_D0_2", "C7_D0_3", "F1_D0_1", "F1_D0_2", "F1_D0_3", "C7_D4_1", "C7_D4_2", "C7_D4_3", "F1_D4_1", "F1_D4_2", "F1_D4_3", "C7_1", "C7_2", "C7_3", "F1_1", "F1_2", "F1_3")
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.bivalent_genes_0118.not_clustered.pdf', width = 7, height = 5)
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5, cluster_cols=F, border_color=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.bivalent_genes_0118.not_clustered.no_border.pdf', width = 7, height = 5)
##new set of genes from qpcr set 0118
sample_order = c("F1_D0_1", "F1_D0_2", "F1_D0_3", "C7_D0_1", "C7_D0_2", "C7_D0_3", "F1_D4_1", "F1_D4_2", "F1_D4_3", "C7_D4_1", "C7_D4_2", "C7_D4_3", "F1_1", "F1_2", "F1_3", "C7_1", "C7_2", "C7_3")
genes <- read.csv('qpcr_ectoderm.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.qpcr_ectoderm_0118.pdf', width = 7, height = 5)
genes <- read.csv('qpcr_endoderm.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.qpcr_endoderm_0118.pdf', width = 7, height = 5)
genes <- read.csv('qpcr_mesendoderm.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.qpcr_mesendoderm_0118.pdf', width = 7, height = 5)
genes <- read.csv('qpcr_pluri.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.qpcr_pluri_0118.pdf', width = 7, height = 5)

genes <- read.csv('qpcr_all_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, cluster_cols=F, fontsize_row=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.qpcr_all_0118.pdf', width = 7, height = 5)



##day4 only
##read in count and metadata
countData1 <- read.table('greb1l_d4_1017.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('greb1l_d4_1017.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="greb1l_rnaseq_d4_1017.norm_counts.csv")
##rlog transform
rld<- rlogTransformation(dds, blind=TRUE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="greb1l_rnaseq_d4_1017.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='greb1l_rnaseq_d4_1017.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "genotype")
ggsave('greb1l_rnaseq_d4_1017.sample_pca.pdf') 
#gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)["genotype"])
pheatmap(mat, annotation_col=df, fontsize_row=10, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d4_1017.25vargene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "WT", "KO"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18006,]
write.csv(resOrdered2DF, file="greb1l_rnaseq_d4_1017.wt_ko.de.csv")
plotCounts(dds, gene = "Greb1l", intgroup=c("genotype"))
dev.copy2pdf(file='greb1l_rnaseq_d4_1017.Greb1l_expression.pdf', width = 7, height = 5)
####from gene list
genes <- read.csv('smrt_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d4_1017.gene_clustering.smrt.pdf', width = 7, height = 5)
genes <- read.csv('ra_reg_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d4_1017.gene_clustering.ra_reg.pdf', width = 7, height = 5)
genes <- read.csv('up_in_3_geo.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=0.1, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d4_1017.gene_clustering.in_3_geo.pdf', width = 7, height = 5)


##day0 only
##read in count and metadata
countData1 <- read.table('greb1l_d0_1017.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('greb1l_d0_1017.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="greb1l_rnaseq_d0_1017.norm_counts.csv")
##rlog transform
rld<- rlogTransformation(dds, blind=TRUE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="greb1l_rnaseq_d0_1017.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='greb1l_rnaseq_d0_1017.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "genotype")
ggsave('greb1l_rnaseq_d0_1017.sample_pca.pdf') 
#gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)["genotype"])
pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='greb1l_rnaseq_d0_1017.25vargene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "WT", "KO"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18035,]
write.csv(resOrdered2DF, file="greb1l_rnaseq_d0_1017.wt_ko.de.csv")
plotCounts(dds, gene = "Greb1l", intgroup=c("genotype"))
dev.copy2pdf(file='greb1l_rnaseq_d0_1017.Greb1l_expression.pdf', width = 7, height = 5)
####from gene list
genes <- read.csv('smrt_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d0_1017.gene_clustering.smrt.pdf', width = 7, height = 5)
genes <- read.csv('ra_reg_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d0_1017.gene_clustering.ra_reg.pdf', width = 7, height = 5)
genes <- read.csv('up_in_3_geo.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=0.1, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d0_1017.gene_clustering.in_3_geo.pdf', width = 7, height = 5)
##bivalent genes 0118
genes <- read.csv('dave_bivalent_genes_0118.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_1017.day0.gene_clustering.bivalent_genes_0118.clustered.pdf', width = 7, height = 5)
sample_order = c("C7_D0_1", "C7_D0_2", "C7_D0_3", "F1_D0_1", "F1_D0_2", "F1_D0_3")
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5, cluster_cols=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.day0.gene_clustering.bivalent_genes_0118.not_clustered.pdf', width = 7, height = 5)
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5, cluster_cols=F, border_color=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.day0.gene_clustering.bivalent_genes_0118.not_clustered.no_border.pdf', width = 7, height = 5)




##day8 only
##read in count and metadata
countData1 <- read.table('greb1l_d8_1017.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('greb1l_d8_1017.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~genotype)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="greb1l_rnaseq_d8_1017.norm_counts.csv")
##rlog transform
rld<- rlogTransformation(dds, blind=TRUE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="greb1l_rnaseq_d8_1017.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='greb1l_rnaseq_d8_1017.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = "genotype")
ggsave('greb1l_rnaseq_d8_1017.sample_pca.pdf') 
#gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)["genotype"])
pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='greb1l_rnaseq_d8_1017.25vargene_clustering.pdf', width = 7, height = 5)
##differeential expression
##do the test
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "WT", "KO"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18035,]
write.csv(resOrdered2DF, file="greb1l_rnaseq_d8_1017.wt_ko.de.csv")
plotCounts(dds, gene = "Greb1l", intgroup=c("genotype"))
dev.copy2pdf(file='greb1l_rnaseq_d8_1017.Greb1l_expression.pdf', width = 7, height = 5)
####from gene list
genes <- read.csv('smrt_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d8_1017.gene_clustering.smrt.pdf', width = 7, height = 5)
genes <- read.csv('ra_reg_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d8_1017.gene_clustering.ra_reg.pdf', width = 7, height = 5)
genes <- read.csv('up_in_3_geo.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["genotype"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=0.1, fontsize_col=5)
dev.copy2pdf(file='greb1l_rnaseq_d8_1017.gene_clustering.in_3_geo.pdf', width = 7, height = 5)

##new set of genes 0318
sample_order = c("F1_D0_1", "F1_D0_2", "F1_D0_3", "C7_D0_1", "C7_D0_2", "C7_D0_3", "F1_D4_1", "F1_D4_2", "F1_D4_3", "C7_D4_1", "C7_D4_2", "C7_D4_3", "F1_1", "F1_2", "F1_3", "C7_1", "C7_2", "C7_3")
gene_order = c("Utf1", "Bub1", "Rif1", "Hspa9", "Dppa3", "Ccna2", "Cdk1", "Chd1", "Kat5", "Cdc42", "Sox15", "Dnmt1", "Chd7", "Esrrb", "Smad1", "Fut4", "Tcf3", "Thap11", "Smad3", "Stat3", "Fgf2", "Fgf4", "Gdf3", "Nppb", "T", "Pthlh", "Gata4", "Pdgfra", "Hopx", "Runx1", "Nkx2-5", "Klf5", "Acta2", "Mixl1", "Abca4", "Bmp4", "Tbx3", "Runx2", "Des", "Plvap", "Hey1", "Cdh2", "Pdgfrb", "Snai2", "Foxf1", "Cdx2", "Cd34", "Gata2", "Hand2", "Igf2", "Cdh5", "Hand1", "Gata6", "Sox7", "Sox17", "Hhex", "Hnf1b", "Foxa2", "Eomes", "Rxrg", "Nodal", "Lefty1", "Gsc", "Cdh1", "Prdm1", "Cdh20", "Elavl3", "Lefty2", "Cxcr4", "Cabp7", "Gata1", "Foxa1", "Hnf4a", "Cplx2", "Pou3f3", "Phox2b", "Otx2", "Foxd3", "Prom1", "Sdc2", "En1", "Myo3b", "Dmbx1", "Sox1", "Ncam1", "Pax6", "Map2", "Pou4f1", "Papln", "Tubb3", "Zbtb16", "Lmx1a", "Col2a1", "Drd4", "Wnt1", "Nr2f2", "Col1a1", "Nos2", "Pax3")
genes <- read.csv('germ_layer_makers_0318.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("genotype", "day")])
pheatmap(newmat, annotation_col=newdf, cluster_cols=F, fontsize_row=3, fontsize_col=8, border_color=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.germ_layer_makers_0318.clusterd.pdf', width = 7, height = 5)
##repeat without clustering genes
newmat <- assay(rld)[gene_order,sample_order] #subset to target genes
newmat <- newmat - rowMeans(newmat)
pheatmap(newmat, annotation_col=newdf, cluster_cols=F, cluster_rows=F, fontsize_row=3, fontsize_col=8, border_color=F)
dev.copy2pdf(file='greb1l_rnaseq_1017.gene_clustering.germ_layer_makers_0318.not_clusterd.pdf', width = 7, height = 5)
