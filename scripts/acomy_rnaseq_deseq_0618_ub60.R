library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/data/atimms/acomy_rnaseq_0618";
setwd(workingDir);


##mouse data
##read in count and metadata
countData1 <- read.table('acomy_rnaseq_0618_mouse.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('acomy_rnaseq_0618_mouse.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ condition)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="acomy_rnaseq_0618_mouse.norm_counts.csv")
##rlog transform
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="acomy_rnaseq_0618_mouse.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("condition"))
ggsave('acomy_rnaseq_0618_mouse.sample_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "day2"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19100,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_mouse.sham_vs_day2.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "day5"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19100,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_mouse.sham_vs_day5.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "day2", "day5"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19100,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_mouse.day2_vs_day5.csv")

##make heatmap for 6 sets of genes 
####from gene list
genes <- read.csv('Arvantis_2d_8d.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Arvantis_2d_8d.pdf', width = 7, height = 10)
genes <- read.csv('Arvantis_2d_SO.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Arvantis_2d_SO.pdf', width = 7, height = 10)
genes <- read.csv('Arvantis_8d_2d.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Arvantis_8d_2d.pdf', width = 7, height = 10)
genes <- read.csv('Arvantis_8d_SO.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Arvantis_8d_SO.pdf', width = 7, height = 10)
genes <- read.csv('Arvantis_SO_2d.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Arvantis_SO_2d.pdf', width = 7, height = 10)
genes <- read.csv('Arvantis_SO_8d.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Arvantis_SO_8d.pdf', width = 7, height = 10)
##and 4 more
genes <- read.csv('Grgic_d0_d2.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Grgic_d0_d2.pdf', width = 7, height = 10)
genes <- read.csv('Grgic_d0_d5.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Grgic_d0_d5.pdf', width = 7, height = 10)
genes <- read.csv('Grgic_d2_d0.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Grgic_d2_d0.pdf', width = 7, height = 10)
genes <- read.csv('Grgic_d5_d0.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.Grgic_d5_d0.pdf', width = 7, height = 10)
##heatmap for dave - dna transcription
genes <- read.csv('DNA_templated_transcription_gene_list.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.DNA_templated_transcription.pdf', width = 7, height = 10)
##heatmap for sam - myofibroblasts tfs
genes <- read.csv('sam_myo_tf.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=10, fontsize_col=10)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.sam_myo_tf.pdf', width = 7, height = 10)
##heatmap for dave - bivalent genes
genes <- read.csv('dave_bivalent_genes_1118.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
ph = pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.bilvalent_gene.pdf', width = 7, height = 10)
##get the row order
ph_order = ph$tree_row$order
##get the corresponding gene names and use to make acomy gene order the same
gene_names = rownames(newmat)[ph_order]
write.csv(gene_names, file="mouse_gene_order.bilvalent.csv")
##now print in same order as acomy (got manually)
gene_order <- c('Hoxd9', 'Irx2', 'Irx1', 'Gbx2', 'Hoxa2', 'Prox1', 'Pax2', 'Sall3', 'Irx5', 'Hoxd3', 'Kazald1', 'Hoxa4', 'Irx3', 'Hoxb4', 'Hoxd4', 'Hoxa10', 'Hoxa6', 'En2', 'Dach1', 'Lmx1b', 'Nr2f2', 'Hoxc9', 'Hoxb9', 'Eif4e3', 'Cntnap1', 'Coasy', 'Hoxb2', 'Hoxa9', 'Hoxd10', 'Hoxa5', 'Hoxd11', 'Hoxa7', 'Hoxc4', 'Hoxa3', 'Hoxb5', 'Hoxb6', 'Hoxb8', 'Pitx2', 'Shh', 'Hoxb7', 'Foxa1', 'Lmo1', 'Htr7', 'Slc22a4', 'Hoxc5', 'Hoxc10', 'Hoxa1', 'Maf', 'Sox6', 'Ebf1', 'Pax9', 'Hoxa11', 'Foxd3', 'Dlx1', 'Nkx6-1', 'Zic1', 'Evx1', 'Hoxd13', 'Pitx3', 'Prok2', 'Irx4', 'Hoxd1', 'Zic2', 'Zic5', 'Six2', 'Pax6', 'Hoxb1', 'Pax5', 'Fgf8', 'Nfatc1', 'Gbf1', 'Hoxc8', 'Hoxc6', 'Nkx2-4', 'Pou4f1', 'Zic4')
##had to delete 3 as weren't in data
gene_order <- c('Hoxd9', 'Irx2', 'Irx1', 'Gbx2', 'Hoxa2', 'Prox1', 'Pax2', 'Sall3', 'Irx5', 'Hoxd3', 'Kazald1', 'Hoxa4', 'Irx3', 'Hoxb4', 'Hoxd4', 'Hoxa10', 'Hoxa6', 'En2', 'Dach1', 'Lmx1b', 'Nr2f2', 'Hoxc9', 'Hoxb9', 'Eif4e3', 'Cntnap1', 'Coasy', 'Hoxb2', 'Hoxa9', 'Hoxd10', 'Hoxa5', 'Hoxd11', 'Hoxa7', 'Hoxc4', 'Hoxa3', 'Hoxb5', 'Hoxb6', 'Hoxb8', 'Pitx2', 'Shh', 'Hoxb7', 'Foxa1', 'Lmo1', 'Htr7', 'Slc22a4', 'Hoxc5', 'Hoxc10', 'Hoxa1', 'Maf', 'Sox6', 'Ebf1', 'Pax9', 'Hoxa11', 'Foxd3', 'Dlx1', 'Nkx6-1', 'Zic1', 'Hoxd13', 'Pitx3', 'Prok2', 'Irx4', 'Hoxd1', 'Zic2', 'Zic5', 'Six2', 'Pax6', 'Hoxb1', 'Pax5', 'Nfatc1', 'Gbf1', 'Hoxc8', 'Hoxc6', 'Pou4f1', 'Zic4')
genes <- read.csv('dave_bivalent_genes_1118.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[gene_order,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
ph = pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, cluster_rows=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse.bilvalent_gene.acomy_order.pdf', width = 7, height = 10)
##take2...
##get mouse gene names in dataset
write.csv(assaygenes, file="all_mouse_genes.csv")
##lets do it again so we have no clustering shown
gene_order <- c('Hoxd9', 'Irx2', 'Irx1', 'Gbx2', 'Hoxa2', 'Slc22a4', 'Gbf1', 'Hoxc8', 'Hoxc5', 'Hoxc10', 'Hoxa1', 'Pitx3', 'Prok2', 'Irx4', 'Six2', 'Hoxb1', 'Hoxd1', 'Zic2', 'Zic5', 'Dlx1', 'Nkx6-1', 'Zic1', 'Hoxd13', 'Pax6', 'Pax5', 'Nfatc1', 'Pou4f1', 'Zic4', 'Pitx2', 'Foxd3', 'Pax9', 'Hoxa11', 'Maf', 'Sox6', 'Hoxa6', 'Dach1', 'Lmx1b', 'En2', 'Ebf1', 'Prox1', 'Pax2', 'Sall3', 'Hoxa4', 'Irx3', 'Hoxb4', 'Hoxd4', 'Hoxd3', 'Kazald1', 'Irx5', 'Nr2f2', 'Hoxb6', 'Hoxb8', 'Hoxc4', 'Hoxa3', 'Hoxb2', 'Hoxb5', 'Htr7', 'Foxa1', 'Lmo1', 'Hoxa5', 'Hoxd11', 'Cntnap1', 'Coasy', 'Hoxa10', 'Hoxa9', 'Hoxd10', 'Hoxc9', 'Hoxb9', 'Hoxa7', 'Shh', 'Hoxb7', 'Eif4e3', 'Hoxc6')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
newmat <- assay(rld)[gene_order,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
ph = pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, cluster_rows=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618.mouse.bilvalent_genes.110618.pdf', width = 7, height = 10)



##acomy data
##read in count and metadata
countData1 <- read.table('acomy_rnaseq_0618_acomy.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('acomy_rnaseq_0618_acomy.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ condition)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="acomy_rnaseq_0618_acomy.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="acomy_rnaseq_0618_acomy.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("condition"))
ggsave('acomy_rnaseq_0618_acomy.sample_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "day2"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:36029,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_acomy.sham_vs_day2.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "day5"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:36029,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_acomy.sham_vs_day5.csv")
##heatmap for dave - dna transcription
genes <- read.csv('DNA_templated_transcription_gene_list.acomy_id.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy.DNA_templated_transcription.pdf', width = 7, height = 10)
##heatmap for sam - myofibroblasts tfs
genes <- read.csv('sam_myo_tf.acomy_id.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=10, fontsize_col=10)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy.sam_myo_tf.pdf', width = 7, height = 10)
##heatmap for dave - bivalent genes
genes <- read.csv('dave_bivalent_genes_1118.acomy_id.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
ph <- pheatmap(newmat, annotation_col=newdf, fontsize_row=7, fontsize_col=7, cluster_cols=F)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy.bilvalent_gene.pdf', width = 7, height = 10)
##get the row order
ph_order = ph$tree_row$order
##get the corresponding gene names
gene_names = rownames(newmat)[ph_order]
write.csv(gene_names, file="acomy_gene_order.bilvalent.csv")
##now print in same order as acomy (got manually)
gene_order <- c('7642_g', '11196_g', '23238_g', '6853_g', '23233_g', '23236_g', '23232_g', '5083_g', '23227_g', '31282_g', '7638_g', '6068_g', '42422_g', '25672_g', '7627_g', '7923_g', '25669_g', '25667_g', '21791_g', '11199_g', '29413_g', '7628_g', '42177_g', '7632_g', '26323_g', '7934_g', '15793_g', '7630_g', '4271_g', '42231_g', '1180_g', '23235_g', '15797_g', '15798_g', '40597_g', '25675_g', '35698_g', '38120_g', '25756_g', '3393_g', '36395_g', '15072_g', '1194_g', '23225_g', '3396_g', '1188_g', '15424_g', '15209_g', '45795_g', '15313_g', '35697_g', '28017_g', '22973_g', '21778_g', '1882_g', '15264_g', '15789_g', '42216_g', '25668_g', '7634_g', '23231_g', '7636_g', '28828_g', '28838_g', '36044_g', '25670_g', '15791_g', '25665_g', '7639_g', '4265_g', '12741_g', '15795_g', '7641_g', '34216_g', '22984_g', '29252_g', '29698_g', '13355_g', '21806_g', '6863_g')
gene_order <- c('7642_g', '11196_g', '23238_g', '6853_g', '23233_g', '23236_g', '23232_g', '5083_g', '23227_g', '31282_g', '7638_g', '6068_g', '42422_g', '25672_g', '7627_g', '7923_g', '25669_g', '25667_g', '21791_g', '11199_g', '29413_g', '7628_g', '42177_g', '7632_g', '26323_g', '7934_g', '15793_g', '7630_g', '4271_g', '42231_g', '23235_g', '15797_g', '15798_g', '40597_g', '35698_g', '38120_g', '25756_g', '15072_g', '23225_g', '3396_g', '45795_g', '15313_g', '35697_g', '28017_g', '22973_g', '21778_g', '1882_g', '15264_g', '15789_g', '42216_g', '25668_g', '7634_g', '23231_g', '7636_g', '28828_g', '28838_g', '36044_g', '25670_g', '15791_g', '25665_g', '7639_g', '4265_g', '15795_g', '7641_g', '34216_g', '22984_g', '29252_g', '29698_g', '13355_g', '21806_g', '6863_g')
#genes <- read.csv('dave_bivalent_genes_1118.acomy_id.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
#genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
#idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[gene_order,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
ph = pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, cluster_rows=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy.bilvalent_gene.mouse_order.pdf', width = 7, height = 10)
##write 
gene_order <- c('25669_g', '21791_g', '21806_g', '6068_g', '7641_g', '13355_g', '42231_g', '6863_g', '28838_g', '25667_g', '42216_g', '7638_g', '28828_g', '23231_g', '25668_g', '7628_g', '7634_g', '7934_g', '15313_g', '26326_g', '28017_g', '15797_g', '23238_g', '11196_g', '22973_g', '22984_g', '23227_g', '7630_g', '25670_g', '7636_g', '25672_g', '7632_g', '15789_g', '7639_g', '23232_g', '23233_g', '23236_g', '31283_g', '7924_g', '23235_g', '4271_g', '29413_g', '42422_g', '34220_g', '15792_g', '15798_g', '7642_g', '29698_g', '29257_g', '36044_g', '4265_g', '7627_g', '38120_g', '25756_g', '43293_g', '35698_g', '15802_g', '45795_g', '11199_g', '21778_g', '25665_g', '15072_g', '15077_g', '3396_g', '1882_g', '23225_g', '5084_g', '6855_g', '42178_g', '15795_g', '15793_g', '15264_g', '35697_g')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
newmat <- assay(rld)[gene_order,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
ph = pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618.acomy.bilvalent_genes.110618.pdf', width = 7, height = 10)
##get the row order
ph_order = ph$tree_row$order
##get the corresponding gene names
gene_names = rownames(newmat)[ph_order]
write.csv(gene_names, file="acomy_gene_order.bilvalent.csv")
##lets do it again so we have no clustering shown
gene_order <- c('25669_g', '21791_g', '21806_g', '6068_g', '7641_g', '34220_g', '42178_g', '15795_g', '15792_g', '15798_g', '7642_g', '45795_g', '11199_g', '21778_g', '3396_g', '23225_g', '25665_g', '15072_g', '15077_g', '25756_g', '43293_g', '35698_g', '15802_g', '1882_g', '5084_g', '6855_g', '15264_g', '35697_g', '31283_g', '38120_g', '4265_g', '7627_g', '29698_g', '29257_g', '7634_g', '15313_g', '26326_g', '7934_g', '36044_g', '13355_g', '42231_g', '6863_g', '7638_g', '28828_g', '23231_g', '25668_g', '25667_g', '42216_g', '28838_g', '28017_g', '23233_g', '23236_g', '15789_g', '7639_g', '23227_g', '23232_g', '42422_g', '4271_g', '29413_g', '7636_g', '25672_g', '22973_g', '22984_g', '7628_g', '7630_g', '25670_g', '15797_g', '23238_g', '7632_g', '7924_g', '23235_g', '11196_g', '15793_g')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
newmat <- assay(rld)[gene_order,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
ph = pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=7, cluster_cols=F, cluster_rows=F, show_rownames=F)
dev.copy2pdf(file='acomy_rnaseq_0618.acomy.bilvalent_genes.110618.pdf', width = 7, height = 10)

##looking at 2 samples that looked switched
##rnemaed to not_known
##read in count and metadata
countData1 <- read.table('acomy_rnaseq_0618_acomy.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('acomy_rnaseq_0618_acomy_test.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ condition)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##rlog transform
rld <- rlog(dds, blind=FALSE)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "day2", "day5"))
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:36029,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_acomy_test.day2_vs_day5.csv")
##take genes from gene list
####from gene list
genes <- read.csv('top20_d2_d5.csv', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(newmat, annotation_col=newdf, fontsize_row=10, fontsize_col=10)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy_test.20_genes_d2_vs_d5.gene_clustering.pdf', width = 7, height = 10)




##mouse data for 'treatment' analysis i.e. two groups sham and treatment
##read in count and metadata
countData1 <- read.table('acomy_rnaseq_0618_mouse_treatment.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('acomy_rnaseq_0618_mouse_treatment.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ condition)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="acomy_rnaseq_0618_mouse_treatment.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="acomy_rnaseq_0618_mouse_treatment.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse_treatment.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("condition"))
ggsave('acomy_rnaseq_0618_mouse_treatment.sample_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='acomy_rnaseq_0618_mouse_treatment.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "treatment"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19100,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_mouse_treatment.sham_vs_treatment.csv")



##acomy data for 'treatment' analysis i.e. two groups sham and treatment
##read in count and metadata
countData1 <- read.table('acomy_rnaseq_0618_acomy_treatment.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('acomy_rnaseq_0618_acomy_treatment.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ condition)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="acomy_rnaseq_0618_acomy_treatment.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="acomy_rnaseq_0618_acomy_treatment.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy_treatment.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("condition"))
ggsave('acomy_rnaseq_0618_acomy_treatment.sample_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy_treatment.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "treatment"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:36029,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_acomy_treatment.sham_vs_treatment.csv")

##acomy data -- with 2 questionable samples removed
##read in count and metadata
countData1 <- read.table('acomy_rnaseq_0618_acomy_reduced.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('acomy_rnaseq_0618_acomy_reduced.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ condition)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)
##write normalized counts
#dds_norm <- estimateSizeFactors(dds)
#count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="acomy_rnaseq_0618_acomy_reduced.norm_counts.csv")
##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="acomy_rnaseq_0618_acomy_reduced.deseq.rlog_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy_reduced.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("condition"))
ggsave('acomy_rnaseq_0618_acomy_reduced.sample_pca.pdf') 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)["condition"])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='acomy_rnaseq_0618_acomy_reduced.25_var.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "day2"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:35368,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_acomy_reduced.sham_vs_day2.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "sham", "day5"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:35368,]
write.csv(resOrdered2DF, file="acomy_rnaseq_0618_acomy_reduced.sham_vs_day5.csv")
