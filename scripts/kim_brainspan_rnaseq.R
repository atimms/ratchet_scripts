library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/home/atimms/kim_brainspan";
setwd(workingDir);

##read in count and metadata
countData1 <- read.table('brainspan_rpkm_0916.txt', header=T, row.names=1)
colData1 <- read.table('brainspan_metadata_0916.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ age + structure)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)


##take genes from kim's list
####from gene list
#genes <- read.csv('65_ASD-Exome_GeneList.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- read.csv('kim_genelist_0916.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize

foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("tissue","individual")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
#dev.copy2pdf(file='dwm0116.deseq.gene_clustering.65_asd.pdf', width = 7, height = 10)
dev.copy2pdf(file='dwm0116.deseq.gene_clustering.kim_gl_0916.pdf', width = 7, height = 10)