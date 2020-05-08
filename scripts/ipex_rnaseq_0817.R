library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/data/atimms/ipex_rnaseq_0817";
setwd(workingDir);

##exp1
countData1 <- read.table('ipex_rnaseq_0817.exp1.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('ipex_rnaseq_0817.exp1.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ group)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##rlog transform
rld <- rlogTransformation(dds, blind=TRUE)
##50 most variable genes
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
write.csv(mat, file="ipex_0817.exp1.gene_clustering.50.csv")
df <- as.data.frame(colData(rld)["group"])
pheatmap(mat, annotation_col=df, fontsize_row=3, fontsize_col=5)
#dev.copy2pdf(file='ipex_0817.exp1.gene_clustering.50.pdf', width = 7, height = 5)
dev.copy2pdf(file='ipex_0817.exp1.gene_clustering.50.pdf')
##100 most variable genes
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),100)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
#write.csv(mat, file="ipex_0817.exp1.gene_clustering.100.csv")
df <- as.data.frame(colData(rld)["group"])
pheatmap(mat, annotation_col=df, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='ipex_0817.exp1.gene_clustering.100new.pdf', width = 7, height = 5)
####from gene list
genes <- read.csv('levings_31_genes.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newmat2 = newmat[genes,]
write.csv(newmat2, file="ipex_0817.exp1.gene_clustering.31_gene_sig.csv")
newdf <- as.data.frame(colData(rld)["group"])
pheatmap(newmat2, annotation_col=newdf, fontsize_row=5, fontsize_col=5, cluster_rows = F)
dev.copy2pdf(file='ipex_0817.exp1.gene_clustering.31_gene_sig.pdf', width = 7, height = 5)
####from gene list
genes <- read.csv('levings_31_genes_minus_il1rn.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newmat2 = newmat[genes,]
write.csv(newmat2, file="ipex_0817.exp1.gene_clustering.30_gene_sig.csv")
newdf <- as.data.frame(colData(rld)["group"])
pheatmap(newmat2, annotation_col=newdf, fontsize_row=5, fontsize_col=5, cluster_rows = F)
dev.copy2pdf(file='ipex_0817.exp1.gene_clustering.30_gene_sig.pdf', width = 7, height = 5)
##pca analysis
plotPCA(rld, intgroup = "group")
ggsave('ipex_0817.exp1.sample_pca.pdf') 

##130 most variable genes
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),130)
mat <- assay(rld)[ topVarGenes, ]
a <- row.names(mat)
write.csv(a, file="ipex_0817.top_130_var_genes.csv")

genes <- read.csv('ipex_topvar99_CD59_minus3_gl.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newmat2 = newmat[genes,]
newdf <- as.data.frame(colData(rld)["group"])
pheatmap(newmat2, annotation_col=newdf, fontsize_row=2, fontsize_col=5)
dev.copy2pdf(file='ipex_0817.exp1.topvar99_CD59_minus3.pdf', width = 7, height = 5)

