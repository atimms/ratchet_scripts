library(WGCNA);
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")


workingDir = "/home/atimms/chris_rnaseq"
setwd(workingDir)

##load data
countData1 <- read.table('css_rnaseq_thyroid_0816.htseq_counts.txt', header=T, row.names=1)
colData1 <- read.table('css_rnaseq_thyroid_0816.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)


##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ group)
dds$group <- factor(dds$group, levels=c("control","case"))
dds

##remove rows of the DESeqDataSet that have no counts, or less than x counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 50, ]
nrow(dds)

##write normalized counts
dds <- estimateSizeFactors(dds)
count_data <- counts(dds, normalized=TRUE)
write.csv(count_data, file="css_rnaseq.norm_counts.csv")



##rlog transform
rld <- rlog(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)

##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- (rld$group )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dev.copy2pdf(file='css_rnaseq_thyroid_0816.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
plotPCA(rld, intgroup = "group")
#ggsave('ipex.sample_pca.pdf') 
ggsave('css_rnaseq_thyroid_0816.sample_pca.pdf') 

##multidimensional scaling plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=group)) + geom_point(size=3)
#ggsave('ipex.sample_mds.pdf') 
ggsave('css_rnaseq_thyroid_0816.sample_mds.pdf') 



##gene clustering
##take 100 most variable genes
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),100)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("group", "sample_number")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#dev.copy2pdf(file='ipex.deseq.gene_clustering.25.pdf', width = 7, height = 5)


####from gene list
genes <- read.csv('65_ASD-Exome_GeneList.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize

genes <- c('AADAT', 'ATXN2', 'BACH2', 'BANK1', 'BTNL2', 'C1QTNF6', 'CAPZB', 'CD40', 'CGA', 'CLIC6', 'CTLA4', 'DAAM1', 'DIO1', 'DIO2', 'DIO3', 'DKK3', 'DUOX1', 'DUOX2', 'DUOXA1', 'DUOXA2', 'EPHA4', 'FBXO15', 'FGF21', 'FGF23', 'FGF7', 'FGFR1', 'FOXE1', 'FZD7', 'GLI2', 'GLIS3', 'GNAQ', 'GRPEL1', 'HACE1', 'HLX', 'IKBKG', 'IL16', 'IL1B', 'INSR', 'ITPK1', 'IYD', 'KALRN', 'LPCAT2', 'MAF', 'MAGI3', 'MBIP', 'MSX1', 'NFIA', 'NFKB1', 'NLRP1', 'NOTCH4', 'NOX4', 'NR3C2', 'NRG1', 'PARP1', 'PAX8', 'PDE10A', 'PDE4D', 'PDE8B', 'PITX1', 'PITX2', 'POR', 'PPP1R9A', 'PRDM11', 'PROP1', 'PTPN22', 'PTPRB', 'RGS12', 'RNASET2', 'SASH1', 'SERPINA7', 'SFRP4', 'SGK1', 'SLA', 'SLAMF6', 'SLC16A10', 'SLC16A2', 'SLC26A4', 'SLC5A5', 'SOCS3', 'SOX9', 'STX16', 'TBX21', 'TG', 'TGFB3', 'THRA', 'THRAP3', 'THRB', 'TNFAIP3', 'TPO', 'TRHR', 'TSHB', 'TSHR', 'TTC21A', 'TTF1', 'VEGFA', 'WFS1', 'WNT5A', 'ZDHHC21', 'ZFAT', 'ZIC2')

foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("group","sample_number")])
pheatmap(newmat, annotation_col=newdf)
dev.copy2pdf(file='css_rnaseq_thyroid_0816.all_samples.100_gene_clustering.pdf', width = 30, height = 40)


