library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

##for zswim
workingDir = "/home/atimms/chris_rnaseq";
setwd(workingDir);
analysis_list <- c("zswim6", "zswim6_and_8_a", "zswim6_and_8_b", "zswim_all")

##for piezo1
workingDir = "/home/atimms/piezo1_deseq";
setwd(workingDir);
analysis_list <- c("piezo1")

for (i in analysis_list) {

  htseq_counts <- paste(i, 'htseq.counts.txt',sep=".")
  metadata <- paste(i, 'metadata.txt',sep=".")
  ##read in count and metadata
  countData1 <- read.table(htseq_counts, header=T, row.names=1)
  colData1 <- read.table(metadata, header=T, row.names=1)

  ##add to deseq, give countdata and metadata and then design information i.e. info on sample types
  dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design =  ~ case_control)
  dds$case_control <- factor(dds$case_control, levels=c("control","case"))
  dds
  
  ##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
  nrow(dds)
  dds <- dds[ rowSums(counts(dds)) > 10, ]
  nrow(dds)

  ##rlog transform
  rld <- rlog(dds, blind=FALSE)
  ##check
  head(assay(rld), 3)
  head(assay(dds),3)

  ##principal components analysis 
  plotPCA(rld, intgroup = c("case_control"))
  #ggsave('lm.sample_pca.pdf') 

  
  ##calculate sample distances from rlog 
  sampleDists <- dist( t( assay(rld) ) )
  sampleDists
  ##and plot as heat map
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- (rld$case_control )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  heat_map <- paste(i, 'heatmap.pdf',sep=".")
  dev.copy2pdf(file=heat_map, width = 7, height = 5)
  
  ##differential expression
  
  ##do the test
  dds <- DESeq(dds)
  ##get results and summary
  (res <- results(dds))
  summary(res)
  ##save differentially expression results
  ##sort results by adjusted p-value
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  ##save results as dataframe and take top 1k results, then write csv file
  resOrderedDF <- as.data.frame(resOrdered)[1:1000,]
  top_1k <- paste(i, 'top_1k_de_genes.csv',sep=".")
  write.csv(resOrderedDF, file=top_1k)
  ##save all genes
  over_resDF <- as.data.frame(resOrdered)[1:35000,]
  all_de <- paste(i, 'all_de_genes.csv',sep=".")
  write.csv(over_resDF, file=all_de)
}
