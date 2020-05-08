library(monocle)
library(Matrix)
library(ggplot2)

setwd('/data/atimms/cherry_single_cell_0717/seurat_de')


##all samples and data
HSMM_expr_matrix <- read.table("retina_seurat_0817.raw_martix.txt",header=T,row.names = 1)
HSMM_sample_sheet <- read.table("retina_seurat_0817.sample_sheet.txt",header=T,row.names = 1)
HSMM_gene_annotation <- read.table("retina_seurat_0817.gene_ann.txt",header=T,row.names = 1)
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                       expressionFamily=negbinomial.size())
HSMM = estimateSizeFactors(HSMM)
HSMM = estimateDispersions(HSMM)
##diff expression, just in seurat cluster
res = differentialGeneTest(HSMM, fullModelFormulaStr='~s_cluster', reducedModelFormulaStr ='~donor' )
write.csv(res, file = "retina_seurat_0817.all_data.s_cluster.de.csv")

###cluster0 vs 6 - single diff test
HSMM_expr_matrix <- read.table("retina_seurat_0817.0_6.matrix.txt",header=T,row.names = 1)
HSMM_sample_sheet <- read.table("retina_seurat_0817.0_6.sample_sheet.txt",header=T,row.names = 1)
HSMM_gene_annotation <- read.table("retina_seurat_0817.gene_ann.txt",header=T,row.names = 1)
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                       expressionFamily=negbinomial.size())
HSMM = estimateSizeFactors(HSMM)
HSMM = estimateDispersions(HSMM)
##diff expression, just in seurat cluster
res = differentialGeneTest(HSMM, fullModelFormulaStr='~s_cluster')
write.csv(res, file = "retina_seurat_0817.0_6.s_cluster.de.csv")





##assess differential expression for multiple tests - no covariate
##test names
tests = c("0_6","0_1", "0_3", "0_4")
head(tests)
##run diff expression
for(i in 1:length(tests)){
  matrix_file <- paste("retina_seurat_0817",tests[i],"matrix.txt",sep=".")
  sample_sheet <- paste("retina_seurat_0817",tests[i],"sample_sheet.txt",sep=".")
  HSMM_expr_matrix <- read.table(matrix_file,header=T,row.names = 1)
  HSMM_sample_sheet <- read.table(sample_sheet,header=T,row.names = 1)
  HSMM_gene_annotation <- read.table("retina_seurat_0817.gene_ann.txt",header=T,row.names = 1)
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                         expressionFamily=negbinomial.size())
  HSMM = estimateSizeFactors(HSMM)
  HSMM = estimateDispersions(HSMM)
  ##diff expression, just in seurat cluster
  res = differentialGeneTest(HSMM, fullModelFormulaStr='~s_cluster', reducedModelFormulaStr ='~donor' )
  write.csv(res, file = "retina_seurat_0817.0_6.s_cluster.de.csv")
  }

##assess differential expression for multiple tests - sample as covariate

##one cluster vs another
##test names
tests = c("0_6","0_1", "0_3", "0_4")
head(tests)
##run diff expression
for(i in 1:length(tests)){
  matrix_file <- paste("retina_seurat_0817",tests[i],"matrix.txt",sep=".")
  sample_sheet <- paste("retina_seurat_0817",tests[i],"sample_sheet.txt",sep=".")
  HSMM_expr_matrix <- read.table(matrix_file,header=T,row.names = 1)
  HSMM_sample_sheet <- read.table(sample_sheet,header=T,row.names = 1)
  HSMM_gene_annotation <- read.table("retina_seurat_0817.gene_ann.txt",header=T,row.names = 1)
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                         expressionFamily=negbinomial.size())
  HSMM = estimateSizeFactors(HSMM)
  HSMM = estimateDispersions(HSMM)
  ##diff expression, just in seurat cluster
  res = differentialGeneTest(HSMM, fullModelFormulaStr='~s_cluster + donor', reducedModelFormulaStr ='~donor' )
  de_file <- paste("retina_seurat_0817",tests[i],"s_cluster.de_with_hu.csv",sep=".")
  write.csv(res, file = de_file)
}

##one cluster vs the rest
##test names
tests = c("0_3.vs_rest","1.vs_rest", "10.vs_rest", "12.vs_rest", "13.vs_rest", "2_9_11.vs_rest", "4.vs_rest", "5.vs_rest", "6.vs_rest", "7.vs_rest", "8.vs_rest")
head(tests)
##run diff expression
for(i in 1:length(tests)){
  sample_sheet <- paste("retina_seurat_0817",tests[i],"sample_sheet.txt",sep=".")
  HSMM_expr_matrix <- read.table("retina_seurat_0817.raw_martix.txt",header=T,row.names = 1)
  HSMM_sample_sheet <- read.table(sample_sheet,header=T,row.names = 1)
  HSMM_gene_annotation <- read.table("retina_seurat_0817.gene_ann.txt",header=T,row.names = 1)
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                         expressionFamily=negbinomial.size())
  HSMM = estimateSizeFactors(HSMM)
  HSMM = estimateDispersions(HSMM)
  ##diff expression, just in seurat cluster
  res = differentialGeneTest(HSMM, fullModelFormulaStr='~in_cluster + donor', reducedModelFormulaStr ='~donor' )
  de_file <- paste("retina_seurat_0817",tests[i],"s_cluster.de_with_hu.csv",sep=".")
  write.csv(res, file = de_file)
}

##get average expression value
## extract expression values
cds_exprs <- exprs(HSMM)
## transpose
cds_exprs <- Matrix::t(cds_exprs)
## compute relative expression?
##cds_exprs <- Matrix::t(cds_exprs)/sizeFactors(HSMM)
## sample info
cds_pData <- pData(HSMM)
## gene names
cds_fData <- fData(HSMM)
## add sample info to expression matrix
cds_exprs = merge(cds_exprs, cds_pData, by="row.names")
## compute mean expression for each gene, grouped by cluster
cds_exprs_mean = aggregate(cds_exprs[, row.names(cds_fData)],
                           list(cds_exprs$s_cluster), mean)
## transpose and then write out file
cds_exprs_mean <- Matrix::t(cds_exprs_mean)
dimnames(cds_exprs_mean)[[2]] = cds_exprs_mean[1,]

cds_exprs_mean = cds_exprs_mean[2:nrow(cds_exprs_mean),]
cds_exprs_mean_char = cds_exprs_mean
cds_exprs_mean = as.numeric(cds_exprs_mean)
attributes(cds_exprs_mean) = attributes(cds_exprs_mean_char)
write.csv(cds_exprs_mean, file = "retina_seurat_1017.mean_expression_per_cluster.csv")

##then make heatmap 
genes <- read.csv('genes_for_heatmap_1017.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
assaygenes <- row.names(cds_exprs_mean) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- cds_exprs_mean[idx,] #subset to target genes
newmat <- newmat + 0.01
newmat <- log2(newmat)
newmat <- newmat - rowMeans(newmat)
##rearrange file
sample_order = c("sc0", "sc3", "sc4", "sc1", "sc2", "sc9", "sc11", "sc5", "sc6", "sc7", "sc8", "sc10", "sc13")
reorder_newmat <- newmat[genes,sample_order]
#newdf <- as.data.frame(cds_exprs_mean)
##cluster genes and sample
pheatmap(newmat, fontsize_row=3, fontsize_col=5)
dev.copy2pdf(file='retina_seurat_1017.heatmap_clustered.pdf', width = 7, height = 5)
##no clustering
pheatmap(reodered_newmat, cluster_rows=F, cluster_cols=F, show_rownames = F)
dev.copy2pdf(file='retina_seurat_1017.heatmap.pdf', width = 7, height = 5)



