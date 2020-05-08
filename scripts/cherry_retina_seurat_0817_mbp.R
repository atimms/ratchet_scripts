library(monocle)
library(Matrix)
library(ggplot2)

setwd('/Users/atimms/Desktop/ngs_data/cherry_single_cell_0717/seurat_de')


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











##other stuff
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "s_cluster")
plot_cell_clusters(HSMM, 1, 2, color = "donor")


cth <- newCellTypeHierarchy()
cds <- clusterCells(HSMM)

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("RP1", "RHO", "PDE6B")))
cds_subset <- HSMM[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~s_cluster")
diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_jitter(cds_subset, grouping = "s_cluster", color_by = "s_cluster",
                  nrow= 1, ncol = NULL, plot_trend = TRUE)
