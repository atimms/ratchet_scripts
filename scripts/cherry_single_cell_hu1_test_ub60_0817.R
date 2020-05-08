library(monocle)
library(ggplot2)
library(Matrix)
library(pheatmap)
library(reshape)

setwd('/data/atimms/cherry_single_cell_0717/hu1_test')

##input files
HSMM_expr_matrix <- read.table("Hu1_retina_1.raw_martix.txt",header=T,row.names = 1)
##turn dataset to matrix and the make sparse matrix
#HSMM_expr_matrix_sparse <- as(as.matrix(HSMM_expr_matrix), "sparseMatrix")
HSMM_sample_sheet <- read.table("Hu1_retina_1.sample_sheet.txt",header=T,row.names = 1)
HSMM_gene_annotation <- read.table("Hu1_retina_1.gene_ann.txt",header=T,row.names = 1)
##create CellDataSet
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
##use sparse matrix ??use
#HSMM <- newCellDataSet(HSMM_expr_matrix_sparse, phenoData = pd, featureData = fd,
#                       expressionFamily=negbinomial.size())
##or use dense matrix
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                       expressionFamily=negbinomial.size())
##useful info for later
HSMM = estimateSizeFactors(HSMM)
HSMM = estimateDispersions(HSMM)

##get number of genes expressed in 50 or more cells (what number to use???)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes)
##look at data on cells, reduce data to those with > 400 genes expressed (just guessing here)
print(head(pData(HSMM)))
valid_cells <- row.names(subset(pData(HSMM),num_genes_expressed > 400))
length(valid_cells)
HSMM <- HSMM[,valid_cells]
##classifying cells by type
RP1_id <- row.names(subset(fData(HSMM), gene_short_name == "RP1"))
APOE_id <- row.names(subset(fData(HSMM), gene_short_name == "APOE"))
GAD1_id <- row.names(subset(fData(HSMM), gene_short_name == "GAD1"))
VSX2_id <- row.names(subset(fData(HSMM), gene_short_name == "VSX2"))
ARR3_id <- row.names(subset(fData(HSMM), gene_short_name == "ARR3"))
LHX1_id <- row.names(subset(fData(HSMM), gene_short_name == "LHX1"))
NEFL_id <- row.names(subset(fData(HSMM), gene_short_name == "NEFL"))
GFAP_id <- row.names(subset(fData(HSMM), gene_short_name == "GFAP"))
HLA_DRA_id <- row.names(subset(fData(HSMM), gene_short_name == "HLA-DRA"))
VWF_id <- row.names(subset(fData(HSMM), gene_short_name == "VWF"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Rod_photoreceptors", classify_func = function(x) { x[RP1_id,] >= 1 })
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.05)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components = 15, num_dim = 15,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM)
plot_cell_clusters(HSMM, 1, 2,3, markers = c("RP1", "APOE", "GAD1", "VSX2" ))













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
