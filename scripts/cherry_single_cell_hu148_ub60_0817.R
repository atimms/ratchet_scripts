library(monocle)
library(ggplot2)
library(Matrix)
library(pheatmap)
library(reshape)

setwd('/data/atimms/cherry_single_cell_0717/hu148')

##input files
HSMM_expr_matrix <- read.table("cherry_hu148_0817.raw_martix.txt",header=T,row.names = 1)
##turn dataset to matrix and the make sparse matrix
#HSMM_expr_matrix_sparse <- as(as.matrix(HSMM_expr_matrix), "sparseMatrix")
HSMM_sample_sheet <- read.table("cherry_hu148_0817.sample_sheet.txt",header=T,row.names = 1)
HSMM_gene_annotation <- read.table("cherry_hu148_0817.gene_ann.txt",header=T,row.names = 1)
##create CellDataSet
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
##use sparse matrix ??use
#HSMM <- newCellDataSet(HSMM_expr_matrix_sparse, phenoData = pd, featureData = fd,
#                       expressionFamily=negbinomial.size())
##or use dense matrix
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                       expressionFamily=negbinomial.size())

##prep input data

##get number of genes expressed in 50 or more cells (what number to use???)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))
length(expressed_genes)
##look at data on cells, reduce data to those with > 400 genes expressed (just guessing here)
print(head(pData(HSMM)))
valid_cells <- row.names(subset(pData(HSMM),num_genes_expressed > 400))
length(valid_cells)
HSMM <- HSMM[,valid_cells]
HSMM = estimateSizeFactors(HSMM)
print(head(pData(HSMM)))

##useful info for later
HSMM = estimateDispersions(HSMM)




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
#cth <- addCellType(cth, "Muller_glia", classify_func = function(x) { x[RP1_id,] < 1 & x[APOE_id,] >= 1 })
cth <- addCellType(cth, "Muller_glia", classify_func = function(x) {x[APOE_id,] >= 1 })
cth <- addCellType(cth, "Amacrine_Cells", classify_func = function(x) { x[GAD1_id,] >= 1 })
cth <- addCellType(cth, "Bipolar_Cells", classify_func = function(x) { x[VSX2_id,] >= 1 })
cth <- addCellType(cth, "Cone_photoreceptors", classify_func = function(x) { x[ARR3_id,] >= 1 })
cth <- addCellType(cth, "Horizontal_Cells", classify_func = function(x) { x[LHX1_id,] >= 1 })
cth <- addCellType(cth, "Ganglion_Cell", classify_func = function(x) { x[NEFL_id,] >= 1 })
cth <- addCellType(cth, "Astrocytes", classify_func = function(x) { x[GFAP_id,] >= 1 })
cth <- addCellType(cth, "Microglia", classify_func = function(x) { x[HLA_DRA_id,] >= 1 })
cth <- addCellType(cth, "Vascular_endothelial", classify_func = function(x) { x[VWF_id,] >= 1 })
HSMM <- classifyCells(HSMM, cth, 0.1)
table(pData(HSMM)$CellType)
pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("cherry_hu148_0817.cell_type_pie.pdf", width = 6, height = 4)


##Clustering cells without marker genes
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.05)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F, maxit = 1000) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T,
                        residualModelFormulaStr = "~donor + num_genes_expressed")
HSMM <- clusterCells(HSMM)
#plot_cell_clusters(HSMM, 1, 2, markers = c("RP1", "APOE", "GAD1", "VSX2" ))
plot_cell_clusters(HSMM, 1, 2)
print(head(pData(HSMM)))
write.csv(pData(HSMM), file="cherry_hu148_0817.cell_info1.csv")
#second attempt
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F, maxit = 1000) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T,
                        residualModelFormulaStr = "~donor")
HSMM <- clusterCells(HSMM)
plot_cell_clusters(HSMM, 1, 2)
print(head(pData(HSMM)))
write.csv(pData(HSMM), file="cherry_hu148_0817.cell_info2.csv")


##not used atm

##Clustering cells using marker genes
marker_diff <- markerDiffTable(HSMM[expressed_genes,],
                               cth,
                               residualModelFormulaStr = "~donor + num_genes_expressed",
                               cores = 1)
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 10))
top_bins <- selectTopMarkers(marker_spec, 3)
write.csv(top_bins, file="cherry_hu148_0817.top_markers.csv")
HSMM = estimateDispersions(HSMM)
##not sure about this bit !!
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 100)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2, norm_method = 'log',
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~donor + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 12)
plot_cell_clusters(HSMM, 1, 2, color = "CellType")
##Imputing cell type
HSMM <- clusterCells(HSMM,
                     frequency_thresh = 0.1,
                     cell_type_hierarchy = cth)
plot_cell_clusters(HSMM, 1, 2, color = "CellType", markers = c("RP1", "APOE", "GAD1", "VSX2" ))

















HSMM <- clusterCells(HSMM,
                     num_clusters = 2,
                     frequency_thresh = 0.1,
                     cell_type_hierarchy = cth)












##Clustering cells without marker genes 
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType", markers = c("MYF5", "ANPEP"))



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
