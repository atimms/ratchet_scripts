library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

setwd('/data/atimms/cherry_single_cell_0717/seurat_analysis_1017')

##load data -- trimmed set and complete
#hu_data <- read.table("retina_seurat_0817.raw_martix.txt",header=T,row.names = 1)
hu_data <- read.table("cherry_hu148_0817.raw_martix.txt",header=T,row.names = 1)

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = hu_data))
dense.size
sparse.size <- object.size(x = hu_data)
sparse.size

hu148 <- CreateSeuratObject(raw.data = hu_data, min.cells = 3, min.genes = 400, 
                           project = "retina_hu148")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = hu148@data), value = TRUE)
percent.mito <- Matrix::colSums(hu148@raw.data[mito.genes, ])/Matrix::colSums(hu148@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
hu148 <- AddMetaData(object = hu148, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = hu148, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = hu148, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = hu148, gene1 = "nUMI", gene2 = "nGene")


##Normalizing the data
hu148 <- NormalizeData(object = hu148, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

##Detection of variable genes across the single cells
hu148 <- FindVariableGenes(object = hu148, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)

