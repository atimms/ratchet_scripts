library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

setwd('/data/atimms/dave_greb1l_sc_jun_1217')

load("All_mouse_eb.RData")
write.csv(df_cell, file = 'greb1l_1217.cell_info.csv')
write.csv(df_gene, file = 'greb1l_1217.gene_info.csv')
gene_count_df = as.data.frame(as.matrix(gene_count))

write.csv(gene_count_df, file = 'greb1l_1217.gene_count.csv')

##load data
seurat_mat <- CreateSeuratObject(raw.data = gene_count_df, min.cells = 3, min.genes = 400, 
                            project = "greb1l_1217")

# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat_mat@data), value = TRUE)
percent.mito <- Matrix::colSums(seurat_mat@raw.data[mito.genes, ])/Matrix::colSums(seurat_mat@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
seurat_mat <- AddMetaData(object = seurat_mat, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = seurat_mat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = seurat_mat, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat_mat, gene1 = "nUMI", gene2 = "nGene")


# filter by mito and umi content
#seurat_mat <- SubsetData(seurat_mat, subset.name = "percent.mito", accept.high = 0.10)
#seurat_mat <- SubsetData(seurat_mat, subset.name = "nUMI", accept.high = 15000)

##didn't use
#seurat_mat <- FilterCells(object = seurat_mat, subset.names = c("nGene", "percent.mito"), 
#                          low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.1))



##Normalizing the data
seurat_mat <- NormalizeData(object = seurat_mat, normalization.method = "LogNormalize", 
                            scale.factor = 10000)

##Detection of variable genes across the single cells
#seurat_mat <- MeanVarPlot(seurat_mat ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
seurat_mat <- FindVariableGenes(object = seurat_mat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = seurat_mat@var.genes)

##Scaling the data and removing unwanted sources of variation
seurat_mat <- ScaleData(object = seurat_mat, vars.to.regress = c("nUMI", "percent.mito"))

##Perform linear dimensional reduction
seurat_mat <- RunPCA(object = seurat_mat, pc.genes = seurat_mat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                     genes.print = 5)
# Examine and visualize PCA results a few different ways
PrintPCA(object = seurat_mat, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
PCAPlot(object = seurat_mat, dim.1 = 1, dim.2 = 2)
seurat_mat <- ProjectPCA(object = seurat_mat, do.print = FALSE)
PCHeatmap(object = seurat_mat, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seurat_mat, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
##elbow plot to visulize pcs
PCElbowPlot(object = seurat_mat)

##Cluster the cells, using pcas from previous i.e. choose how many you want using in dims.use
seurat_mat <- FindClusters(object = seurat_mat, reduction.type = "pca", dims.use = 1:10, 
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)
##Run Non-linear dimensional reduction (tSNE)
seurat_mat <- RunTSNE(object = seurat_mat, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = seurat_mat)
