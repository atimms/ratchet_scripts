library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

setwd('/data/atimms/cahan_sc_eb_1017')

##day6

##load data
hu_data <- read.csv("cahan_eb_d6.csv",header=T,row.names = 1)

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = hu_data))
dense.size
sparse.size <- object.size(x = hu_data)
sparse.size
seurat_mat <- CreateSeuratObject(raw.data = hu_data, min.cells = 3, min.genes = 400, 
                                 project = "cahan_eb_d6")

# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat_mat@data), value = TRUE)
percent.mito <- Matrix::colSums(seurat_mat@raw.data[mito.genes, ])/Matrix::colSums(seurat_mat@raw.data)

# filter by mito and umi content
#seurat_mat <- FilterCells(object = seurat_mat, subset.names = c("nGene", "percent.mito"), 
#                          low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.1))

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
#GenePlot(object = seurat_mat, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat_mat, gene1 = "nUMI", gene2 = "nGene")

##Normalizing the data
seurat_mat <- NormalizeData(object = seurat_mat, normalization.method = "LogNormalize", 
                            scale.factor = 10000)

##Detection of variable genes across the single cells
seurat_mat <- FindVariableGenes(object = seurat_mat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = seurat_mat@var.genes)
##Scaling the data and removing unwanted sources of variation
seurat_mat <- ScaleData(object = seurat_mat, vars.to.regress = c("nUMI", "percent.mito"))
##run PC analysis
seurat_mat <- RunPCA(object = seurat_mat, pc.genes = seurat_mat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                     genes.print = 5)
# Examine and visualize PCA results a few different ways
PrintPCA(object = seurat_mat, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
PCAPlot(object = seurat_mat, dim.1 = 1, dim.2 = 2)
VizPCA(object = seurat_mat, pcs.use = 1:2)
##elbow plot to decide with PCs to use
PCElbowPlot(object = seurat_mat)
##cluster the cells
seurat_mat <- FindClusters(object = seurat_mat, reduction.type = "pca", dims.use = 1:10, 
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)
##run tsne
seurat_mat <- RunTSNE(object = seurat_mat, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = seurat_mat, do.label=T)



##so remove cluster 5 and reanlyze

##look at data then remove cluster 5
head(seurat_mat@ident)
sm_minus5 <- SubsetData(object = seurat_mat, ident.remove = 5)
##cluster the cells -- added force.recalc =true
sm_minus5 <- FindClusters(object = sm_minus5, reduction.type = "pca", dims.use = 1:10, 
                           resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
##run tsne
sm_minus5 <- RunTSNE(object = sm_minus5, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = sm_minus5, do.label=T)
##and print
pdf("cahan_d6_minus_cluster5_tSNE.pdf")
TSNEPlot(sm_minus5, do.label=T)
dev.off()

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
all.markers <- FindAllMarkers(object = seurat_mat, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
write.csv(all.markers, file="cahan_d6_minus_cluster5.all_cluster_markers.csv")


