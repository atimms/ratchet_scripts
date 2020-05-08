library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

setwd('/data/atimms/kim_splitseq_0318')

##load data -- trimmed set and complete
hu_data <- read.csv("digital_gene_expression.csv", header=T)
hu_data <- t(hu_data)
hu_data <- as.data.frame(hu_data)
head(hu_data)

##load data
seurat_mat <- CreateSeuratObject(raw.data = hu_data, min.cells = 3, min.genes = 20, 
                                 project = "kim_splitseq_0318")

##make sparse matrix
class(x = seurat_mat@raw.data)
seurat_mat <- MakeSparse(object = seurat_mat)
class(x = seurat_mat@raw.data)


# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT.", x = rownames(x = seurat_mat@data), value = TRUE)
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

##numbers to use?
seurat_mat <- FilterCells(object = seurat_mat, subset.names = c("nGene", "percent.mito"), 
                          low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.1))

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

##takes a long time - JackStraw gives significance to PCs
#pbmc <- JackStraw(object = seurat_mat, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(object = seurat_mat, PCs = 1:12)
##elbow plot to visulize pcs
PCElbowPlot(object = seurat_mat)

##Cluster the cells, using pcas from previous i.e. choose how many you want using in dims.use
seurat_mat <- FindClusters(object = seurat_mat, reduction.type = "pca", dims.use = 1:10, 
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)
##Run Non-linear dimensional reduction (tSNE)
seurat_mat <- RunTSNE(object = seurat_mat, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = seurat_mat, do.label = T)
##and repeat for the pdf
pdf("kim_splitseq_0318.tsne.all_data.pdf")
TSNEPlot(object = seurat_mat, do.label = T)
dev.off()
save(seurat_mat, file = "kim_splitseq_0318a.Robj")

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
all.markers <- FindAllMarkers(object = seurat_mat, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
write.csv(all.markers, file="splitseq_0318.all_cluster_markers.csv")

##induivdual test
#cluster0.markers <- FindMarkers(object = seurat_mat, ident.1 = 0, only.pos = TRUE, min.pct = 0.8, thresh.use = 0.5)
#head(cluster0.markers)

##vizulize genes
VlnPlot(object = seurat_mat, features.plot = c("FSTL4", "ERBB4"))
FeaturePlot(object = seurat_mat, features.plot = c("FSTL4", "ERBB4"), cols.use = c("grey", "blue"), reduction.use = "tsne")
##and print put



