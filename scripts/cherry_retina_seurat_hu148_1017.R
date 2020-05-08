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

seurat_mat <- CreateSeuratObject(raw.data = hu_data, min.cells = 3, min.genes = 400, 
                           project = "retina_hu148")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
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

seurat_mat <- FilterCells(object = seurat_mat, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.1))

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
TSNEPlot(object = seurat_mat)
#save(seurat_mat, file = "~/Projects/datasets/pbmc3k/pbmc_tutorial.Robj")

cluster8.markers <- FindMarkers(object = seurat_mat, ident.1 = 8, min.pct = 0.25)
print(x = head(x = cluster8.markers, n = 5))

VlnPlot(object = seurat_mat, features.plot = c("NEFL", "NEFM"))

FeaturePlot(object = seurat_mat, features.plot = c("NEFL", "NEFM", "AHNAK2", "NEFH", 
                                             "UCHL1"), cols.use = c("grey", "blue"), reduction.use = "tsne")

##from gene list 050818

##and repeat for the pdf
pdf("cherry_hu148_0817.vln.genelist_050818.pdf")
VlnPlot(object = seurat_mat, features.plot = c("PHGDH", "CPS1", "PSPH", "ALDH1L1", "TMEM161B", "TMEM161-AS1", "MEF2C", "LINC00461", "MIR9-1", "MIR9-2", "MIR9-3", "SPTLC1", "SPTLC2", "SPTLC3", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6", "DES", "CERT", "SK1", "SK2", "NR2E3", "ARR3", "ONECUT2", "VSX2", "CLU", "GAD1", "NEFM", "RHO", "RP1", "PDE6H", "OPN1LW", "OPN1SW", "CABP5", "RASA1", "CCNH", "ALDH3A2"))
dev.off()

gene_list = c('MEF2C', 'TMEM161B', 'LINC00461', 'NR2E3', 'ARR3', 'ONECUT2', 'VSX2', 'CLU', 'GAD1', 'NEFM')
gene_list = c("PHGDH", "CPS1", "PSPH", "ALDH1L1", "TMEM161B", "MEF2C", "LINC00461", "SPTLC1", "SPTLC2", "SPTLC3", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6", "DES", "NR2E3", "ARR3", "ONECUT2", "VSX2", "CLU", "GAD1", "NEFM", "RHO", "RP1", "PDE6H", "OPN1LW", "OPN1SW", "CABP5", "RASA1", "CCNH", "ALDH3A2")
VlnPlot(object = seurat_mat, gene_list, point.size.use=NA, nCol=11)

