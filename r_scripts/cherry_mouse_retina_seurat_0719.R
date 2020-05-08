library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

setwd('/data/atimms/cherry_single_cell_0717/mouse_retina_sc_1017')

##need this to use umap
#reticulate::py_install(packages ='umap-learn')

##load data
mouse_data <- read.table("GSE63472_P14Retina_merged_digital_expression.txt",header=T,row.names = 1)

# Initialize the Seurat object with the raw (non-normalized data).
seurat_mat <- CreateSeuratObject(counts = mouse_data, project = "mouse_retina_GSE63472", min.cells = 3, min.features = 200)
seurat_mat


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat_mat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_mat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

##normalize data
seurat_mat <- NormalizeData(seurat_mat)

##find variable genes, default is top 2k
seurat_mat <- FindVariableFeatures(seurat_mat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_mat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_mat)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))
# Scaling the data
all.genes <- rownames(seurat_mat)
seurat_mat <- ScaleData(seurat_mat, features = all.genes)
#Perform linear dimensional reduction
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(object = seurat_mat))
# Examine and visualize PCA results a few different ways
print(seurat_mat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat_mat, dims = 1:2, reduction = "pca")
DimPlot(seurat_mat, reduction = "pca")

#ranking of principle components based on the percentage of variance explained by each one
ElbowPlot(seurat_mat)

#Cluster the cells, using 10 pcs in this case
seurat_mat <- FindNeighbors(seurat_mat, dims = 1:10)
seurat_mat <- FindClusters(seurat_mat, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(seurat_mat), 5)


#Run non-linear dimensional reduction (UMAP/tSNE)
seurat_mat <- RunUMAP(seurat_mat, dims = 1:10)
DimPlot(seurat_mat, reduction = "umap")
dev.copy2pdf(file="cherry_mouse_seurat_0719.umap.pdf", width = 10)
DimPlot(seurat_mat, reduction = "umap", group.by = "orig.ident")
dev.copy2pdf(file="cherry_mouse_seurat_0719.umap.sample_id.pdf", width = 10)
seurat_mat <- RunTSNE(seurat_mat, dims = 1:10)
DimPlot(seurat_mat, reduction = "tsne")
dev.copy2pdf(file="cherry_mouse_seurat_0719.tsne.pdf", width = 10)
DimPlot(seurat_mat, reduction = "tsne", group.by = "orig.ident")
dev.copy2pdf(file="cherry_mouse_seurat_0719.tsne.sample_id.pdf", width = 10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_mat.markers <- FindAllMarkers(seurat_mat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(seurat_mat.markers, file="cherry_mouse_seurat_0719.cluster_markers.csv")

##save file
saveRDS(seurat_mat, file = "mouse_retina_sc_0719.rds")

##and then load
seurat_mat <- readRDS(file ="mouse_retina_sc_0719.rds")

# violin plots
VlnPlot(seurat_mat, features = c("HTRA1", "APOE"))
dev.copy2pdf(file="cherry_mouse_seurat_0719.violin_plot.HTRA1_APOE.pdf")
VlnPlot(seurat_mat, features = c("HTRA1", "APOE", "ARMS2", "PLEKHA1", "ONECUT2", "OPN1LW"))
dev.copy2pdf(file="cherry_mouse_seurat_0719.violin_plot.genes_070919.pdf")
##plot genes on tsne/umap
FeaturePlot(seurat_mat, features = c("HTRA1", "APOE"))
dev.copy2pdf(file="cherry_mouse_seurat_0719.feature_plot.HTRA1_APOE.pdf", width = 10, height = 5)
FeaturePlot(seurat_mat, features = c("HTRA1", "APOE", "ARMS2", "PLEKHA1", "ONECUT2", "OPN1LW"))
dev.copy2pdf(file="cherry_mouse_seurat_0719.feature_plot.genes_070919.pdf", width = 10, height = 15)








