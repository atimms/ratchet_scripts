# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA')

# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)

#Import 10x data and convert each dataset to Seurat object (take from individual sample outputs, not cellranger's aggr output). Define sample with project= "samplename"
d53.data = Read10X(data.dir = './d53/outs/filtered_feature_bc_matrix/')
d53 = CreateSeuratObject(counts = d53.data, project = "d53", min.cells = 3, min.features = 200)

human = d53
#Store mitochondrial percentage in the Seurat object metadata
human[["percent.mt"]] <- PercentageFeatureSet(human, pattern = "^MT-")

#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.qc.pdf", width=20)

#Filter the data
human_f1 <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)
human_f2 <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 12)

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human_f1 = NormalizeData(human_f1)
human_f1 <- FindVariableFeatures(human_f1, selection.method = "vst", nfeatures = 2000)
human_f1 = ScaleData(human_f1)
human_f1 = RunPCA(human_f1)
DimPlot(human_f1, reduction='pca')
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.pca.pdf", width=20)
ElbowPlot(human_f1, ndims=30)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.elbow_plot.pdf", width=20)

#Run SCTransform, set var.to.regress to percent.mt
human_f1 <- SCTransform(human_f1, vars.to.regress = "percent.mt", verbose = FALSE)

#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters).
##change parameters --20 dims and 0.8 resolution
human_f1 <- RunPCA(human_f1, verbose = FALSE)
human_f1 <- RunUMAP(human_f1, dims = 1:20, verbose = FALSE)
human_f1 <- FindNeighbors(human_f1, dims = 1:20, verbose=FALSE)
human_f1 <- FindClusters(human_f1, resolution = 0.8)
#Plot UMAP
DimPlot(human_f1, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.UMAP.20dim_res0.8.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f1$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.counts_cluster.20dim_res0.8.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.markers.20dim_res0.8.csv')

##change parameters --20 dims and 0.4 resolution
human_f1 <- RunPCA(human_f1, verbose = FALSE)
human_f1 <- RunUMAP(human_f1, dims = 1:20, verbose = FALSE)
human_f1 <- FindNeighbors(human_f1, dims = 1:20, verbose=FALSE)
human_f1 <- FindClusters(human_f1, resolution = 0.4)
#Plot UMAP
DimPlot(human_f1, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.UMAP.20dim_res0.4.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f1$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.counts_cluster.20dim_res0.4.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.markers.20dim_res0.4.csv')

##change parameters --30 dims and 0.8 resolution
human_f1 <- RunPCA(human_f1, verbose = FALSE)
human_f1 <- RunUMAP(human_f1, dims = 1:30, verbose = FALSE)
human_f1 <- FindNeighbors(human_f1, dims = 1:30, verbose=FALSE)
human_f1 <- FindClusters(human_f1, resolution = 0.8)
#Plot UMAP
DimPlot(human_f1, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.UMAP.30dim_res0.8.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f1$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.counts_cluster.30dim_res0.8.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.markers.30dim_res0.8.csv')

##change parameters --30 dims and 0.4 resolution
human_f1 <- RunPCA(human_f1, verbose = FALSE)
human_f1 <- RunUMAP(human_f1, dims = 1:30, verbose = FALSE)
human_f1 <- FindNeighbors(human_f1, dims = 1:30, verbose=FALSE)
human_f1 <- FindClusters(human_f1, resolution = 0.4)
#Plot UMAP
DimPlot(human_f1, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.UMAP.30dim_res0.4.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f1$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.counts_cluster.30dim_res0.4.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f1.markers.30dim_res0.4.csv')

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human_f2 = NormalizeData(human_f2)
human_f2 <- FindVariableFeatures(human_f2, selection.method = "vst", nfeatures = 2000)
human_f2 = ScaleData(human_f2)
human_f2 = RunPCA(human_f2)
DimPlot(human_f2, reduction='pca')
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.pca.pdf", width=20)
ElbowPlot(human_f2, ndims=30)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.elbow_plot.pdf", width=20)

#Run SCTransform, set var.to.regress to percent.mt
human_f2 <- SCTransform(human_f2, vars.to.regress = "percent.mt", verbose = FALSE)

#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters).
##change parameters --20 dims and 0.8 resolution
human_f2 <- RunPCA(human_f2, verbose = FALSE)
human_f2 <- RunUMAP(human_f2, dims = 1:20, verbose = FALSE)
human_f2 <- FindNeighbors(human_f2, dims = 1:20, verbose=FALSE)
human_f2 <- FindClusters(human_f2, resolution = 0.8)
#Plot UMAP
DimPlot(human_f2, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.UMAP.20dim_res0.8.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.counts_cluster.20dim_res0.8.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.markers.20dim_res0.8.csv')

##change parameters --20 dims and 0.4 resolution
human_f2 <- RunPCA(human_f2, verbose = FALSE)
human_f2 <- RunUMAP(human_f2, dims = 1:20, verbose = FALSE)
human_f2 <- FindNeighbors(human_f2, dims = 1:20, verbose=FALSE)
human_f2 <- FindClusters(human_f2, resolution = 0.4)
#Plot UMAP
DimPlot(human_f2, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.UMAP.20dim_res0.4.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.counts_cluster.20dim_res0.4.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.markers.20dim_res0.4.csv')

##change parameters --30 dims and 0.8 resolution
human_f2 <- RunPCA(human_f2, verbose = FALSE)
human_f2 <- RunUMAP(human_f2, dims = 1:30, verbose = FALSE)
human_f2 <- FindNeighbors(human_f2, dims = 1:30, verbose=FALSE)
human_f2 <- FindClusters(human_f2, resolution = 0.8)
#Plot UMAP
DimPlot(human_f2, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.UMAP.30dim_res0.8.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.counts_cluster.30dim_res0.8.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.markers.30dim_res0.8.csv')

##change parameters --30 dims and 0.4 resolution
human_f2 <- RunPCA(human_f2, verbose = FALSE)
human_f2 <- RunUMAP(human_f2, dims = 1:30, verbose = FALSE)
human_f2 <- FindNeighbors(human_f2, dims = 1:30, verbose=FALSE)
human_f2 <- FindClusters(human_f2, resolution = 0.4)
#Plot UMAP
DimPlot(human_f2, reduction='umap', label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.UMAP.30dim_res0.4.pdf", width=10)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2$seurat_clusters)
write.csv(counts_cluster_sample, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.counts_cluster.30dim_res0.4.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_f2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d50s/human_scrnaseq_d50s_0520.f2.markers.30dim_res0.4.csv')



#Save object
saveRDS(human_f1, file = "./seurat_analysis/d50s/human_f1.rds")
#saveRDS(human_f2, file = "./seurat_analysis/d50s/human_f2.rds")

