# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA')

# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Make sure Harmony library is installed
#library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)

#Import 10x data and convert each dataset to Seurat object (take from individual sample outputs, not cellranger's aggr output). Define sample with project= "samplename"
d74.data = Read10X(data.dir = './d74/outs/filtered_feature_bc_matrix/')
d74= CreateSeuratObject(counts = d74.data, project = "d74", min.cells = 3, min.features = 200)
d78.data = Read10X(data.dir = './d78/outs/filtered_feature_bc_matrix/')
d78= CreateSeuratObject(counts = d78.data, project = "d78", min.cells = 3, min.features = 200)


# Merge into one single Seurat object
human=merge(d74, y=c(d78))
#Validate the merge by checking number of cells per group
table(human$orig.ident)

#Store mitochondrial percentage in the Seurat object metadata
human[["percent.mt"]] <- PercentageFeatureSet(human, pattern = "^MT-")

#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.qc.pdf", width=20)

#Add sample and condition information explicitly into the metadata (as oppsoed to storing in 'orig.ident') for future downstream analysis
#Add sample info (simply copying 'orig.ident')
human$sample <- plyr::mapvalues(
  x = human$orig.ident, 
  from = c('d74', 'd78'), 
  to = c('d74', 'd78')
)

#Validate new metadata columns by checking that number of cells per sample/phenotype adds up
table(human$sample)

#Run same QC metrics by new metadata columns to ensure it is the same as original QC metrics
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by='sample', ncol = 3, pt.size=0.1)

#Filter the data
human <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
human_f2 <- subset(human, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)
table(human$sample)
table(human_f2$sample)


#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human=NormalizeData(human)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)
human=ScaleData(human)
human = RunPCA(human)
DimPlot(human, reduction='pca', group.by='sample')
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.pca.pdf", width=20)
ElbowPlot(human, ndims=30)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.elbow_plot.pdf", width=20)

#If samples do not cluster together in PC space (which is the case here), then there is no need to run harmony (and likely no need to use the integrated method either); the merged analysis should do.
#Determine the number of dimensions to use in downstream analyses based on the point at which the Elbow Plot becomes flat (ok to be conservative)

#Save merged and filtered dataset as an R object to be able to re-load it for various future analyses without needing to perform the previous computations
saveRDS(human, file = "./seurat_analysis/d70s/human.rds")
##and load if needed
#human <- readRDS(file = "./seurat_analysis/d70s/human.rds")

#Run SCTransform, set var.to.regress to percent.mt
human_merged <- SCTransform(human, vars.to.regress = "percent.mt", verbose = FALSE)

##resolution - default i.e. 0.8
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
human_merged <- RunPCA(human_merged, verbose = FALSE)
human_merged <- RunUMAP(human_merged, dims = 1:30, verbose = FALSE)
human_merged <- FindNeighbors(human_merged, dims = 1:30, verbose=FALSE)
human_merged <- FindClusters(human_merged, resolution = 0.8)
#Plot UMAP
DimPlot(human_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.UMAP.res0.8.sample_split.pdf", width = 20)
DimPlot(human_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.UMAP.res0.8.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_merged$seurat_clusters, human_merged$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.counts_cluster_sample.res0.8.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.markers.res0.8.csv')

##resolution 0.4 and 30 dims
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
human_merged <- RunPCA(human_merged, verbose = FALSE)
human_merged <- RunUMAP(human_merged, dims = 1:30, verbose = FALSE)
human_merged <- FindNeighbors(human_merged, dims = 1:30, verbose=FALSE)
human_merged <- FindClusters(human_merged, resolution = 0.4)
#Plot UMAP
DimPlot(human_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.UMAP.res0.4.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_merged$seurat_clusters, human_merged$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.counts_cluster_sample.res0.4.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.markers.res0.4.csv')

##resolution 0.4 and 20 dims
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
human_merged <- RunPCA(human_merged, verbose = FALSE)
human_merged <- RunUMAP(human_merged, dims = 1:20, verbose = FALSE)
human_merged <- FindNeighbors(human_merged, dims = 1:20, verbose=FALSE)
human_merged <- FindClusters(human_merged, resolution = 0.4)
#Plot UMAP
DimPlot(human_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.UMAP.dim20_res0.4.sample_split.pdf", width = 20)
DimPlot(human_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.UMAP.dim20_res0.4.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_merged$seurat_clusters, human_merged$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.counts_cluster_sample.dim20_res0.4.csv')
#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.merged.markers.dim20_res0.4.csv')

#Save object
saveRDS(human_merged, file = "./seurat_analysis/d70s/human_merged.rds")


##harmony analysis...
##load if needed
human <- readRDS(file = "./seurat_analysis/d70s/human.rds")

#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
human_harmony <- RunHarmony(object = human, group.by.vars = 'sample')
DimPlot(human_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.harmony_plot.pdf", width = 20)

# Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
human_harmony <- RunUMAP(human_harmony, dims = 1:30, reduction = 'harmony')
human_harmony <- FindNeighbors(human_harmony, reduction = 'harmony', dims = 1:30)

##change resolutions - default i.e. 0.8
human_harmony <- FindClusters(human_harmony, resolution = 0.8)
#Plot UMAP
DimPlot(human_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.UMAP.res0.8.sample_split.pdf", width = 20)
DimPlot(human_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.UMAP.res0.8.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_harmony$seurat_clusters, human_harmony$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.counts_cluster_sample.res0.8.csv')
#Find all markers that define each cluster
human_harmony.markers <- FindAllMarkers(human_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_harmony.markers, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.markers.res0.8.csv')

##change resolutions - 0.4
human_harmony <- FindClusters(human_harmony, resolution = 0.4)
#Plot UMAP
DimPlot(human_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.UMAP.res0.4.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_harmony$seurat_clusters, human_harmony$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.counts_cluster_sample.res0.4.csv')
#Find all markers that define each cluster
human_harmony.markers <- FindAllMarkers(human_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_harmony.markers, file='./seurat_analysis/d70s/human_scrnaseq_d70s_0520.harmony.markers.res0.4.csv')

#Save object to avoid needing to re-run previous computations
saveRDS(human_harmony, file = "./seurat_analysis/d70s/human_harmony.rds")

##integrated analysis



##looking at filter2

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human_f2=NormalizeData(human_f2)
human_f2 <- FindVariableFeatures(human_f2, selection.method = "vst", nfeatures = 2000)
human_f2=ScaleData(human_f2)
human_f2 = RunPCA(human_f2)
DimPlot(human_f2, reduction='pca', group.by='sample')
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.pca.pdf", width=20)
ElbowPlot(human_f2, ndims=30)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.elbow_plot.pdf", width=20)

#If samples do not cluster together in PC space (which is the case here), then there is no need to run harmony (and likely no need to use the integrated method either); the merged analysis should do.
#Determine the number of dimensions to use in downstream analyses based on the point at which the Elbow Plot becomes flat (ok to be conservative)

#Save merged and filtered dataset as an R object to be able to re-load it for various future analyses without needing to perform the previous computations
saveRDS(human_f2, file = "./seurat_analysis/d70s/human_f2.rds")
##and load if needed
#human_f2 <- readRDS(file = "./seurat_analysis/d70s/human_f2.rds")

#Run SCTransform, set var.to.regress to percent.mt
human_f2_merged <- SCTransform(human_f2, vars.to.regress = "percent.mt", verbose = FALSE)

##resolution - default i.e. 0.8
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
human_f2_merged <- RunPCA(human_f2_merged, verbose = FALSE)
human_f2_merged <- RunUMAP(human_f2_merged, dims = 1:30, verbose = FALSE)
human_f2_merged <- FindNeighbors(human_f2_merged, dims = 1:30, verbose=FALSE)
human_f2_merged <- FindClusters(human_f2_merged, resolution = 0.8)
#Plot UMAP
DimPlot(human_f2_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.UMAP.res0.8.sample_split.pdf", width = 20)
DimPlot(human_f2_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.UMAP.res0.8.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2_merged$seurat_clusters, human_f2_merged$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.counts_cluster_sample.res0.8.csv')
#Find all markers that define each cluster
human_f2_merged.markers <- FindAllMarkers(human_f2_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_f2_merged.markers, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.markers.res0.8.csv')

##resolution 0.4 and 30 dims
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
human_f2_merged <- RunPCA(human_f2_merged, verbose = FALSE)
human_f2_merged <- RunUMAP(human_f2_merged, dims = 1:30, verbose = FALSE)
human_f2_merged <- FindNeighbors(human_f2_merged, dims = 1:30, verbose=FALSE)
human_f2_merged <- FindClusters(human_f2_merged, resolution = 0.4)
#Plot UMAP
DimPlot(human_f2_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human_f2_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.UMAP.res0.4.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2_merged$seurat_clusters, human_f2_merged$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.counts_cluster_sample.res0.4.csv')
#Find all markers that define each cluster
human_f2_merged.markers <- FindAllMarkers(human_f2_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_f2_merged.markers, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.merged.markers.res0.4.csv')

##harmony analysis...
##load if needed
#human_f2 <- readRDS(file = "./seurat_analysis/d70s/human_f2.rds")

#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
human_f2_harmony <- RunHarmony(object = human_f2, group.by.vars = 'sample')
DimPlot(human_f2_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.harmony_plot.pdf", width = 20)

# Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
human_f2_harmony <- RunUMAP(human_f2_harmony, dims = 1:30, reduction = 'harmony')
human_f2_harmony <- FindNeighbors(human_f2_harmony, reduction = 'harmony', dims = 1:30)

##change resolutions - default i.e. 0.8
human_f2_harmony <- FindClusters(human_f2_harmony, resolution = 0.8)
#Plot UMAP
DimPlot(human_f2_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.UMAP.res0.8.sample_split.pdf", width = 20)
DimPlot(human_f2_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.UMAP.res0.8.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2_harmony$seurat_clusters, human_f2_harmony$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.counts_cluster_sample.res0.8.csv')
#Find all markers that define each cluster
human_f2_harmony.markers <- FindAllMarkers(human_f2_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_f2_harmony.markers, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.markers.res0.8.csv')

##change resolutions - 0.4
human_f2_harmony <- FindClusters(human_f2_harmony, resolution = 0.4)
#Plot UMAP
DimPlot(human_f2_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human_f2_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.UMAP.res0.4.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_f2_harmony$seurat_clusters, human_f2_harmony$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.counts_cluster_sample.res0.4.csv')
#Find all markers that define each cluster
human_f2_harmony.markers <- FindAllMarkers(human_f2_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_f2_harmony.markers, file='./seurat_analysis/d70s/human_f2_scrnaseq_d70s_0520.harmony.markers.res0.4.csv')

#Save object to avoid needing to re-run previous computations
saveRDS(human_f2_harmony, file = "./seurat_analysis/d70s/human_f2_harmony.rds")

