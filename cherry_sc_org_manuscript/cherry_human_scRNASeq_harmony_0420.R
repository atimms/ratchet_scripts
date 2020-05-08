# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Make sure Harmony library is installed
#library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)

# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA')

##if data isn't loaded
human_precheck <- readRDS(file = "./seurat_analysis/human_precheck.rds")

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)


#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
human_harmony <- RunHarmony(object = human_precheck, group.by.vars = 'sample')
DimPlot(human_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.harmony.harmony_plot.pdf", width = 20)

# Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
human_harmony <- RunUMAP(human_harmony, dims = 1:30, reduction = 'harmony')
human_harmony <- FindNeighbors(human_harmony, reduction = 'harmony', dims = 1:30)
#human_harmony <- FindClusters(human_harmony, resolution = 0.4)

##change resolutions - default
human_harmony <- FindClusters(human_harmony, group.singletons = TRUE)
#Plot UMAP
DimPlot(human_harmony, reduction='umap', split.by='sample', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.harmony.UMAP_res0.8.sample_split.pdf", width = 20)
DimPlot(human_harmony, reduction='umap', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.harmony.UMAP_res0.8.pdf", width=20)
DimPlot(human_harmony, reduction='umap', group.by = "sample", pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.harmony.UMAP_res0.8.sample.pdf", width=20)

#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_harmony$seurat_clusters, human_harmony$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/human_scrnaseq_0420.harmony.counts_cluster_sample.res0.8.csv')
#Identify the number of cells in each cluster between genotypes
counts_cluster_time =  table(human_harmony$seurat_clusters, human_harmony$time)
write.csv(counts_cluster_sample, file='./seurat_analysis/human_scrnaseq_0420.harmony.counts_cluster_time.res0.8.csv')

#Find all markers that define each cluster
human_harmony.markers <- FindAllMarkers(human_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_harmony.markers, file='./seurat_analysis/human_scrnaseq_0420.harmony.markers.res0.8.csv')

#Perform differential expression
human_harmony_time = human_harmony
Idents(human_harmony_time) <- 'time'
human_harmony.diffexp <- FindAllMarkers(human_harmony_time, min.pct=0.1, logfc.threshold = 0.25)
write.csv(human_harmony.diffexp, file='./seurat_analysis/human_scrnaseq_0420.harmony.DE.res0.8.csv')


##change resolutions - res0.4
human_harmony <- FindClusters(human_harmony, resolution = 0.4, group.singletons = TRUE)
#Plot UMAP
DimPlot(human_harmony, reduction='umap', split.by='sample', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.harmony.UMAP_res0.4.sample_split.pdf", width = 20)
DimPlot(human_harmony, reduction='umap', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.harmony.UMAP_res0.4.pdf", width = 20)
DimPlot(human_harmony, reduction='umap', group.by = "sample", pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.harmony.UMAP_res0.4.sample.pdf", width = 20)

#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_harmony$seurat_clusters, human_harmony$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/human_scrnaseq_0420.harmony.counts_cluster_sample.res0.4.csv')
#Identify the number of cells in each cluster between genotypes
counts_cluster_time =  table(human_harmony$seurat_clusters, human_harmony$time)
write.csv(counts_cluster_sample, file='./seurat_analysis/human_scrnaseq_0420.harmony.counts_cluster_time.res0.4.csv')

#Find all markers that define each cluster
human_harmony.markers <- FindAllMarkers(human_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_harmony.markers, file='./seurat_analysis/human_scrnaseq_0420.harmony.markers.res0.4.csv')

#Perform differential expression
human_harmony_time = human_harmony
Idents(human_harmony_time) <- 'time'
human_harmony.diffexp <- FindAllMarkers(human_harmony_time, min.pct=0.1, logfc.threshold = 0.25)
write.csv(human_harmony.diffexp, file='./seurat_analysis/human_scrnaseq_0420.harmony.DE.res0.4.csv')

#Save object to avoid needing to re-run previous computations
saveRDS(human_harmony, file = "./seurat_analysis/human_harmony.rds")

