# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA')

##if data isn't loaded
human <- readRDS(file = "./seurat_analysis/human.rds")

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)


#Run SCTransform, set var.to.regress to percent.mt
human_merged <- SCTransform(human, vars.to.regress = "percent.mt", verbose = FALSE)

#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
human_merged <- RunPCA(human_merged, verbose = FALSE)
human_merged <- RunUMAP(human_merged, dims = 1:30, verbose = FALSE)
human_merged <- FindNeighbors(human_merged, dims = 1:30, verbose=FALSE)
human_merged <- FindClusters(human_merged, resolution = 0.5)

#Plot UMAP
DimPlot(human_merged, reduction='umap', split.by='sample', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.merged.UMAP_res0.5.sample_split.pdf", width = 20)
DimPlot(human_merged, reduction='umap', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.merged.UMAP_res0.5.pdf", width=20)
DimPlot(human_merged, reduction='umap', group.by = "sample", pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.merged.UMAP_res0.5.sample.pdf", width=20)

#Save object
saveRDS(human_merged, file = "./seurat_analysis/human_merged.rds")

#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_merged$seurat_clusters, human_merged$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/human_scrnaseq_0420.merged.counts_cluster_sample.res0.5.csv')
#Identify the number of cells in each cluster between genotypes
counts_cluster_time = table(human_merged$seurat_clusters, human_merged$time)
write.csv(counts_cluster_time, file='./seurat_analysis/human_scrnaseq_0420.merged.counts_cluster_time.res0.5.csv')

#Find all markers that define each cluster
human_merged.markers <- FindAllMarkers(human_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged.markers, file='./seurat_analysis/human_scrnaseq_0420.merged.markers.res0.5.csv')

#In order to perform differential expression, use Idents to change the metadata column to phenotype
human_merged_pheno = human_merged
Idents(human_merged_pheno) <- 'time'
human_merged.diffexp <- FindAllMarkers(human_merged_pheno, min.pct=0.1, logfc.threshold = 0.25)
write.csv(human_merged.diffexp, file='./seurat_analysis/human_scrnaseq_0420.merged.DE.res0.5.csv')




#Default assay after SCTransform is SCT. Change to RNA and re-run FindMarkers --- not used
DefaultAssay(human_merged) <- 'RNA'
human_merged_RNA.markers <- FindAllMarkers(human_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_merged_RNA.markers, file='human_markers_merged_RNA.csv')

#Run differential expression
human_merged_phenoRNA = human_merged
Idents(human_merged_phenoRNA) <- 'pheno'
human_merged_RNA.diffexp <- FindAllMarkers(human_merged_phenoRNA, min.pct=0.1, logfc.threshold = 0.25)
write.csv(human_merged_RNA.diffexp, file='human_DE_merged_RNA.csv')
