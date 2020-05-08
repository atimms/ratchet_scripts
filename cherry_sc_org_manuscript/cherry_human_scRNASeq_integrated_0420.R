# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA')

##if data isn't loaded
human_precheck <- readRDS(file = "./seurat_analysis/human_precheck.rds")

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)


#First, split the merged Seurat object by sample
human.list <- SplitObject(human_precheck, split.by = 'sample')

#Run SCTransform on each dataset with "var.to.regress=percent.mt"
for (i in 1:length(human.list)) {
  human.list[[i]] <- SCTransform(human.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}

#Select features for downstream integration
human.features <- SelectIntegrationFeatures(object.list = human.list, nfeatures = 3000)
human.list <- PrepSCTIntegration(object.list = human.list, anchor.features = human.features, 
                                 verbose = FALSE)

#Find anchors and integrate data
human.anchors <- FindIntegrationAnchors(object.list = human.list, normalization.method = "SCT", 
                                        anchor.features = human.features, verbose = FALSE)
human_integrated <- IntegrateData(anchorset = human.anchors, normalization.method = "SCT", 
                                  verbose = FALSE)


#Run standard Seurat workflow (PCA, UMAP, FindClusters, FindNeighbors) but with assay set to 'integrated'
human_integrated = RunPCA(human_integrated, assay='integrated', verbose=FALSE)
human_integrated <- RunUMAP(human_integrated, assay='integrated', dims = 1:30)
human_integrated <- FindNeighbors(human_integrated, assay='integrated', dims = 1:30, verbose=FALSE)
human_integrated <- FindClusters(human_integrated, assay = 'integrated', resolution = 0.5)

#Plot UMAP
DimPlot(human_integrated, reduction='umap', split.by='sample', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.integrated.UMAP_res0.5.sample_split.pdf", width = 20)
DimPlot(human_integrated, reduction='umap', pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.integrated.UMAP_res0.5.pdf", width=20)
DimPlot(human_integrated, reduction='umap', group.by = "sample", pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.integrated.UMAP_res0.5.sample.pdf", width=20)


#Save object to avoid needing to re-run previous computations
saveRDS(human_integrated, file = "./seurat_analysis/human_integrated.rds")

#Identify the number of cells in each cluster between samples/time
counts_cluster_sample = table(human_integrated$seurat_clusters, human_integrated$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/human_scrnaseq_0420.integrated.counts_cluster_sample.res0.5.csv')
counts_cluster_time = table(human_integrated$seurat_clusters, human_integrated$time)
write.csv(counts_cluster_time, file='./seurat_analysis/human_scrnaseq_0420.integrated.counts_cluster_time.res0.5.csv')

#Change Deafualt Assay to SCT or RNA for downstream marker and expression analyses. Let's start with SCT.
DefaultAssay(human_integrated) <- 'SCT'

#Find all markers that define each cluster
human_integrated_SCT.markers <- FindAllMarkers(human_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_integrated_SCT.markers, file='./seurat_analysis/human_scrnaseq_0420.integrated.markers.res0.5.csv')

#Perform differential expression
human_integrated_pheno = human_integrated
Idents(human_integrated_pheno) <- 'time'
human_integrated_SCT.diffexp <- FindAllMarkers(human_integrated_pheno, min.pct=0.1, logfc.threshold = 0.25)
write.csv(human_integrated_SCT.diffexp, file='./seurat_analysis/human_scrnaseq_0420.merged.DE.res0.5.csv')










#Re-do expression/marker analyses in RNA assay ---- not used
DefaultAssay(human_integrated) <- 'RNA'
#Normalize RNA data for visualization purposes
human_integrated <- NormalizeData(human_integrated, verbose = FALSE)
#Check that normalization worked by looking at expression of genes of interest
FeaturePlot(human_integrated, features='LINC00461')
DotPlot(human_integrated, features='LINC00461', group.by='pheno')

#Find all markers that define each cluster
human_integrated_RNA.markers <- FindAllMarkers(human_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_integrated_RNA.markers, file='human_markers_integrated_RNA.csv')

#Run DE
human_integrated_pheno = human_integrated
Idents(human_integrated_pheno) <- 'pheno'
human_integrated_RNA.diffexp <- FindAllMarkers(human_integrated_pheno, min.pct=0.1, logfc.threshold = 0.25)
write.csv(human_integrated_RNA.diffexp, file='human_DE_integrated_RNA.csv')


