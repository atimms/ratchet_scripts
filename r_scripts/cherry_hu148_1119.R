#need to get latest version of seurat i.e.
#install.packages('Seurat')
library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
setwd('/data/atimms/cherry_single_cell_0717/seurat_analysis_1017')

##load human data
human_data <- read.table("cherry_hu148_0817.raw_martix.txt",header=T,row.names = 1)
human_retina <- CreateSeuratObject(counts = human_data, project = "human_retina", min.cells = 3, min.features = 200)
human_retina
# store mitochondrial percentage in object meta data
human_retina <- PercentageFeatureSet(human_retina, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
VlnPlot(human_retina, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# run sctransform
human_retina <- SCTransform(human_retina, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
human_retina <- RunPCA(human_retina, verbose = FALSE)
human_retina <- RunUMAP(human_retina, dims = 1:30, verbose = FALSE)
human_retina <- FindNeighbors(human_retina, dims = 1:30, verbose = FALSE)
human_retina <- FindClusters(human_retina, verbose = FALSE)
DimPlot(human_retina, label = TRUE) 
dev.copy2pdf(file="cherry_human148_seurat_0919.umap_cluster.pdf", width = 20)


# genes as violin plots.
#VlnPlot(human_retina, features = c('MEIS1', 'MEIS2', 'MEIS3', 'IRX1', 'IRX2', 'IRX3', 'IRX4', 'IRX5', 'IRX6', 'MKX', 'PBX1', 'PBX2', 'PBX3', 'PBX4', 'PKNOX1', 'PKNOX2', 'TGIF1', 'TGIF2', 'TGIF2LX', 'TGIF2LY', 'ASCL1', 'NFIA', 'NFIB', 'NFIX', 'APOE', 'CLU', 'RLBP', 'AQP4'), 
#        pt.size = 0.2, ncol = 4)
VlnPlot(human_retina, features = c('MEIS1', 'MEIS2', 'MEIS3'), 
        pt.size = 0.2, ncol = 4)
dev.copy2pdf(file="cherry_human148_seurat_0919.MEIS.violin.pdf")
# genes on the sctransform embedding.
#FeaturePlot(human_retina, features = c('MEIS1', 'MEIS2', 'MEIS3', 'IRX1', 'IRX2', 'IRX3', 'IRX4', 'IRX5', 'IRX6', 'MKX', 'PBX1', 'PBX2', 'PBX3', 'PBX4', 'PKNOX1', 'PKNOX2', 'TGIF1', 'TGIF2', 'TGIF2LX', 'TGIF2LY', 'ASCL1', 'NFIA', 'NFIB', 'NFIX', 'APOE', 'CLU', 'RLBP', 'AQP4'), 
 #           pt.size = 0.2, ncol = 4)
FeaturePlot(human_retina, features = c('MEIS1', 'MEIS2', 'MEIS3'))
dev.copy2pdf(file="cherry_human148_seurat_0919.MEIS.clustered.pdf", width = 10)

FeaturePlot(human_retina, features = c('AQP4', 'APOE'))
dev.copy2pdf(file="cherry_human148_seurat_0919.AQP4_APOE.clustered.pdf", width = 10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
retina.integrated.markers <- FindAllMarkers(human_retina, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(retina.integrated.markers, file="cherry_human148_seurat_0919.cluster_markers.csv")

##save file
saveRDS(human_retina, file = "cherry_human148_seurat_0919.rds")


