library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(stringr)
setwd('/data/atimms/beier_kidney_10x_0320')

# Load the 10x datasets
Kidney_mouse.data <- Read10X(data.dir = "/data/atimms/beier_kidney_10x_0320/Kidney_mouse_filtered_feature_bc_matrix")
Kidney_mouse <- CreateSeuratObject(counts = Kidney_mouse.data, project = "Kidney_mouse", min.cells = 3, min.features = 200)

# store mitochondrial percentage in object meta data (mt- is for mouse)
Kidney_mouse[["percent.mt"]] <- PercentageFeatureSet(Kidney_mouse, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(Kidney_mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.copy2pdf(file="beier_kidney_0520.qc.pdf", width = 10)

# filter the data
Kidney_mouse <- subset(Kidney_mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
table(Kidney_mouse$orig.ident)

# run sctransform
Kidney_mouse <- SCTransform(Kidney_mouse, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
Kidney_mouse <- RunPCA(Kidney_mouse, verbose = FALSE)
Kidney_mouse <- RunUMAP(Kidney_mouse, dims = 1:30, verbose = FALSE)
Kidney_mouse <- FindNeighbors(Kidney_mouse, dims = 1:30, verbose = FALSE)
Kidney_mouse <- FindClusters(Kidney_mouse, verbose = FALSE)
DimPlot(Kidney_mouse, label = TRUE) 
dev.copy2pdf(file="beier_kidney_0520.umap_cluster.pdf", width = 20)


# find markers for every cluster compared to all remaining cells, report only the positive ones
kidney.integrated.markers <- FindAllMarkers(Kidney_mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(kidney.integrated.markers, file="beier_kidney_0520.cluster_markers.csv")

##counts per cluster??
count.by.cluster = table(Idents(Kidney_mouse), Kidney_mouse$orig.ident)
write.csv(count.by.cluster, file="beier_kidney_0520.cell_count_by_cluster.csv")

##make heatmap
##scale data using all genes in RNA assay
DefaultAssay(Kidney_mouse) <- 'RNA'
Kidney_mouse <- ScaleData(object = Kidney_mouse, features = rownames(Kidney_mouse))
##genes to graph
genes <- read.csv('genelist_0520.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
genes <- genes$genes #vectorize
DoHeatmap(Kidney_mouse, assay = "RNA", features = genes, size = 3)
dev.copy2pdf(file="beier_kidney_0520.heatmap.pdf", width = 20)


