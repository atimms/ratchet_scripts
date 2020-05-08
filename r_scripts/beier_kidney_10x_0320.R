library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(stringr)
setwd('/data/atimms/beier_kidney_10x_0320')

###two new samples pointless not used

# Load the 10x datasets for 0319
Kidney_mouse.data <- Read10X(data.dir = "/data/atimms/beier_kidney_10x_0320/Kidney_mouse_filtered_feature_bc_matrix")
Kidney_mouse <- CreateSeuratObject(counts = Kidney_mouse.data, project = "Kidney_mouse", min.cells = 3, min.features = 200)
Kidney_cystic.data <- Read10X(data.dir = "/data/atimms/beier_kidney_10x_0320/Kidney_cystic_filtered_feature_bc_matrix")
Kidney_cystic <- CreateSeuratObject(counts = Kidney_cystic.data, project = "Kidney_cystic", min.cells = 3, min.features = 200)

##simple merge of data (should use sctranform to so this properly)
kidney_merge <- merge(Kidney_mouse, y = c(Kidney_cystic), add.cell.ids = c("Kidney_mouse", "Kidney_cystic"), project = "kidney_merge")
head(colnames(kidney_merge))
##number of cells
table(kidney_merge$orig.ident)


# store mitochondrial percentage in object meta data (mt- is for mouse)
kidney_merge[["percent.mt"]] <- PercentageFeatureSet(kidney_merge, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(kidney_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.copy2pdf(file="beier_kidney_0320.qc.pdf", width = 10)

# filter the data
kidney_merge <- subset(kidney_merge, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)
table(kidney_merge$orig.ident)

# run sctransform
kidney_merge <- SCTransform(kidney_merge, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
kidney_merge <- RunPCA(kidney_merge, verbose = FALSE)
kidney_merge <- RunUMAP(kidney_merge, dims = 1:30, verbose = FALSE)
kidney_merge <- FindNeighbors(kidney_merge, dims = 1:30, verbose = FALSE)
kidney_merge <- FindClusters(kidney_merge, verbose = FALSE)
DimPlot(kidney_merge, label = TRUE) 
dev.copy2pdf(file="beier_kidney_0320.umap_cluster.pdf", width = 20)
p1 <- DimPlot(kidney_merge, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(kidney_merge, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.copy2pdf(file="beier_kidney_0320.umap_cluster_split.pdf", width = 20)
DimPlot(kidney_merge, reduction = "umap", split.by = "orig.ident")
dev.copy2pdf(file="beier_kidney_0320.umap_cluster_split2.pdf", width = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
kidney_merge.integrated.markers <- FindAllMarkers(kidney_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(kidney_merge.integrated.markers, file="beier_kidney_0320.cluster_markers.csv")
##change default assay to RNA abd then normalize, do I need to scale too?
DefaultAssay(kidney_merge) <- "RNA"
kidney_merge <- NormalizeData(kidney_merge)
kidney_merge.integrated.markers <- FindAllMarkers(kidney_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(kidney_merge.integrated.markers, file="beier_kidney_0320.cluster_markers_rna.csv")

##counts per cluster??
count.by.cluster = table(Idents(kidney_merge), kidney_merge$orig.ident)
write.csv(count.by.cluster, file="beier_kidney_0320.cell_count_by_cluster.csv")


##markers by condition
##not used
#kidney_merge$celltype.id <- paste(Idents(kidney_merge), kidney_merge$orig.ident, sep = "_")
#kidney_merge$celltype <- Idents(kidney_merge)
##use orig.ident rather than cluster
Idents(kidney_merge) <- "orig.ident"
marker.by.condition <- FindAllMarkers(kidney_merge, assay=RNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(marker.by.condition, n = 15)
write.csv(marker.by.condition, file="beier_kidney_0320.sample_markers.csv")
##change default assay to RNA
DefaultAssay(kidney_merge) <- "RNA"
marker.by.condition <- FindAllMarkers(kidney_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(marker.by.condition, n = 15)
write.csv(marker.by.condition, file="beier_kidney_0320.sample_markers_rna.csv")

