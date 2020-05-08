library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(stringr)
setwd('/data/atimms/cherry_10x_organoid_0320')

# Load the 10x datasets for 0319
Mut26.data <- Read10X(data.dir = "/data/atimms/cherry_10x_organoid_0320/Mut2-6_filtered_feature_bc_matrix")
Mut26 <- CreateSeuratObject(counts = Mut26.data, project = "Mut26", min.cells = 3, min.features = 200)
Mut38.data <- Read10X(data.dir = "/data/atimms/cherry_10x_organoid_0320/Mut3-8_filtered_feature_bc_matrix")
Mut38 <- CreateSeuratObject(counts = Mut38.data, project = "Mut38", min.cells = 3, min.features = 200)
Wt21.data <- Read10X(data.dir = "/data/atimms/cherry_10x_organoid_0320/Wt2-1_filtered_feature_bc_matrix")
Wt21 <- CreateSeuratObject(counts = Wt21.data, project = "Wt21", min.cells = 3, min.features = 200)
Wt31.data <- Read10X(data.dir = "/data/atimms/cherry_10x_organoid_0320/Wt3-1_filtered_feature_bc_matrix")
Wt31 <- CreateSeuratObject(counts = Wt31.data, project = "Wt31", min.cells = 3, min.features = 200)

##simple merge of data (should use sctranform to so this properly)
organoid_merge <- merge(Mut26, y = c(Mut38, Wt21, Wt31), add.cell.ids = c("Mut26", "Mut38", "Wt21", "Wt31"), project = "organoid_merge")
head(colnames(organoid_merge))
unique(sapply(X = strsplit(colnames(organoid_merge), split = "_"), FUN = "[", 1))
##number of cells
table(organoid_merge$orig.ident)

# store mitochondrial percentage in object meta data (mt- is for mouse)
organoid_merge[["percent.mt"]] <- PercentageFeatureSet(organoid_merge, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(organoid_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.copy2pdf(file="organoid_merge_0320.qc.pdf", width = 10)

# filter the data
organoid_merge <- subset(organoid_merge, subset = nFeature_RNA > 400 & nFeature_RNA < 2000 & percent.mt < 10)

# run sctransform
organoid_merge <- SCTransform(organoid_merge, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
organoid_merge <- RunPCA(organoid_merge, verbose = FALSE)
organoid_merge <- RunUMAP(organoid_merge, dims = 1:30, verbose = FALSE)
organoid_merge <- FindNeighbors(organoid_merge, dims = 1:30, verbose = FALSE)
organoid_merge <- FindClusters(organoid_merge, verbose = FALSE)
DimPlot(organoid_merge, label = TRUE) 
dev.copy2pdf(file="organoid_merge_0320.umap_cluster.pdf", width = 20)
p1 <- DimPlot(organoid_merge, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(organoid_merge, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.copy2pdf(file="organoid_merge_0320.umap_cluster_split.pdf", width = 20)

##number of cells
table(organoid_merge$orig.ident)

###not setup yet..
# find markers for every cluster compared to all remaining cells, report only the positive ones
organoid_5wk_0320.integrated.markers <- FindAllMarkers(organoid_5wk_0320, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(organoid_5wk_0320.integrated.markers, file="organoid_5wk_0320.cluster_markers.res0.4.csv")




