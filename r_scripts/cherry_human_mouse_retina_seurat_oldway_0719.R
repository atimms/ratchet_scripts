library(Seurat)
library(cowplot)
library(ggplot2)
setwd('/data/atimms/cherry_single_cell_0717/mouse_human_retina_seurat_0719')


mouse_data <- read.table("GSE63472_P14Retina_merged_digital_expression.shared_genes.txt",header=T,row.names = 1)
human_data <- read.table("cherry_hu148_0817.raw_martix.shared_genes.txt",header=T,row.names = 1)

# Set up mouse object
mouse <- CreateSeuratObject(counts = mouse_data, project = "mouse_retina", min.cells = 5)
mouse$species <- "mouse"
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^MT-")
mouse <- subset(mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
mouse <- NormalizeData(mouse, verbose = FALSE)
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
# Visualize QC metrics as a violin plot
#VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Set up human object
human <- CreateSeuratObject(counts = human_data, project =  "human_retina", min.cells = 5)
human$species <- "human"
human[["percent.mt"]] <- PercentageFeatureSet(human, pattern = "^MT-")
human <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
human <- NormalizeData(human, verbose = FALSE)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)
# Visualize QC metrics as a violin plot
#VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Perform integration
retina.anchors <- FindIntegrationAnchors(object.list = list(mouse, human), dims = 1:20)
retina.combined <- IntegrateData(anchorset = retina.anchors, dims = 1:20)

#Perform an integrated analysis
DefaultAssay(retina.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
retina.combined <- ScaleData(retina.combined, verbose = FALSE)
retina.combined <- RunPCA(retina.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
retina.combined <- RunUMAP(retina.combined, reduction = "pca", dims = 1:20)
retina.combined <- FindNeighbors(retina.combined, reduction = "pca", dims = 1:20)
##which is best resolution, 0.5 is default
retina.combined <- FindClusters(retina.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(retina.combined, reduction = "umap", group.by = "species")
p1
dev.copy2pdf(file="cherry_human_mouse_seurat_0719.umap_species.pdf", width = 10)
p2 <- DimPlot(retina.combined, reduction = "umap", label = TRUE)
p2
dev.copy2pdf(file="cherry_human_mouse_seurat_0719.umap_cluster.pdf", width = 10)
plot_grid(p1, p2)
DimPlot(retina.combined, reduction = "umap", split.by = "species")
dev.copy2pdf(file="cherry_human_mouse_seurat_0719.umap_cluster_species.pdf", width = 10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
retina.combined.markers <- FindAllMarkers(retina.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(retina.combined.markers, file="cherry_human_mouse_seurat_0719.cluster_markers.csv")

##save file
saveRDS(retina.combined, file = "mouse_human_retina_sc_0719.rds")

##and then load
retina.combined <- readRDS(file ="mouse_human_retina_sc_0719.rds")

##name clusters
new.cluster.ids <- c("rods", "bipolars (ON-type)", "cones", "bipolars (ON-type)", "amacrines (GABAergic)", "amacrines (glycinergic)", "bipolars (OFF-type)", "müller glial cells", "rods", "amacrines (GABAerigc)", "ganglions", "bipolars (ON-type)", "horizontals", "vasculature (pericytes, endothelium?)", "amacrines (glycinergic)", "amacrines (ChAT, NPY, GABAergic)", "cones", "vasculature (pericytes, endothelium?)", "microglia", "astrocytes")
names(new.cluster.ids) <- levels(retina.combined)
retina.combined <- RenameIdents(retina.combined, new.cluster.ids)
DimPlot(retina.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.copy2pdf(file="cherry_human_mouse_seurat_0719.umap_named.pdf", width = 10)

rods <- subset(retina.combined, idents = "rods")
Idents(rods) <- "species"
avg.rods <- log1p(AverageExpression(rods, verbose = FALSE)$RNA)
avg.rods$gene <- rownames(avg.rods)
p1 <- ggplot(avg.rods, aes(mouse, human)) + geom_point() + ggtitle("rods")
p1
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

retina.combined$celltype.species <- paste(Idents(retina.combined), retina.combined$species, sep = "_")
retina.combined$celltype <- Idents(retina.combined)
Idents(retina.combined) <- "celltype.species"
b.rods.response <- FindMarkers(retina.combined, ident.1 = "rods_mouse", ident.2 = "rods_human", verbose = FALSE)
head(b.rods.response, n = 15)
b.cones.response <- FindMarkers(retina.combined, ident.1 = "cones_mouse", ident.2 = "cones_human", verbose = FALSE)
head(b.cones.response, n = 15)

Idents(retina.combined) <- factor(Idents(retina.combined), levels = c("rods", "bipolars (ON-type)", "cones", "amacrines (GABAergic)", "amacrines (glycinergic)", "bipolars (OFF-type)", "müller glial cells", "ganglions", "horizontals", "amacrines (ChAT, NPY, GABAergic)", "microglia", "astrocytes"))
markers.to.plot <- c("PPP1R42", "ARL4A", "GRAMD2", "RAMP3")
DotPlot(retina.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "species") + RotatedAxis()

DotPlot(retina.combined, features = markers.to.plot)
