#need to get latest version of seurat i.e.
#install.packages('Seurat')
library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
setwd('/data/atimms/cherry_single_cell_0717/mouse_human_retina_seurat_0719')

##set for integrating the datasets
options(future.globals.maxSize = 4000 * 1024^2)

##need this to use umap
#reticulate::py_install(packages ='umap-learn')

##to get 2x matrices that shared genes
##use python script mouse_human_retina_get_shared_genes_0719.py

##load mouse data
mouse_data <- read.table("GSE63472_P14Retina_merged_digital_expression.shared_genes.txt",header=T,row.names = 1)
mouse_retina <- CreateSeuratObject(counts = mouse_data, project = "mouse_retina_GSE63472", min.cells = 3, min.features = 200)
mouse_retina
mouse_retina$species <- "mouse"
# store mitochondrial percentage in object meta data
mouse_retina <- PercentageFeatureSet(mouse_retina, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
VlnPlot(mouse_retina, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# run sctransform
mouse_retina <- SCTransform(mouse_retina, vars.to.regress = "percent.mt", verbose = FALSE)
warnings()

##load human data
human_data <- read.table("cherry_hu148_0817.raw_martix.shared_genes.txt",header=T,row.names = 1)
human_retina <- CreateSeuratObject(counts = human_data, project = "human_retina_GSE63472", min.cells = 3, min.features = 200)
human_retina
human_retina$species <- "human"
# store mitochondrial percentage in object meta data
human_retina <- PercentageFeatureSet(human_retina, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
VlnPlot(human_retina, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# run sctransform
human_retina <- SCTransform(human_retina, vars.to.regress = "percent.mt", verbose = FALSE)

#select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
retina.list <- list(mouse_retina, human_retina)
retina.features <- SelectIntegrationFeatures(object.list = retina.list, nfeatures = 3000)
retina.list <- PrepSCTIntegration(object.list = retina.list, anchor.features = retina.features, 
                                    verbose = FALSE)
##identify anchors and integrate the datasets.
retina.anchors <- FindIntegrationAnchors(object.list = retina.list, normalization.method = "SCT", 
                                           anchor.features = retina.features, verbose = FALSE)
retina.integrated <- IntegrateData(anchorset = retina.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
##Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
retina.integrated <- RunPCA(retina.integrated, verbose = FALSE)
retina.integrated <- RunUMAP(retina.integrated, dims = 1:30)
retina.integrated <- FindNeighbors(retina.integrated, reduction = "pca", dims = 1:30)
##which is best resolution, 0.5 is default
retina.integrated <- FindClusters(retina.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(retina.integrated, reduction = "umap", group.by = "species")
p1
dev.copy2pdf(file="cherry_human_mouse_seurat_0919.umap_species.pdf", width = 10)
p2 <- DimPlot(retina.integrated, reduction = "umap", label = TRUE)
p2
dev.copy2pdf(file="cherry_human_mouse_seurat_0919.umap_cluster.pdf", width = 10)
plot_grid(p1, p2)
DimPlot(retina.integrated, reduction = "umap", split.by = "species")
dev.copy2pdf(file="cherry_human_mouse_seurat_0919.umap_cluster_species.pdf", width = 10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
retina.integrated.markers <- FindAllMarkers(retina.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(retina.integrated.markers, file="cherry_human_mouse_seurat_0919.cluster_markers.csv")


##save file
saveRDS(retina.integrated, file = "mouse_human_retina_sc_0919.rds")

##and then load
retina.integrated <- readRDS(file ="mouse_human_retina_sc_0919.rds")

##featueplot for 2 genes different in mouse and human
#FeaturePlot(retina.integrated, features = c("VIM", "PDE6H", "APOE", "RIMS1", "PCP2", "SPARC"), min.cutoff = "q9", split.by = "species")
FeaturePlot(retina.integrated, features = c("HTRA1","APOE", "CABP5"), min.cutoff = "q9", split.by = "species")
dev.copy2pdf(file="cherry_human_mouse_seurat_0919.HTRA1_CABP5.FP_species.pdf", width = 10)



##name clusters
new.cluster.ids <- c('rods1', 'rods2', 'amacrine', 'bipolar on-type', 'amacrine subtype?', 'bipolar rod-type', 'müller glia', 'cones', 'ganglions', 'vascular endothelium', 'horizontals', 'amacrine chat-type')
names(new.cluster.ids) <- levels(retina.integrated)
retina.integrated <- RenameIdents(retina.integrated, new.cluster.ids)
DimPlot(retina.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.copy2pdf(file="cherry_human_mouse_seurat_0919.umap_named.pdf", width = 10)


#DotPlot of 
Idents(retina.integrated) <- factor(Idents(retina.integrated), levels = c('rods1', 'rods2', 'amacrine', 'bipolar on-type', 'amacrine subtype?', 'bipolar rod-type', 'müller glia', 'cones', 'ganglions', 'vascular endothelium', 'horizontals', 'amacrine chat-type'))
markers.to.plot <- c("HTRA1", "APOE", "CABP5")
DotPlot(retina.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "species") + RotatedAxis()
dev.copy2pdf(file="cherry_human_mouse_seurat_0919.HTRA1_CABP5.dotplot.pdf")

