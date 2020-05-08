library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(stringr)
setwd('/data/atimms/cherry_single_cell_0717/mouse_human_retina_seurat_0719')


##load mouse data
mouse_data <- read.table("GSE63472_P14Retina_merged_digital_expression.txt",header=T,row.names = 1)
mouse_retina <- CreateSeuratObject(counts = mouse_data, project = "mouse_retina_GSE63472", min.cells = 3, min.features = 200)
mouse_retina
# store mitochondrial percentage in object meta data (mt- is for mouse)
mouse_retina[["percent.mt"]] <- PercentageFeatureSet(mouse_retina, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(mouse_retina, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.copy2pdf(file="mouse_retina_1219.qc.pdf", width = 10)

# filter the data
mouse_retina <- subset(mouse_retina, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# run sctransform
mouse_retina <- SCTransform(mouse_retina, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
mouse_retina <- RunPCA(mouse_retina, verbose = FALSE)
mouse_retina <- RunUMAP(mouse_retina, dims = 1:30, verbose = FALSE)
mouse_retina <- FindNeighbors(mouse_retina, dims = 1:30, verbose = FALSE)
mouse_retina <- FindClusters(mouse_retina, verbose = FALSE, resolution = 0.4)
DimPlot(mouse_retina, label = TRUE) 
dev.copy2pdf(file="mouse_retina_1219.umap_cluster.res0.4.pdf", width = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
mouse_retina.integrated.markers <- FindAllMarkers(mouse_retina, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mouse_retina.integrated.markers, file="mouse_retina_1219.cluster_markers.res0.4.csv")

mouse_retina <- FindClusters(mouse_retina, verbose = FALSE, resolution = 0.3)
DimPlot(mouse_retina, label = TRUE) 
dev.copy2pdf(file="mouse_retina_1219.umap_cluster.res0.3.pdf", width = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
mouse_retina.integrated.markers <- FindAllMarkers(mouse_retina, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mouse_retina.integrated.markers, file="mouse_retina_1219.cluster_markers.res0.3.csv")

mouse_retina <- FindClusters(mouse_retina, verbose = FALSE, resolution = 0.2)
DimPlot(mouse_retina, label = TRUE) 
dev.copy2pdf(file="mouse_retina_1219.umap_cluster.res0.2.pdf", width = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
mouse_retina.integrated.markers <- FindAllMarkers(mouse_retina, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mouse_retina.integrated.markers, file="mouse_retina_1219.cluster_markers.res0.2.csv")

# How many cells are in each cluster
mouse_cell_counts = table(Idents(mouse_retina))
write.csv(mouse_cell_counts, file="mouse_retina_1219.cell_counts.res0.2.csv")

# rename clusters
mouse_retina_renamed <- RenameIdents(object = mouse_retina, '0' = 'Rod photoreceptors', '1' = 'Rod photoreceptors', '2' = 'Bipolar cells',
                            '3' = 'Cone photoreceptors', '4' = 'Amacrine cells', '5' = 'Bipolar cells', '6' = 'Amacrine cells',
                            '7' = 'Müller glia', '8' = 'Ganglion cells', '9' = '---', '10' = '--', '11' = 'Microglia',
                            '12' = 'Rod photoreceptors', '13' = 'Rod photoreceptors', '14' = 'Amacrine cells',
                            '15' = 'Horizontal cells', '16' = 'Microglia', '17' = 'Müller glia/astrocytes')
DimPlot(mouse_retina_renamed, label = TRUE) 
dev.copy2pdf(file="mouse_retina_1219_renamed.umap_cluster.res0.2.pdf", width = 10, height = 10)

# violin plots 
VlnPlot(mouse_retina_renamed, features = c("CABP5", "HTRA1"))
dev.copy2pdf(file="mouse_retina_1219_renamed.violin_plot.CABP5_HTRA1.pdf", width = 10)

# feature counts
FeaturePlot(mouse_retina_renamed, features = c("CABP5"))
dev.copy2pdf(file="mouse_retina_1219_renamed.feature_plot.CABP5.pdf", width = 10, height = 10)
FeaturePlot(mouse_retina_renamed, features = c("HTRA1"))
dev.copy2pdf(file="mouse_retina_1219_renamed.feature_plot.HTRA1.pdf", width = 10, height = 10)



##load human data
human_data <- read.table("cherry_hu148_0817.raw_martix.txt",header=T,row.names = 1)
human_retina <- CreateSeuratObject(counts = human_data, project = "human_148", min.cells = 3, min.features = 200)
human_retina
# store mitochondrial percentage in object meta data (mt- is for human)
human_retina <- PercentageFeatureSet(human_retina, pattern = "^mt-", col.name = "percent.mt")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
human_retina[["percent.mt"]] <- PercentageFeatureSet(human_retina, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(human_retina, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.copy2pdf(file="human_retina_1219.qc.pdf", width = 10)

# filter the data
human_retina <- subset(human_retina, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)

# run sctransform
human_retina <- SCTransform(human_retina, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
human_retina <- RunPCA(human_retina, verbose = FALSE)
human_retina <- RunUMAP(human_retina, dims = 1:30, verbose = FALSE)
human_retina <- FindNeighbors(human_retina, dims = 1:30, verbose = FALSE)
human_retina <- FindClusters(human_retina, verbose = FALSE)
DimPlot(human_retina, label = TRUE) 
dev.copy2pdf(file="human_retina_1219.umap_cluster.pdf", width = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
human_retina.integrated.markers <- FindAllMarkers(human_retina, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_retina.integrated.markers, file="human_retina_1219.cluster_markers.csv")

# How many cells are in each cluster
human_cell_counts = table(Idents(human_retina))
write.csv(human_cell_counts, file="human_retina_1219.cell_counts.csv")

# rename clusters
human_retina_renamed <- RenameIdents(object = human_retina, '0' = 'Rod photoreceptors', '1' = 'Rod photoreceptors', '2' = 'Müller glia',
                                     '3' = 'Bipolar cells', '4' = 'Amacrine cells', '5' = 'Bipolar cells - on', '6' = 'Rod photoreceptors',
                                     '7' = 'Horizontal cells', '8' = 'Mitochondrial genes', '9' = 'Cone photoreceptors', '10' = 'Ganglion cells', '11' = 'Rod photoreceptors',
                                     '12' = 'Müller glia/astrocytes', '13' = 'Amacrine cells', '14' = 'Rod photoreceptors',
                                     '15' = 'Amacrine cells', '16' = 'Microglia', '17' = 'Pericytes', '18' = '--', '19' = 'Rod photoreceptors')
DimPlot(human_retina_renamed, label = TRUE) 
dev.copy2pdf(file="human_retina_1219_renamed.umap_cluster.pdf", width = 10, height = 10)

DimPlot(human_retina_renamed, group.by = "orig.ident") 
dev.copy2pdf(file="human_retina_1219_renamed.umap_cluster.sample_id.pdf", width = 10, height = 10)



# violin plots 
VlnPlot(human_retina_renamed, features = c("CABP5", "HTRA1"))
dev.copy2pdf(file="human_retina_1219_renamed.violin_plot.CABP5_HTRA1.pdf", width = 10)

# feature counts
FeaturePlot(human_retina_renamed, features = c("CABP5"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.CABP5.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("HTRA1"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.HTRA1.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("RP1"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.RP1.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("PDE6H"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.PDE6H.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("OPN1LW"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.OPN1LW.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("OPN1SW"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.OPN1SW.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("NR2E3"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.NR2E3.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("ARR3"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.ARR3.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("ONECUT2"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.ONECUT2.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("VSX2"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.VSX2.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("CLU"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.CLU.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("GAD1"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.GAD1.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("NEFM"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.NEFM.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("PBX1", "PBX2", "PBX3", "PBX4"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.PBX_genes.pdf", width = 10, height = 10)
FeaturePlot(human_retina_renamed, features = c("PKNOX1", "PKNOX2"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.PKNOX_genes.pdf", width = 10, height = 5)
FeaturePlot(human_retina_renamed, features = c("ASCL1", "NFIA", "NFIB", "NFIX", "AQP4"))
dev.copy2pdf(file="human_retina_1219_renamed.feature_plot.other_genes.pdf", width = 10, height = 15)

