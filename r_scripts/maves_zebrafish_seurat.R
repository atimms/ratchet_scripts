library(Seurat)
library(ggplot2)

workingDir = "/data/atimms/maves_zf_single_cell_0619";
setwd(workingDir);


#Import Cluster names
Cluster_names <- read.csv(file = "GSE112294_ClusterNames.csv")
#write.csv(Cluster_names, file = "/archive/dobyns_w/sc_data_1118/GSE112294_ClusterNames.csv")

#Read in hpf04 data
hpf04_clustID <- read.table(file = "GSM3067189_04hpf_clustID.txt")
head(hpf04_clustID)
hpf04 <- read.csv(file = "GSM3067189_04hpf.csv", row.names = 1)
hpf04 <- as.data.frame(hpf04)
hpf04[1:5,1:5]

#Seurat Pipeline
hpf04_seurat <- CreateSeuratObject(counts = hpf04, project = "hpf04_zebrafish", meta.data = hpf04_clustID)
hpf04_seurat <- NormalizeData(hpf04_seurat)
hpf04_seurat <- FindVariableFeatures(hpf04_seurat)
hpf04_seurat <- ScaleData(hpf04_seurat)
hpf04_seurat <- RunPCA(hpf04_seurat, npcs = 50)
ElbowPlot(hpf04_seurat, ndims = 50)
hpf04_seurat <- FindNeighbors(hpf04_seurat, dims = 1:20)
hpf04_seurat <- FindClusters(hpf04_seurat, resolution = 1)
hpf04_seurat <- RunTSNE(hpf04_seurat, dims = 1:20)
DimPlot(hpf04_seurat, reduction = "tsne", label = T) + NoLegend()

#Add clusterID
hpf04_seurat$clustID <- hpf04_clustID$V1
DimPlot(hpf04_seurat, reduction = "tsne", label = F, group.by = "clustID")
table(hpf04_seurat@meta.data$clustID)

#Save hpf04 data set
saveRDS(hpf04_seurat, file = "hpf04_seurat.rds")
readRDS(file = "hpf04_seurat.rds")

#Load in hpf24 data set
hpf24_clustID <- read.table(file = "/data/zthoms/GSM3067195_24hpf_clustID.txt")
head(hpf24_clustID)
hpf24 <- read.csv(file = "/data/zthoms/GSM3067195_24hpf.csv", row.names = 1)
hpf24 <- as.data.frame(hpf24)
hpf24[1:5,1:5]

#Seurat Pipeline
hpf24_seurat <- CreateSeuratObject(counts = hpf24, project = "hpf24_zebrafish", meta.data = hpf24_clustID)
hpf24_seurat <- NormalizeData(hpf24_seurat)
hpf24_seurat <- FindVariableFeatures(hpf24_seurat)
hpf24_seurat <- ScaleData(hpf24_seurat)
hpf24_seurat <- RunPCA(hpf24_seurat, npcs = 50)
ElbowPlot(hpf24_seurat, ndims = 50)
hpf24_seurat <- FindNeighbors(hpf24_seurat, dims = 1:40)
hpf24_seurat <- FindClusters(hpf24_seurat, resolution = 2)
hpf24_seurat <- RunTSNE(hpf24_seurat, dims = 1:40)
DimPlot(hpf24_seurat, reduction = "tsne", label = T) + NoLegend()

#Add ClusterID
hpf24_seurat$clustID <- hpf24_clustID$V1
DimPlot(hpf24_seurat, reduction = "tsne", label = F, group.by = "clustID")
VlnPlot(hpf24_seurat, features = "nFeature_RNA", group.by = "orig.ident", log = F)
#hpf24_seurat <- SetIdent(hpf24_seurat, value = "clustID")
table(hpf24_seurat@meta.data$clustID)

#SaveRDS file
saveRDS(hpf24_seurat, file = "/data/zthoms/hpf24_seurat.rds")
hpf24_seurat <- readRDS(file = "/data/zthoms/hpf24_seurat.rds")

FeaturePlot(hpf24_seurat, features = c("aspm", "gfap"))


###Below this was investigating for a different PI###
hpf24_subset <- SubsetData(hpf24_seurat, ident.use = c(130, 137, 142, 147, 157,165,168,170))
hpf24_subset <- NormalizeData(hpf24_subset)
hpf24_subset <- FindVariableFeatures(hpf24_subset)
hpf24_subset <- ScaleData(hpf24_subset)
hpf24_subset <- RunPCA(hpf24_subset, npcs = 50)
#hpf24_subset <- RunPCA(hpf24_subset, npcs = 50, features = VariableFeatures(hpf24_subset), verbose = F)
ElbowPlot(hpf24_subset, ndims = 50)
hpf24_subset <- FindNeighbors(hpf24_subset, dims = 1:30)
hpf24_subset <- FindClusters(hpf24_subset, resolution = 2)
hpf24_subset <- RunTSNE(hpf24_subset, dims = 1:30)
DimPlot(hpf24_subset, reduction = "tsne", label = F)
DimPlot(hpf24_subset, reduction = "tsne", label = T, group.by = "RNA_snn_res.2")
hpf24_subset <- SetIdent(hpf24_subset, value = "clustID")
de_genes_24_subset <- FindAllMarkers(hpf24_subset, only.pos = T)
write.csv(de_genes_24_subset, file = "/archive/dobyns_w/sc_data_1118/zebrafish_de_genes.csv")
cluster_names <- c("neural - ventral hindbrain", "neural - dorsal spinal cord", "tailbud - spinal cord",
  "neural - dorsal hindbrain", "differentiating neurons - eomesa", "neural - hindbrain roofplate",
  "neural - telencephalon", "hindbrain - gsx1")
names(x = cluster_names) <- levels(x = hpf24_subset)
hpf24_subset <- RenameIdents(hpf24_subset, cluster_names)
avg_exp_24_subset <- AverageExpression(hpf24_subset)
write.csv(avg_exp_24_subset, file = "/archive/dobyns_w/sc_data_1118/zebrafish_avg_exp.csv")
FeaturePlot(hpf24_subset, features = c("aspm", "gfap"))
saveRDS(hpf24_subset, file = "/data/zthoms/hpf24_subset.rds")
hpf24_subset <- readRDS(file = "/data/zthoms/hpf24_subset.rds")

tele_subset <- SubsetData(hpf24_subset, ident.use = "neural - telencephalon")
tele_subset <- NormalizeData(tele_subset)
tele_subset <- FindVariableFeatures(tele_subset)
tele_subset <- ScaleData(tele_subset)
tele_subset <- RunPCA(tele_subset, npcs = 50)
#tele_subset <- RunPCA(tele_subset, npcs = 50, features = VariableFeatures(tele_subset), verbose = F)
ElbowPlot(tele_subset, ndims = 50)
tele_subset <- FindNeighbors(tele_subset, dims = 1:30)
tele_subset <- FindClusters(tele_subset, resolution = 1)
tele_subset <- RunTSNE(tele_subset, dims = 1:30)
DimPlot(tele_subset, reduction = "tsne", label = F)
FeaturePlot(tele_subset, features = c("aspm", "elavl3"))
de_genes_tele <- FindAllMarkers(tele_subset, only.pos = T)
write.csv(de_genes_tele, file = "/archive/dobyns_w/sc_data_1118/de_genes_tele.csv")
