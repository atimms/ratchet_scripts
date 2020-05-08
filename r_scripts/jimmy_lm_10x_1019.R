library(dplyr)
library(Seurat)
setwd('/data/atimms/jimmy_lm_10x_1019')

##need this to use umap
#reticulate::py_install(packages ='umap-learn')

# Load the 10x dataset
lm1.data <- Read10X(data.dir = "/data/atimms/jimmy_lm_10x_1019/10x_data")
# Initialize the Seurat object with the raw (non-normalized data).
lm_obj <- CreateSeuratObject(counts = lm1.data, project = "lm1", min.cells = 3, min.features = 200)
lm_obj

# store mitochondrial percentage in object meta data
lm_obj <- PercentageFeatureSet(lm_obj, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
lm_obj <- SCTransform(lm_obj, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
lm_obj <- RunPCA(lm_obj, verbose = FALSE)
lm_obj <- RunUMAP(lm_obj, dims = 1:30, verbose = FALSE)

lm_obj <- FindNeighbors(lm_obj, dims = 1:30, verbose = FALSE)
lm_obj <- FindClusters(lm_obj, verbose = FALSE)
DimPlot(lm_obj, label = TRUE) + NoLegend()
dev.copy2pdf(file="jimmy_lm1_10x_1019.umap.pdf", width = 10)


# find markers for every cluster compared to all remaining cells, report only the positive ones
lm_obj.markers <- FindAllMarkers(lm_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(lm_obj.markers, file="jimmy_lm1_10x_1019.cluster_markers.csv")

##save file
saveRDS(lm_obj, file = "jimmy_lm1_10x_1019.rds")

##and then load
retina.integrated <- readRDS(file ="jimmy_lm1_10x_1019.rds")