library(dplyr)
library(Seurat)
setwd('/data/atimms/denise_cr_0919')

##need this to use umap
#reticulate::py_install(packages ='umap-learn')

# Load the 10x dataset
dcf.data <- Read10X(data.dir = "/data/atimms/denise_cr_0919/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
dcf <- CreateSeuratObject(counts = dcf.data, project = "dc_filtered", min.cells = 3, min.features = 200)
dcf

# store mitochondrial percentage in object meta data
dcf <- PercentageFeatureSet(dcf, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
dcf <- SCTransform(dcf, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
dcf <- RunPCA(dcf, verbose = FALSE)
dcf <- RunUMAP(dcf, dims = 1:30, verbose = FALSE)

dcf <- FindNeighbors(dcf, dims = 1:30, verbose = FALSE)
dcf <- FindClusters(dcf, verbose = FALSE)
DimPlot(dcf, label = TRUE) + NoLegend()
dev.copy2pdf(file="denise_10x_1019.umap.pdf", width = 10)


# find markers for every cluster compared to all remaining cells, report only the positive ones
lm_obj.markers <- FindAllMarkers(lm_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(lm_obj.markers, file="denise_10x1019.cluster_markers.csv")

##save file
saveRDS(lm_obj, file = "denise_10x_1019.rds")

##and then load
retina.integrated <- readRDS(file ="denise_10x_1019.rds")