library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
set.seed(1234)

##info from:
#https://satijalab.org/signac/articles/integration.html
#https://github.com/timoast/signac/issues/89
#https://github.com/timoast/signac/issues/94

## Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC')

## load data if needed - objects with combined peaks assay and the peaks from merged analysis
d53 <- readRDS(file = "./signac_analysis/d53_combPeaks.rds")
d74 <- readRDS(file = "./signac_analysis/d74_combPeaks.rds")
d78 <- readRDS(file = "./signac_analysis/d78_combPeaks.rds")
Hu5 <- readRDS(file = "./signac_analysis/Hu5_combPeaks.rds")
Hu7 <- readRDS(file = "./signac_analysis/Hu7_combPeaks.rds")
Hu8 <- readRDS(file = "./signac_analysis/Hu8_combPeaks.rds")
combined.peaks <- readRDS(file = "./signac_analysis/combined_peaks.rds")

##get peaks to use for integration
##just get peaks from object (because the combined peaks are active assay)
peaks.use <- sample(rownames(Hu8), size = 10000, replace = FALSE) 

# find integration anchors 
anchors <- FindIntegrationAnchors(object.list = list(d53,d74,d78,
  Hu5,Hu7,Hu8), anchor.features = peaks.use,
  assay = c('peaks', 'peaks','peaks', 'peaks','peaks', 'peaks'), k.filter = NA)

# integrate data and create a new merged object
integrated <- IntegrateData(anchorset = anchors, dims = 2:30, preserve.order = TRUE)

# we now have a "corrected" TF-IDF matrix??, and can run LSI again on this corrected matrix
integrated <- RunSVD(object = integrated, n = 30, reduction.name = 'integratedLSI')
##what dims to use for umap?, test..
integrated <- RunUMAP(object = integrated, dims = 2:30, reduction = 'integratedLSI')
DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_sample.pdf", width=20)
DimPlot(integrated, split.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_sample.pdf", width=20)
integrated <- FindNeighbors(object = integrated, reduction = 'integratedLSI', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3,resolution = 0.4)
DimPlot(object = integrated, pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_cluster.pdf", width=20)

##save file
saveRDS(integrated, file = "./signac_analysis/integrated.rds")

