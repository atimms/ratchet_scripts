library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
set.seed(1234)

## Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC')

## load data if needed
d53_subset <- readRDS(file = "./signac_analysis/d53_subset.rds")
d74_subset <- readRDS(file = "./signac_analysis/d74_subset.rds")
d78_subset <- readRDS(file = "./signac_analysis/d78_subset.rds")
Hu5_subset <- readRDS(file = "./signac_analysis/Hu5_subset.rds")
Hu7_subset <- readRDS(file = "./signac_analysis/Hu7_subset.rds")
Hu8_subset <- readRDS(file = "./signac_analysis/Hu8_subset.rds")

##merge data
## create a common peak set
combined.peaks <- UnifyPeaks(object.list = list(d53_subset,d74_subset,d78_subset,Hu5_subset,Hu7_subset,Hu8_subset), mode = "reduce")
combined.peaks

## quantify peaks in each dataset, and add peaks assay to object
d53.counts <- FeatureMatrix(fragments = GetFragments(d53_subset), features = combined.peaks,
                                    sep = c(":", "-"), cells = colnames(d53_subset))
d53_subset[['peaks']] <- CreateAssayObject(counts = d53.counts)
d74.counts <- FeatureMatrix(fragments = GetFragments(d74_subset), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(d74_subset))
d74_subset[['peaks']] <- CreateAssayObject(counts = d74.counts)
d78.counts <- FeatureMatrix(fragments = GetFragments(d78_subset), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(d78_subset))
d78_subset[['peaks']] <- CreateAssayObject(counts = d78.counts)
Hu5.counts <- FeatureMatrix(fragments = GetFragments(Hu5_subset), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(Hu5_subset))
Hu5_subset[['peaks']] <- CreateAssayObject(counts = Hu5.counts)
Hu7.counts <- FeatureMatrix(fragments = GetFragments(Hu7_subset), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(Hu7_subset))
Hu7_subset[['peaks']] <- CreateAssayObject(counts = Hu7.counts)
Hu8.counts <- FeatureMatrix(fragments = GetFragments(Hu8_subset), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(Hu8_subset))
Hu8_subset[['peaks']] <- CreateAssayObject(counts = Hu8.counts)

# add information to identify dataset of origin
d53_subset$dataset <- 'd53'
d74_subset$dataset <- 'd74'
d78_subset$dataset <- 'd78'
Hu5_subset$dataset <- 'Hu5'
Hu7_subset$dataset <- 'Hu7'
Hu8_subset$dataset <- 'Hu8'

# make sure to change to the assay containing common peaks
DefaultAssay(d53_subset) <- "peaks" 
DefaultAssay(d74_subset) <- "peaks" 
DefaultAssay(d78_subset) <- "peaks" 
DefaultAssay(Hu5_subset) <- "peaks" 
DefaultAssay(Hu7_subset) <- "peaks" 
DefaultAssay(Hu8_subset) <- "peaks" 

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(x = d53_subset, y = list(d74_subset, d78_subset,Hu5_subset,Hu7_subset,Hu8_subset), add.cell.ids = c("d53", "d74", "d78", "Hu5", "Hu7", "Hu8"))

# make sure to change to the assay containing common peaks
DefaultAssay(combined) <- "peaks" 
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined, reduction.key = 'LSI_',
                   reduction.name = 'lsi', irlba.work = 400)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)
DimPlot(object = combined, label = TRUE)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.merged.UMAP_by_cluster.pdf", width=20)
DimPlot(combined, split.by = 'dataset', pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.merged.UMAP_split_by_sample.pdf", width=20)

##merged fragments using unix commands, so add to merged object
combined <- SetFragments(combined, "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments/fragments.tsv.gz")


# save individual samples with peaks assay and combined object
saveRDS(d53_subset, file = "./signac_analysis/d53_combPeaks.rds")
saveRDS(d74_subset, file = "./signac_analysis/d74_combPeaks.rds")
saveRDS(d78_subset, file = "./signac_analysis/d78_combPeaks.rds")
saveRDS(Hu5_subset, file = "./signac_analysis/Hu5_combPeaks.rds")
saveRDS(Hu7_subset, file = "./signac_analysis/Hu7_combPeaks.rds")
saveRDS(Hu8_subset, file = "./signac_analysis/Hu8_combPeaks.rds")
saveRDS(combined, file = "./signac_analysis/combined.rds")
saveRDS(combined.peaks, file = "./signac_analysis/combined_peaks.rds")
