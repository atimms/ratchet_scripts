library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
set.seed(1234)

##info from:
#https://satijalab.org/signac/articles/integration.html
#https://github.com/timoast/signac/issues/89
#https://github.com/timoast/signac/issues/94

## Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC')

## load data if needed - objects with combined peaks assay and the peaks from merged analysis
combined <- readRDS(file = "./signac_analysis/combined.rds")


# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
combined_split <- SplitObject(combined, split.by = "dataset")

##get peaks to use for integration
##just get peaks from object (because the combined peaks are active assay)
peaks.use <- sample(rownames(combined), size = 10000, replace = FALSE) 

# find integration anchors 
anchors <- FindIntegrationAnchors(object.list = combined_split, anchor.features = peaks.use, k.filter = NA)

##save/load data
#saveRDS(peaks.use, file = "./signac_analysis/peaks.use.rds")
#saveRDS(anchors, file = "./signac_analysis/anchors.rds")
#saveRDS(combined_split, file = "./signac_analysis/combined_split.rds")
combined_split <- readRDS(file = "./signac_analysis/combined_split.rds")
anchors <- readRDS(file = "./signac_analysis/anchors.rds")
peaks.use <- readRDS(file = "./signac_analysis/peaks.use.rds")

# integrate data and create a new merged object
integrated <- IntegrateData(anchorset = anchors, dims = 2:30, preserve.order = TRUE)

# we now have a "corrected" TF-IDF matrix??, and can run LSI again on this corrected matrix
integrated <- RunSVD(object = integrated, n = 30, reduction.name = 'integratedLSI')
##what dims to use for umap?, test..
integrated <- RunUMAP(object = integrated, dims = 2:30, reduction = 'integratedLSI')
DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_sample.pdf", width=20)
integrated <- FindNeighbors(object = integrated, reduction = 'integratedLSI', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(object = integrated, pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_cluster.res0.4.pdf", width=20)
DimPlot(integrated, split.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_sample.res0.4.pdf", width=20)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution = 0.8)
DimPlot(object = integrated, pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_cluster.res0.8.pdf", width=20)
DimPlot(integrated, split.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.UMAP_by_sample.res0.8.pdf", width=20)


##save/load integrated data
#saveRDS(integrated, file = "./signac_analysis/integrated.rds")
#integrated <- readRDS(file = "./signac_analysis/integrated.rds")

##make gene matrix and add gene names instead of peaks

# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
fragment.path = "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments/fragments.tsv.gz"
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(integrated),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)

##visualize canonical marker 
DefaultAssay(integrated) <- 'RNA'
FeaturePlot(object = integrated, features = c('SLC1A3', 'PROM1', 'RLBP1', 'ONECUT1', 'VSX1', 'SLC6A9', 'POU4F2', 'NGFR', 'C1QA', 'GJA4', 'BEST1', 'VIM', 'ATOH7', 'PRDM13', 'DLL3'), pt.size = 0.1, max.cutoff = 'q95', ncol = 4)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.integrated.feature_plot1.res0.8.pdf", width=20)

##save/load integrated data
saveRDS(integrated, file = "./signac_analysis/integrated.rds")
integrated <- readRDS(file = "./signac_analysis/integrated.rds")
