library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(harmony)
set.seed(1234)

##info from:
##https://satijalab.org/signac/articles/pbmc_vignette.html
##https://satijalab.org/signac/articles/integration.html

## Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC')

## load data if needed - from merged data
combined <- readRDS(file = "./signac_analysis/combined.rds")

hm.integrated <- RunHarmony(object = combined, group.by.vars = 'dataset',
  reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)

# re-compute the UMAP using corrected LSI embeddings
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
#DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1, label = TRUE)
#dev.copy2pdf(file="./signac_analysis/human_scatac_0420.harmony.UMAP_by_sample.pdf", width=20)
#DimPlot(hm.integrated, split.by = 'dataset', pt.size = 0.1, label = TRUE)
#dev.copy2pdf(file="./signac_analysis/human_scatac_0420.harmony.UMAP_split_by_sample.pdf", width=20)
##find clusters for 2 different resolutions
hm.integrated <- FindNeighbors(hm.integrated, reduction = 'harmony', dims = 2:30)
##res 0.8
hm.integrated <- FindClusters(hm.integrated, verbose = FALSE, algorithm = 3, resolution = 0.8)
DimPlot(object = hm.integrated, pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.harmony_res0.8.UMAP_by_cluster.pdf", width=20)
DimPlot(hm.integrated, split.by = 'dataset', pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.harmony_res0.8.UMAP_split_by_sample.pdf", width=20)
##res 0.4
hm.integrated <- FindClusters(hm.integrated, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(object = hm.integrated, pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.harmony_res0.4.UMAP_by_cluster.pdf", width=20)
DimPlot(hm.integrated, split.by = 'dataset', pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.harmony_res0.4.UMAP_split_by_sample.pdf", width=20)



##save/load file
#saveRDS(hm.integrated, file = "./signac_analysis/harmony_integrated.rds")
#hm.integrated <- readRDS(file = "./signac_analysis/harmony_integrated.rds")

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
  cells = colnames(hm.integrated),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
hm.integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
hm.integrated <- NormalizeData(
  object = hm.integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm.integrated$nCount_RNA)
)

##visualize canonical marker 
DefaultAssay(hm.integrated) <- 'RNA'
FeaturePlot(object = hm.integrated, features = c('SLC1A3', 'PROM1', 'RLBP1', 'ONECUT1', 'VSX1', 'SLC6A9', 'POU4F2', 'NGFR', 'C1QA', 'GJA4', 'BEST1', 'VIM', 'ATOH7', 'PRDM13', 'DLL3'), pt.size = 0.1, max.cutoff = 'q95', ncol = 4)
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.harmony.feature_plot1.res0.4.pdf", width=20)

##save/load file
#saveRDS(hm.integrated, file = "./signac_analysis/harmony_integrated.rds")
hm.integrated <- readRDS(file = "./signac_analysis/harmony_integrated.rds")

##find markers
# switch back to working with peaks instead of gene activities
DefaultAssay(hm.integrated) <- 'peaks'
harmony.markers <- FindAllMarkers(hm.integrated, only.pos = TRUE, min.pct = 0.25, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(harmony.markers, file='./signac_analysis/human_scatac_0420.harmony.marker.res0.4.csv')
head(harmony.markers)

## get genes closest to marker genes, row names are weird in a few places
##temp way to lead results
#harmony.markers <- read.csv('./signac_analysis/human_scatac_0420.harmony.marker.res0.4.csv')
##list of peaks
open_regions <- harmony.markers$gene
head(open_regions)

##get gene locations.. code copied from earlier in script
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

closest_genes <- ClosestFeature(regions = open_regions, annotation = gene.ranges,
                                       sep = c(':', '-'))
write.csv(closest_genes, file='./signac_analysis/human_scatac_0420.harmony.marker_gene_ids.res0.4.csv')

