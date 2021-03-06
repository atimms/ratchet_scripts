library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomeInfoDb)
library(ggplot2)
library(harmony)
set.seed(1234)

# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC')

#create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')


# define a convenient function to load all the data and create a Seurat object
create_obj <- function(dir, sample_id) {
  count.path <- list.files(path = dir, pattern = "filtered_peak_bc_matrix.h5", full.names = TRUE)
  fragment.path <- list.files(path = dir, pattern = "fragments.tsv.gz", full.names = TRUE)[1]
  counts <- Read10X_h5(count.path)
  md.path <- list.files(path = dir, pattern = "singlecell.csv", full.names = TRUE)
  #md <- read.table(file = md.path, stringsAsFactors = FALSE, sep = ",", header = TRUE)
  md <- read.csv(file = md.path, header = TRUE, row.names = 1)
  obj <- CreateSeuratObject(counts = counts, assay = "ATAC", meta.data = md, project = sample_id)
  obj <- SetFragments(obj, file = fragment.path)
  ##qc metrics
  obj <- NucleosomeSignal(object = obj)
  obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
  obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments
  ##fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength
  obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
  ##look for enrichment aroung TSS
  obj <- TSSEnrichment(object = obj, tss.positions = tss.ranges[1:2000])
  obj$high.tss <- ifelse(obj$TSS.enrichment > 2, 'High', 'Low')
  return(obj)
}

##create object for each sample and gather data
d74 <- create_obj("./d74/outs", "d74")
d78 <- create_obj("./d78/outs", "d78")

##graph QC and subset, then normalize etc
##d74
VlnPlot(object = d74, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.d74.qc.pdf", width=20)
FragmentHistogram(object = d74, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.d74.frag_by_ns.pdf", width=20)
TSSPlot(d74, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.d74.tss_enrichment.pdf", width=20)
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
d74_subset <- subset(d74, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
d74
d74_subset
##Normalization and linear dimensional reduction
d74_subset <- RunTFIDF(d74_subset)
d74_subset <- FindTopFeatures(d74_subset, min.cutoff = 'q0')
d74_subset <- RunSVD(object = d74_subset, assay = 'ATAC',
                     reduction.key = 'LSI_', reduction.name = 'lsi')
##correlation between depth and dims (may need to remove dim1 from analysis)
DepthCor(d74_subset)
##Non-linear dimension reduction and clustering
d74_subset <- RunUMAP(object = d74_subset, reduction = 'lsi', dims = 2:30)
d74_subset <- FindNeighbors(object = d74_subset, reduction = 'lsi', dims = 2:30)
d74_subset <- FindClusters(object = d74_subset, verbose = FALSE, algorithm = 3)
DimPlot(object = d74_subset, label = TRUE) + NoLegend()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.UMAP.pdf", width=20)

##d78
VlnPlot(object = d78, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.d78.qc.pdf", width=20)
FragmentHistogram(object = d78, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.d78.frag_by_ns.pdf", width=20)
TSSPlot(d78, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.d78.tss_enrichment.pdf", width=20)
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
d78_subset <- subset(d78, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
d78
d78_subset
##Normalization and linear dimensional reduction
d78_subset <- RunTFIDF(d78_subset)
d78_subset <- FindTopFeatures(d78_subset, min.cutoff = 'q0')
d78_subset <- RunSVD(object = d78_subset, assay = 'ATAC',
                     reduction.key = 'LSI_', reduction.name = 'lsi')
##correlation between depth and dims (may need to remove dim1 from analysis)
DepthCor(d78_subset)
##Non-linear dimension reduction and clustering
d78_subset <- RunUMAP(object = d78_subset, reduction = 'lsi', dims = 2:30)
d78_subset <- FindNeighbors(object = d78_subset, reduction = 'lsi', dims = 2:30)
d78_subset <- FindClusters(object = d78_subset, verbose = FALSE, algorithm = 3)
DimPlot(object = d78_subset, label = TRUE) + NoLegend()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0620.UMAP.pdf", width=20)


##merge data
## create a common peak set
combined.peaks <- UnifyPeaks(object.list = list(d74_subset,d78_subset), mode = "reduce")
combined.peaks

## quantify peaks in each dataset, and add peaks assay to object
d74.counts <- FeatureMatrix(fragments = GetFragments(d74_subset), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(d74_subset))
d74_subset[['peaks']] <- CreateAssayObject(counts = d74.counts)
d78.counts <- FeatureMatrix(fragments = GetFragments(d78_subset), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(d78_subset))
d78_subset[['peaks']] <- CreateAssayObject(counts = d78.counts)

# add information to identify dataset of origin
d74_subset$dataset <- 'd74'
d78_subset$dataset <- 'd78'

# make sure to change to the assay containing common peaks
DefaultAssay(d74_subset) <- "peaks" 
DefaultAssay(d78_subset) <- "peaks" 

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(x = d74_subset, y = list(d78_subset), add.cell.ids = c("d74", "d78"))

# make sure to change to the assay containing common peaks
DefaultAssay(combined) <- "peaks" 
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined, reduction.key = 'LSI_',
                   reduction.name = 'lsi', irlba.work = 400)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(object = combined, label = TRUE)
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0520.merged.res_0.4.UMAP_by_cluster.pdf", width=20)
DimPlot(combined, split.by = 'dataset', pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0520.merged.res_0.4.UMAP_split_by_sample.pdf", width=20)

##merged fragments using unix commands, so add to merged object
combined <- SetFragments(combined, "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments_d70s/fragments.tsv.gz")

# save individual samples with peaks assay and combined object
saveRDS(d74_subset, file = "./signac_analysis/d70s/d74_combPeaks.rds")
saveRDS(d78_subset, file = "./signac_analysis/d70s/d78_combPeaks.rds")
saveRDS(combined, file = "./signac_analysis/d70s/combined.rds")
saveRDS(combined.peaks, file = "./signac_analysis/d70s/combined_peaks.rds")

##load file if needed
#combined <- readRDS(file = "./signac_analysis/d70s/combined.rds")

##find markers
# switch back to working with peaks instead of gene activities
DefaultAssay(combined) <- 'peaks'
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(combined.markers, file='./signac_analysis/d70s/human_scatac_d70s_0520.merged.markers.csv')
head(combined.markers)

## get genes closest to marker genes, row names are weird in a few places
##temp way to lead results
#harmony.markers <- read.csv('./signac_analysis/human_scatac_0420.harmony.marker.res0.4.csv')
##list of peaks
open_regions <- combined.markers$gene
head(open_regions)

##get gene locations.. code copied from earlier in script
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')
closest_genes <- ClosestFeature(regions = open_regions, annotation = gene.ranges,
                                sep = c(':', '-'))
write.csv(closest_genes, file='./signac_analysis/d70s/human_scatac_d70s_0520.merged.marker_gene_ids.csv')


##make gene matrix and add gene names instead of peaks

# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
fragment.path = "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments_d70s/fragments.tsv.gz"
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(combined),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)

##make dotplots
DefaultAssay(combined) <- 'RNA'
#Dot plot - the size of the dot = % of cells and color represents the average expression
markers <- c('RCVRN', 'RHO', 'CRX', 'ARR3', 'GNAT2', 'VSX2', 'LHX4', 'TRPM1', 'GRM6', 'SLC1A3', 'RLBP1', 'PAX6', 'LHX1', 'ONECUT2', 'TFAP2B', 'GAD1', 'SLC6A9', 'RBPMS', 'NEFM', 'GFAP', 'CD74', 'P2RY12', 'BEST1', 'RPE65', 'SFRP2')
DotPlot(combined, features = markers) + RotatedAxis()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0520.merged.known_marker.dotplot.res0.4.pdf", width = 20)



##run harmony

## load data if needed - from merged data
combined <- readRDS(file = "./signac_analysis/d70s/combined.rds")

hm.integrated <- RunHarmony(object = combined, group.by.vars = 'dataset',
                            reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)

# re-compute the UMAP using corrected LSI embeddings
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')

##find clusters for 2 different resolutions
hm.integrated <- FindNeighbors(hm.integrated, reduction = 'harmony', dims = 2:30)
##res 0.4
hm.integrated <- FindClusters(hm.integrated, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(object = hm.integrated, pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0520.harmony_res0.4.UMAP_by_cluster.pdf", width=20)
DimPlot(hm.integrated, split.by = 'dataset', pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0520.harmony_res0.4.UMAP_split_by_sample.pdf", width=20)

##save/load file
saveRDS(hm.integrated, file = "./signac_analysis/d70s/harmony_integrated.rds")
#hm.integrated <- readRDS(file = "./signac_analysis/d70s/harmony_integrated.rds")

##make gene matrix and add gene names instead of peaks

# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
fragment.path = "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments_d70s/fragments.tsv.gz"
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
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0520.harmony.feature_plot1.res0.4.pdf", width=20)

##save/load file
saveRDS(hm.integrated, file = "./signac_analysis/d70s/harmony_integrated.rds")
hm.integrated <- readRDS(file = "./signac_analysis/d70s/harmony_integrated.rds")

##find markers
# switch back to working with peaks instead of gene activities
DefaultAssay(hm.integrated) <- 'peaks'
harmony.markers <- FindAllMarkers(hm.integrated, only.pos = TRUE, min.pct = 0.25, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(harmony.markers, file='./signac_analysis/d70s/human_scatac_d70s_0520.harmony.marker.res0.4.csv')
head(harmony.markers)

## get genes closest to marker genes, row names are weird in a few places
##temp way to lead results
#harmony.markers <- read.csv('./signac_analysis/d70s/human_scatac_d70s_0520.harmony.marker.res0.4.csv')
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
write.csv(closest_genes, file='./signac_analysis/d70s/human_scatac_d70s_0520.harmony.marker_gene_ids.res0.4.csv')

##make dotplots
DefaultAssay(hm.integrated) <- 'RNA'
#Dot plot - the size of the dot = % of cells and color represents the average expression
markers <- c('RCVRN', 'RHO', 'CRX', 'ARR3', 'GNAT2', 'VSX2', 'LHX4', 'TRPM1', 'GRM6', 'SLC1A3', 'RLBP1', 'PAX6', 'LHX1', 'ONECUT2', 'TFAP2B', 'GAD1', 'SLC6A9', 'RBPMS', 'NEFM', 'GFAP', 'CD74', 'P2RY12', 'BEST1', 'RPE65', 'SFRP2')
DotPlot(hm.integrated, features = markers) + RotatedAxis()
dev.copy2pdf(file="./signac_analysis/d70s/human_scatac_d70s_0520.harmony.known_marker.dotplot.res0.4.pdf", width = 20)




