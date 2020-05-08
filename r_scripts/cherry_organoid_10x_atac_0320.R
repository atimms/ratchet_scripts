library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
set.seed(1234)
setwd('/data/atimms/cherry_organoid_10x_atac_0320')

#create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

# define a convenient function to load all the data and create a Seurat object
create_obj <- function(dir) {
  count.path <- list.files(path = dir, pattern = "*_filtered_peak_bc_matrix.h5", full.names = TRUE)
  fragment.path <- list.files(path = dir, pattern = "*_fragments.tsv.gz", full.names = TRUE)[1]
  counts <- Read10X_h5(count.path)
  md.path <- list.files(path = dir, pattern = "*_singlecell.csv", full.names = TRUE)
  #md <- read.table(file = md.path, stringsAsFactors = FALSE, sep = ",", header = TRUE)
  md <- read.csv(file = md.path, header = TRUE, row.names = 1)
  obj <- CreateSeuratObject(counts = counts, assay = "ATAC", meta.data = md)
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


##make seurat object and gather data
org_IPSC_c4 <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/IPSC_c4")
##graph QC etc
VlnPlot(object = org_IPSC_c4, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.IPSC_c4.qc.pdf")
FragmentHistogram(object = org_IPSC_c4, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.IPSC_c4.frag_by_ns.pdf")
TSSPlot(org_IPSC_c4, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.IPSC_c4.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_IPSC_c4_subset <- subset(org_IPSC_c4, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_IPSC_c4
org_IPSC_c4_subset

##make seurat object and gather data
org_IPSC_c5_1 <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/IPSC_c5_1")
##graph QC etc
VlnPlot(object = org_IPSC_c5_1, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.IPSC_c5_1.qc.pdf")
FragmentHistogram(object = org_IPSC_c5_1, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.IPSC_c5_1.frag_by_ns.pdf")
TSSPlot(org_IPSC_c5_1, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.IPSC_c5_1.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_IPSC_c5_1_subset <- subset(org_IPSC_c5_1, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_IPSC_c5_1
org_IPSC_c5_1_subset

##make seurat object and gather data
org_5wk <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/5wk")
##graph QC etc
VlnPlot(object = org_5wk, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.5wk.qc.pdf")
FragmentHistogram(object = org_5wk, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.5wk.frag_by_ns.pdf")
TSSPlot(org_5wk, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.5wk.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_5wk_subset <- subset(org_5wk, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_5wk
org_5wk_subset
##make seurat object and gather data
org_5wk_c5_1 <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/5wk_c5_1")
##graph QC etc
VlnPlot(object = org_5wk_c5_1, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.5wk_c5_1.qc.pdf")
FragmentHistogram(object = org_5wk_c5_1, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.5wk_c5_1.frag_by_ns.pdf")
TSSPlot(org_5wk_c5_1, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.5wk_c5_1.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_5wk_c5_1_subset <- subset(org_5wk_c5_1, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_5wk_c5_1
org_5wk_c5_1_subset

##make seurat object and gather data
org_20wk <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/20wk")
##graph QC etc
VlnPlot(object = org_20wk, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.20wk.qc.pdf")
FragmentHistogram(object = org_20wk, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.20wk.frag_by_ns.pdf")
TSSPlot(org_20wk, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.20wk.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_20wk_subset <- subset(org_20wk, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_20wk
org_20wk_subset

##make seurat object and gather data
org_20wk_c5_1 <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/20wk_c5_1")
##graph QC etc
VlnPlot(object = org_20wk_c5_1, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.20wk_c5_1.qc.pdf")
FragmentHistogram(object = org_20wk_c5_1, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.20wk_c5_1.frag_by_ns.pdf")
TSSPlot(org_20wk_c5_1, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.20wk_c5_1.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_20wk_c5_1_subset <- subset(org_20wk_c5_1, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_20wk_c5_1
org_20wk_c5_1_subset

##make seurat object and gather data
org_28_1 <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/28_1")
##graph QC etc
VlnPlot(object = org_28_1, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.28_1.qc.pdf")
FragmentHistogram(object = org_28_1, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.28_1.frag_by_ns.pdf")
TSSPlot(org_28_1, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.28_1.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_28_1_subset <- subset(org_28_1, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_28_1
org_28_1_subset

##make seurat object and gather data
org_28_2 <- create_obj("/data/atimms/cherry_organoid_10x_atac_0320/28_2")
##graph QC etc
VlnPlot(object = org_28_2, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.28_2.qc.pdf")
FragmentHistogram(object = org_28_2, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.28_2.frag_by_ns.pdf")
TSSPlot(org_28_2, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.28_2.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_28_2_subset <- subset(org_28_2, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
org_28_2
org_28_2_subset

##merge data
## create a common peak set
combined.peaks <- UnifyPeaks(object.list = list(org_IPSC_c4_subset,org_IPSC_c5_1_subset,org_5wk_subset,org_5wk_c5_1_subset,org_20wk_subset,org_20wk_c5_1_subset,org_28_1_subset,org_28_2_subset), mode = "reduce")
combined.peaks
## quantify peaks in each dataset
org_IPSC_c4.counts <- FeatureMatrix(fragments = GetFragments(org_IPSC_c4_subset), features = combined.peaks,
                                    sep = c(":", "-"), cells = colnames(org_IPSC_c4_subset))
org_IPSC_c5_1.counts <- FeatureMatrix(fragments = GetFragments(org_IPSC_c5_1_subset), features = combined.peaks,
                                      sep = c(":", "-"), cells = colnames(org_IPSC_c5_1_subset))
org_5wk.counts <- FeatureMatrix(fragments = GetFragments(org_5wk_subset), features = combined.peaks,
                                sep = c(":", "-"), cells = colnames(org_5wk_subset))
org_5wk_c5_1.counts <- FeatureMatrix(fragments = GetFragments(org_5wk_c5_1_subset), features = combined.peaks,
                                     sep = c(":", "-"), cells = colnames(org_5wk_c5_1_subset))
org_20wk.counts <- FeatureMatrix(fragments = GetFragments(org_20wk_subset), features = combined.peaks,
                                 sep = c(":", "-"), cells = colnames(org_20wk_subset))
org_20wk_c5_1.counts <- FeatureMatrix(fragments = GetFragments(org_20wk_c5_1_subset), features = combined.peaks,
                                      sep = c(":", "-"), cells = colnames(org_20wk_c5_1_subset))
org_28_1.counts <- FeatureMatrix(fragments = GetFragments(org_28_1_subset), features = combined.peaks,
                                 sep = c(":", "-"), cells = colnames(org_28_1_subset))
org_28_2.counts <- FeatureMatrix(fragments = GetFragments(org_28_2_subset), features = combined.peaks,
                                 sep = c(":", "-"), cells = colnames(org_28_2_subset))

org_IPSC_c4_subset[['peaks']] <- CreateAssayObject(counts = org_IPSC_c4.counts)
org_IPSC_c5_1_subset[['peaks']] <- CreateAssayObject(counts = org_IPSC_c5_1.counts)
org_5wk_subset[['peaks']] <- CreateAssayObject(counts = org_5wk.counts)
org_5wk_c5_1_subset[['peaks']] <- CreateAssayObject(counts = org_5wk_c5_1.counts)
org_20wk_subset[['peaks']] <- CreateAssayObject(counts = org_20wk.counts)
org_20wk_c5_1_subset[['peaks']] <- CreateAssayObject(counts = org_20wk_c5_1.counts)
org_28_1_subset[['peaks']] <- CreateAssayObject(counts = org_28_1.counts)
org_28_2_subset[['peaks']] <- CreateAssayObject(counts = org_28_2.counts)

##Merge objects... is it working?
# add information to identify dataset of origin
org_IPSC_c4_subset$dataset <- 'org_IPSC_c4'
org_IPSC_c5_1_subset$dataset <- 'org_IPSC_c5_1'
org_5wk_subset$dataset <- 'org_5wk'
org_5wk_c5_1_subset$dataset <- 'org_5wk_c5_1'
org_20wk_subset$dataset <- 'org_20wk'
org_20wk_c5_1_subset$dataset <- 'org_20wk_c5_1'
org_28_1_subset$dataset <- 'org_28_1'
org_28_2_subset$dataset <- 'org_28_2'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(x = org_IPSC_c4_subset, y = list(org_IPSC_c5_1_subset,org_5wk_subset,org_5wk_c5_1_subset,org_20wk_subset,org_20wk_c5_1_subset,org_28_1_subset,org_28_2_subset), add.cell.ids = c("org_IPSC_c4","org_IPSC_c5_1","org_5wk","org_5wk_c5_1","org_20wk","org_20wk_c5_1","org_28_1","org_28_2"))
# make sure to change to the assay containing common peaks
DefaultAssay(combined) <- "peaks" 
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined, reduction.key = 'LSI_',
                   reduction.name = 'lsi', irlba.work = 400)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="cherry_atac_0320.organoid_merged.umap_by_sample.pdf")
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)
DimPlot(object = combined, label = TRUE) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.organoid_merged.umap.pdf")

##save/load merged object..
saveRDS(combined, file = "cherry_atac_0320.organoid_merged.rds")
combined <- readRDS(file = "cherry_atac_0320.organoid_merged.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
organoid.integrated.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.2, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(organoid.integrated.markers, file="cherry_atac_0320.organoid_merged.cluster_markers.csv")
## get genes closest to marker genes, row names are weird in a few places
open_regions <- organoid.integrated.markers$gene
#open_regions <- c("chr19:48043130-48045460", "chr5:30197484-30202939")
closest_genes <- ClosestFeature(regions = open_regions, annotation = EnsDb.Hsapiens.v86,
  sep = c(':', '-'))
write.csv(closest_genes, file="cherry_atac_0320.organoid_merged.cluster_marker_closest_genes.csv")


