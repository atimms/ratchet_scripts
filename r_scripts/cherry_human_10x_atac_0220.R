library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
set.seed(1234)
setwd('/data/atimms/cherry_human_10x_atac_0220')

##install EnsDb.Hsapiens.v86
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("EnsDb.Hsapiens.v86")

##analysis sourced from:
##https://satijalab.org/signac/articles/pbmc_vignette.html
##https://satijalab.org/signac/articles/integration.html
##need 4x files; metadata, matrix, fragments and fragments index

##Human5
##load files
hu5_counts <- Read10X_h5(filename = "Hu5_filtered_peak_bc_matrix.h5")
hu5_metadata <- read.csv(file = "Hu5_singlecell.csv", header = TRUE, row.names = 1)
#hu5_data <- CreateSeuratObject(counts = hu5_counts, assay = 'peaks', project = 'ATAC', min.cells = 1, meta.data = hu5_metadata)
hu5_data <- CreateSeuratObject(counts = hu5_counts, assay = 'ATAC', project = 'hu5', min.cells = 1, meta.data = hu5_metadata)
##load fragment file and then filter for cells that are in the matrix
fragment.path <- 'Hu5_fragments.tsv.gz'
hu5_data <- SetFragments(object = hu5_data, file = fragment.path)
hu5_data
##filter fragment file -- not using as combining samples
#fragment_file_filtered <- "Hu5_filtered_fragments.tsv"
#FilterFragments(fragment.path = fragment.path, cells = colnames(hu5_data), output.path = fragment_file_filtered)
#hu5_data <- SetFragments(object = hu5_data, file = paste0(fragment_file_filtered, '.bgz'))
##qc metrics
hu5_data <- NucleosomeSignal(object = hu5_data)
hu5_data$pct_reads_in_peaks <- hu5_data$peak_region_fragments / hu5_data$passed_filters * 100
hu5_data$blacklist_ratio <- hu5_data$blacklist_region_fragments / hu5_data$peak_region_fragments
VlnPlot(object = hu5_data, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.Hu5.qc.pdf")
##fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength
hu5_data$nucleosome_group <- ifelse(hu5_data$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
FragmentHistogram(object = hu5_data, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.Hu5.frag_by_ns.pdf")
##look for enrichment aroung TSS
# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
# to save time use the first 2000 TSSs
hu5_data <- TSSEnrichment(object = hu5_data, tss.positions = tss.ranges[1:2000])
hu5_data$high.tss <- ifelse(hu5_data$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(hu5_data, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.Hu5.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
hu5_subset_data <- subset(hu5_data, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
hu5_data
hu5_subset_data
##Normalization and linear dimensional reduction
#hu5_subset_data$sample <- 'hu5'
#hu5_subset_data <- RunTFIDF(hu5_subset_data)
#hu5_subset_data <- FindTopFeatures(hu5_subset_data, min.cutoff = 50)
#hu5_subset_data <- RunSVD(hu5_subset_data, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
##Non-linear dimension reduction and clustering
#hu5_subset_data <- RunUMAP(hu5_subset_data, reduction = 'lsi', dims = 2:30)
#hu5_subset_data <- FindNeighbors(object = hu5_subset_data, reduction = 'lsi', dims = 2:30)
#hu5_subset_data <- FindClusters(object = hu5_subset_data, verbose = FALSE, algorithm = 3)
#DimPlot(object = hu5_subset_data, label = TRUE) + NoLegend()
#dev.copy2pdf(file="cherry_atac_0320.Hu5.UMAP.pdf")


##Hu7
##load files
hu7_counts <- Read10X_h5(filename = "Hu7_filtered_peak_bc_matrix.h5")
hu7_metadata <- read.csv(file = "Hu7_singlecell.csv", header = TRUE, row.names = 1)
#hu7_data <- CreateSeuratObject(counts = hu7_counts, assay = 'peaks', project = 'ATAC', min.cells = 1, meta.data = hu7_metadata)
hu7_data <- CreateSeuratObject(counts = hu7_counts, assay = 'ATAC', project = 'hu7', min.cells = 1, meta.data = hu7_metadata)
##load fragment file and then filter for cells that are in the matrix
fragment.path <- 'Hu7_fragments.tsv.gz'
hu7_data <- SetFragments(object = hu7_data, file = fragment.path)
hu7_data
##filter fragment file -- not using as combining samples
#fragment_file_filtered <- "Hu7_filtered_fragments.tsv"
#FilterFragments(fragment.path = fragment.path, cells = colnames(hu7_data), output.path = fragment_file_filtered)
#hu7_data <- SetFragments(object = hu7_data, file = paste0(fragment_file_filtered, '.bgz'))
##qc metrics
hu7_data <- NucleosomeSignal(object = hu7_data)
hu7_data$pct_reads_in_peaks <- hu7_data$peak_region_fragments / hu7_data$passed_filters * 100
hu7_data$blacklist_ratio <- hu7_data$blacklist_region_fragments / hu7_data$peak_region_fragments
VlnPlot(object = hu7_data, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.Hu7.qc.pdf")
##fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength
hu7_data$nucleosome_group <- ifelse(hu7_data$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
FragmentHistogram(object = hu7_data, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.Hu7.frag_by_ns.pdf")
##look for enrichment aroung TSS
# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
# to save time use the first 2000 TSSs
hu7_data <- TSSEnrichment(object = hu7_data, tss.positions = tss.ranges[1:2000])
hu7_data$high.tss <- ifelse(hu7_data$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(hu7_data, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.Hu7.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
hu7_subset_data <- subset(hu7_data, subset = peak_region_fragments > 2000 & peak_region_fragments < 40000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
hu7_data
hu7_subset_data
##Normalization and linear dimensional reduction
#hu7_subset_data$sample <- 'hu7'
#hu7_subset_data <- RunTFIDF(hu7_subset_data)
#hu7_subset_data <- FindTopFeatures(hu7_subset_data, min.cutoff = 50)
#hu7_subset_data <- RunSVD(hu7_subset_data, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
##Non-linear dimension reduction and clustering
#hu7_subset_data <- RunUMAP(hu7_subset_data, reduction = 'lsi', dims = 2:30)
#hu7_subset_data <- FindNeighbors(object = hu7_subset_data, reduction = 'lsi', dims = 2:30)
#hu7_subset_data <- FindClusters(object = hu7_subset_data, verbose = FALSE, algorithm = 3)
#DimPlot(object = hu7_subset_data, label = TRUE) + NoLegend()
#dev.copy2pdf(file="cherry_atac_0320.Hu7.UMAP.pdf")

##Hu8
##load files
hu8_counts <- Read10X_h5(filename = "Hu8_filtered_peak_bc_matrix.h5")
hu8_metadata <- read.csv(file = "Hu8_singlecell.csv", header = TRUE, row.names = 1)
#hu8_data <- CreateSeuratObject(counts = hu8_counts, assay = 'peaks', project = 'ATAC', min.cells = 1, meta.data = hu8_metadata)
hu8_data <- CreateSeuratObject(counts = hu8_counts, assay = 'ATAC', project = 'hu8', min.cells = 1, meta.data = hu8_metadata)
##load fragment file and then filter for cells that are in the matrix
fragment.path <- 'Hu8_fragments.tsv.gz'
hu8_data <- SetFragments(object = hu8_data, file = fragment.path)
hu8_data
##filter fragment file -- not using as combining samples
#fragment_file_filtered <- "Hu8_filtered_fragments.tsv"
#FilterFragments(fragment.path = fragment.path, cells = colnames(hu8_data), output.path = fragment_file_filtered)
#hu8_data <- SetFragments(object = hu8_data, file = paste0(fragment_file_filtered, '.bgz'))
##qc metrics
hu8_data <- NucleosomeSignal(object = hu8_data)
hu8_data$pct_reads_in_peaks <- hu8_data$peak_region_fragments / hu8_data$passed_filters * 100
hu8_data$blacklist_ratio <- hu8_data$blacklist_region_fragments / hu8_data$peak_region_fragments
VlnPlot(object = hu8_data, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.Hu8.qc.pdf")
##fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength
hu8_data$nucleosome_group <- ifelse(hu8_data$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
FragmentHistogram(object = hu8_data, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.Hu8.frag_by_ns.pdf")
##look for enrichment aroung TSS
# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
# to save time use the first 2000 TSSs
hu8_data <- TSSEnrichment(object = hu8_data, tss.positions = tss.ranges[1:2000])
hu8_data$high.tss <- ifelse(hu8_data$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(hu8_data, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.Hu8.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
hu8_subset_data <- subset(hu8_data, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
hu8_data
hu8_subset_data
##Normalization and linear dimensional reduction
hu8_subset_data$sample <- 'hu8'
hu8_subset_data <- RunTFIDF(hu8_subset_data)
hu8_subset_data <- FindTopFeatures(hu8_subset_data, min.cutoff = 50)
hu8_subset_data <- RunSVD(hu8_subset_data, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
##Non-linear dimension reduction and clustering
hu8_subset_data <- RunUMAP(hu8_subset_data, reduction = 'lsi', dims = 2:30)
hu8_subset_data <- FindNeighbors(object = hu8_subset_data, reduction = 'lsi', dims = 2:30)
hu8_subset_data <- FindClusters(object = hu8_subset_data, verbose = FALSE, algorithm = 3)
DimPlot(object = hu8_subset_data, label = TRUE) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.Hu8.UMAP.pdf")
# find markers for every cluster compared to all remaining cells, report only the positive ones
hu.integrated.markers <- FindAllMarkers(hu8_subset_data, only.pos = TRUE, min.pct = 0.4, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(hu.integrated.markers, file="cherry_atac_0320.hu8.cluster_markers.csv")
## get genes closest to marker genes, row names are weird in a few places
open_regions <- hu.integrated.markers$gene
closest_genes <- ClosestFeature(regions = open_regions, annotation = EnsDb.Hsapiens.v86,
                                sep = c(':', '-'))
write.csv(closest_genes, file="cherry_atac_0320.hu8.cluster_marker_closest_genes.csv")


## create a common peak set
combined.peaks <- UnifyPeaks(object.list = list(hu5_subset_data, hu7_subset_data, hu8_subset_data), mode = "reduce")
combined.peaks
## quantify peaks in each dataset
hu5.counts <- FeatureMatrix(fragments = GetFragments(hu5_subset_data), features = combined.peaks,
  sep = c(":", "-"), cells = colnames(hu5_subset_data))
hu7.counts <- FeatureMatrix(fragments = GetFragments(hu7_subset_data), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(hu7_subset_data))
hu8.counts <- FeatureMatrix(fragments = GetFragments(hu8_subset_data), features = combined.peaks,
                            sep = c(":", "-"), cells = colnames(hu8_subset_data))
hu5_subset_data[['peaks']] <- CreateAssayObject(counts = hu5.counts)
hu7_subset_data[['peaks']] <- CreateAssayObject(counts = hu7.counts)
hu8_subset_data[['peaks']] <- CreateAssayObject(counts = hu8.counts)
##Merge objects... is it working?
# add information to identify dataset of origin
hu5_subset_data$dataset <- 'hu5'
hu7_subset_data$dataset <- 'hu7'
hu8_subset_data$dataset <- 'hu8'

##merge data
# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(x = hu5_subset_data, y = list(hu7_subset_data, hu8_subset_data), add.cell.ids = c("hu5", "hu7", "hu8"))
# make sure to change to the assay containing common peaks
DefaultAssay(combined) <- "peaks" 
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff ="q75")
combined <- RunSVD(combined, reduction.key = 'LSI_',
  reduction.name = 'lsi', irlba.work = 400)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="cherry_atac_0320.hu_merged.umap_by_sample.pdf")
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3,resolution = 0.8)
DimPlot(object = combined, label = TRUE, pt.size = 0.1)
dev.copy2pdf(file="cherry_atac_0320.hu_merged.umap.pdf")

##save/load merged object..
#saveRDS(combined, file = "cherry_atac_0320.hu_merged.rds")
#combined <- readRDS(file = "cherry_atac_0320.hu_merged.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
hu.integrated.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.4, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(hu.integrated.markers, file="cherry_atac_0320.hu_merged.cluster_markers.csv")
## get genes closest to marker genes, row names are weird in a few places
open_regions <- hu.integrated.markers$gene
closest_genes <- ClosestFeature(regions = open_regions, annotation = EnsDb.Hsapiens.v86,
                                sep = c(':', '-'))
write.csv(closest_genes, file="cherry_atac_0320.hu_merged.cluster_marker_closest_genes.csv")


##integrate data
##make sure using correct assay as default 
DefaultAssay(hu5_subset_data) <- "peaks" 
DefaultAssay(hu7_subset_data) <- "peaks" 
DefaultAssay(hu8_subset_data) <- "peaks" 

# find integration anchors between 10x and sci-ATAC
anchors <- FindIntegrationAnchors(object.list = list(hu5_subset_data,hu7_subset_data, hu8_subset_data),
  k.filter = NA)
# integrate data and create a new merged object
integrated <- IntegrateData(anchorset = anchors,preserve.order = TRUE)
# we now have a "corrected" TF-IDF matrix??, and can run LSI again on this corrected matrix
integrated <- RunSVD(object = integrated, n = 30, reduction.name = 'integratedLSI')
##what dims to use for umap?, test..
integrated <- RunUMAP(object = integrated,dims = 2:30, reduction = 'integratedLSI')
DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)
dev.copy2pdf(file="cherry_atac_0320.hu_integrated.umap_by_sample.pdf")
integrated <- FindNeighbors(object = integrated, reduction = 'integratedLSI', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3,resolution = 0.4)
DimPlot(object = integrated, pt.size = 0.1)
dev.copy2pdf(file="cherry_atac_0320.hu_integrated.umap.pdf")


DefaultAssay(integrated) <- "peaks" 

# find markers for every cluster compared to all remaining cells, report only the positive ones
hu.integrated.markers.second <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.4, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(hu.integrated.markers.second, file="cherry_atac_0320.hu_integrated.cluster_markers.csv")
## get genes closest to marker genes, row names are weird in a few places
open_regions.second <- hu.integrated.markers.second$gene
closest_genes.second <- ClosestFeature(regions = open_regions.second, annotation = EnsDb.Hsapiens.v86,
                                sep = c(':', '-'))
write.csv(closest_genes.second, file="cherry_atac_0320.hu_integrated.cluster_marker_closest_genes.csv")



