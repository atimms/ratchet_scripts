library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
set.seed(1234)
setwd('/data/atimms/cherry_human_10x_atac_0220')



##Human data combined using aggr
##load files
human_counts <- Read10X_h5(filename = "human_filtered_peak_bc_matrix.h5")
human_metadata <- read.csv(file = "human_singlecell.csv", header = TRUE, row.names = 1)
#human_data <- CreateSeuratObject(counts = human_counts, assay = 'peaks', project = 'ATAC', min.cells = 1, meta.data = human_metadata)
human_data <- CreateSeuratObject(counts = human_counts, assay = 'ATAC', project = 'human', min.cells = 1, meta.data = human_metadata, names.field = 2, names.delim = "-")
##load fragment file and then filter for cells that are in the matrix
fragment.path <- 'human_fragments.tsv.gz'
human_data <- SetFragments(object = human_data, file = fragment.path)
human_data
##filter fragment file
fragment_file_filtered <- "human_filtered_fragments.tsv"
FilterFragments(fragment.path = fragment.path, cells = colnames(human_data), output.path = fragment_file_filtered)
human_data <- SetFragments(object = human_data, file = paste0(fragment_file_filtered, '.bgz'))
##qc metrics
human_data <- NucleosomeSignal(object = human_data)
human_data$pct_reads_in_peaks <- human_data$peak_region_fragments / human_data$passed_filters * 100
human_data$blacklist_ratio <- human_data$blacklist_region_fragments / human_data$peak_region_fragments
VlnPlot(object = human_data, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.human.qc.pdf")
##fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength
human_data$nucleosome_group <- ifelse(human_data$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
FragmentHistogram(object = human_data, group.by = 'nucleosome_group')
dev.copy2pdf(file="cherry_atac_0320.human.frag_by_ns.pdf")
##look for enrichment aroung TSS
# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
# to save time use the first 2000 TSSs
human_data <- TSSEnrichment(object = human_data, tss.positions = tss.ranges[1:2000])
human_data$high.tss <- ifelse(human_data$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(human_data, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.human.tss_enrichment.pdf")
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
human_subset_data <- subset(human_data, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 10 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
human_data
human_subset_data
##Normalization and linear dimensional reduction
human_subset_data$sample <- 'human'
human_subset_data <- RunTFIDF(human_subset_data)
human_subset_data <- FindTopFeatures(human_subset_data, min.cutoff = 50)
human_subset_data <- RunSVD(human_subset_data, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
##Non-linear dimension reduction and clustering
human_subset_data <- RunUMAP(human_subset_data, reduction = 'lsi', dims = 2:30)
human_subset_data <- FindNeighbors(object = human_subset_data, reduction = 'lsi', dims = 2:30)
human_subset_data <- FindClusters(object = human_subset_data, verbose = FALSE, algorithm = 3)
DimPlot(object = human_subset_data, label = TRUE) + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.human.UMAP.pdf")
DimPlot(object = human_subset_data, label = TRUE, group.by = "orig.ident") + NoLegend()
dev.copy2pdf(file="cherry_atac_0320.human.UMAP_by_sample.pdf")
