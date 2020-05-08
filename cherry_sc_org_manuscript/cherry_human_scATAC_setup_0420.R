library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomeInfoDb)
library(ggplot2)
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
d53 <- create_obj("./d53/outs", "d53")
d74 <- create_obj("./d74/outs", "d74")
d78 <- create_obj("./d78/outs", "d78")
Hu5 <- create_obj("./Hu5/outs", "Hu5")
Hu7 <- create_obj("./Hu7/outs", "Hu7")
Hu8 <- create_obj("./Hu8/outs", "Hu8")

##graph QC and subset, then normalize etc
##d53
VlnPlot(object = d53, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d53.qc.pdf", width=20)
FragmentHistogram(object = d53, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d53.frag_by_ns.pdf", width=20)
TSSPlot(d53, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d53.tss_enrichment.pdf", width=20)
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
d53_subset <- subset(d53, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
d53
d53_subset
##Normalization and linear dimensional reduction
d53_subset <- RunTFIDF(d53_subset)
d53_subset <- FindTopFeatures(d53_subset, min.cutoff = 'q0')
d53_subset <- RunSVD(object = d53_subset, assay = 'ATAC',
                     reduction.key = 'LSI_', reduction.name = 'lsi')
##correlation between depth and dims (may need to remove dim1 from analysis)
DepthCor(d53_subset)
##Non-linear dimension reduction and clustering
d53_subset <- RunUMAP(object = d53_subset, reduction = 'lsi', dims = 2:30)
d53_subset <- FindNeighbors(object = d53_subset, reduction = 'lsi', dims = 2:30)
d53_subset <- FindClusters(object = d53_subset, verbose = FALSE, algorithm = 3)
DimPlot(object = d53_subset, label = TRUE) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d53.UMAP.pdf", width=20)
##save file
saveRDS(d53_subset, file = "./signac_analysis/d53_subset.rds")

##d74
VlnPlot(object = d74, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d74.qc.pdf", width=20)
FragmentHistogram(object = d74, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d74.frag_by_ns.pdf", width=20)
TSSPlot(d74, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d74.tss_enrichment.pdf", width=20)
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
DepthCor(d74_subset,n = 10)
##Non-linear dimension reduction and clustering
d74_subset <- RunUMAP(object = d74_subset, reduction = 'lsi', dims = 2:30)
d74_subset <- FindNeighbors(object = d74_subset, reduction = 'lsi', dims = 2:30)
d74_subset <- FindClusters(object = d74_subset, verbose = FALSE, algorithm = 3)
DimPlot(object = d74_subset, label = TRUE) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d74.UMAP.pdf", width=20)
##save file
saveRDS(d74_subset, file = "./signac_analysis/d74_subset.rds")

##d78
VlnPlot(object = d78, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d78.qc.pdf", width=20)
FragmentHistogram(object = d78, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d78.frag_by_ns.pdf", width=20)
TSSPlot(d78, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d78.tss_enrichment.pdf", width=20)
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
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.d78.UMAP.pdf", width=20)
##save file
saveRDS(d78_subset, file = "./signac_analysis/d78_subset.rds")

##Hu5
VlnPlot(object = Hu5, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu5.qc.pdf", width=20)
FragmentHistogram(object = Hu5, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu5.frag_by_ns.pdf", width=20)
TSSPlot(Hu5, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu5.tss_enrichment.pdf", width=20)
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
Hu5_subset <- subset(Hu5, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
Hu5
Hu5_subset
##Normalization and linear dimensional reduction
Hu5_subset <- RunTFIDF(Hu5_subset)
Hu5_subset <- FindTopFeatures(Hu5_subset, min.cutoff = 'q0')
Hu5_subset <- RunSVD(object = Hu5_subset, assay = 'ATAC',
                     reduction.key = 'LSI_', reduction.name = 'lsi')
##correlation between depth and dims (may need to remove dim1 from analysis)
DepthCor(Hu5_subset)
##Non-linear dimension reduction and clustering
Hu5_subset <- RunUMAP(object = Hu5_subset, reduction = 'lsi', dims = 2:30)
Hu5_subset <- FindNeighbors(object = Hu5_subset, reduction = 'lsi', dims = 2:30)
Hu5_subset <- FindClusters(object = Hu5_subset, verbose = FALSE, algorithm = 3)
DimPlot(object = Hu5_subset, label = TRUE) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu5.UMAP.pdf", width=20)
##save file
saveRDS(Hu5_subset, file = "./signac_analysis/Hu5_subset.rds")

##Hu7
VlnPlot(object = Hu7, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu7.qc.pdf", width=20)
FragmentHistogram(object = Hu7, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu7.frag_by_ns.pdf", width=20)
TSSPlot(Hu7, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu7.tss_enrichment.pdf", width=20)
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
Hu7_subset <- subset(Hu7, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
Hu7
Hu7_subset
##Normalization and linear dimensional reduction
Hu7_subset <- RunTFIDF(Hu7_subset)
Hu7_subset <- FindTopFeatures(Hu7_subset, min.cutoff = 'q0')
Hu7_subset <- RunSVD(object = Hu7_subset, assay = 'ATAC',
                     reduction.key = 'LSI_', reduction.name = 'lsi')
##correlation between depth and dims (may need to remove dim1 from analysis)
DepthCor(Hu7_subset)
##Non-linear dimension reduction and clustering
Hu7_subset <- RunUMAP(object = Hu7_subset, reduction = 'lsi', dims = 2:30)
Hu7_subset <- FindNeighbors(object = Hu7_subset, reduction = 'lsi', dims = 2:30)
Hu7_subset <- FindClusters(object = Hu7_subset, verbose = FALSE, algorithm = 3)
DimPlot(object = Hu7_subset, label = TRUE) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu7.UMAP.pdf", width=20)
##save file
saveRDS(Hu7_subset, file = "./signac_analysis/Hu7_subset.rds")

##Hu8
VlnPlot(object = Hu8, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu8.qc.pdf", width=20)
FragmentHistogram(object = Hu8, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu8.frag_by_ns.pdf", width=20)
TSSPlot(Hu8, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu8.tss_enrichment.pdf", width=20)
##filter using qc metrics
#pbmc <- subset(pbmc, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
Hu8_subset <- subset(Hu8, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
Hu8
Hu8_subset
##Normalization and linear dimensional reduction
Hu8_subset <- RunTFIDF(Hu8_subset)
Hu8_subset <- FindTopFeatures(Hu8_subset, min.cutoff = 'q0')
Hu8_subset <- RunSVD(object = Hu8_subset, assay = 'ATAC',
                     reduction.key = 'LSI_', reduction.name = 'lsi')
##correlation between depth and dims (may need to remove dim1 from analysis)
DepthCor(Hu8_subset)
##Non-linear dimension reduction and clustering
Hu8_subset <- RunUMAP(object = Hu8_subset, reduction = 'lsi', dims = 2:30)
Hu8_subset <- FindNeighbors(object = Hu8_subset, reduction = 'lsi', dims = 2:30)
Hu8_subset <- FindClusters(object = Hu8_subset, verbose = FALSE, algorithm = 3)
DimPlot(object = Hu8_subset, label = TRUE) + NoLegend()
dev.copy2pdf(file="./signac_analysis/human_scatac_0420.Hu8.UMAP.pdf", width=20)
##save file
saveRDS(Hu8_subset, file = "./signac_analysis/Hu8_subset.rds")




