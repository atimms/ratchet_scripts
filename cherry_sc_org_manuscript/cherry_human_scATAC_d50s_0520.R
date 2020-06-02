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


##graph QC and subset, then normalize etc
##d53
VlnPlot(object = d53, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/d50s/human_scatac_d50s_0520.d53.qc.pdf", width=20)
FragmentHistogram(object = d53, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/d50s/human_scatac_d50s_0520.d53.frag_by_ns.pdf", width=20)
TSSPlot(d53, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/d50s/human_scatac_d50s_0520.d53.tss_enrichment.pdf", width=20)
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
dev.copy2pdf(file="./signac_analysis/d50s/human_scatac_d50s_0520.UMAP.pdf", width=20)

# create granges object with TSS positions, different than previous
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
fragment.path = "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/d53/outs/fragments.tsv.gz"
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(d53_subset),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
d53_subset[['RNA']] <- CreateAssayObject(counts = gene.activities)
d53_subset <- NormalizeData(
  object = d53_subset,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(d53_subset$nCount_RNA)
)

##visualize canonical marker 
DefaultAssay(d53_subset) <- 'RNA'
FeaturePlot(object = d53_subset, features = c('SLC1A3', 'PROM1', 'RLBP1', 'ONECUT1', 'VSX1', 'SLC6A9', 'POU4F2', 'NGFR', 'C1QA', 'GJA4', 'BEST1', 'VIM', 'ATOH7', 'PRDM13', 'DLL3'), pt.size = 0.1, max.cutoff = 'q95', ncol = 4)
dev.copy2pdf(file="./signac_analysis/d50s/human_scatac_d50s_0520.feature_plot1.pdf", width=20)

##save/load file
saveRDS(d53_subset, file = "./signac_analysis/d50s/d53_subset.rds")
d53_subset <- readRDS(file = "./signac_analysis/d50s/d53_subset.rds")


##Find differentially accessible peaks between clusters and integrate RNA

