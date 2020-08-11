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
Hu5 <- create_obj("./Hu5/outs", "Hu5")
Hu7 <- create_obj("./Hu7/outs", "Hu7")
Hu8 <- create_obj("./Hu8/outs", "Hu8")

##graph QC and subset, then normalize etc
##Hu5
VlnPlot(object = Hu5, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu5.qc.pdf", width=20)
FragmentHistogram(object = Hu5, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu5.frag_by_ns.pdf", width=20)
TSSPlot(Hu5, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu5.tss_enrichment.pdf", width=20)
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
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.UMAP.pdf", width=20)

##Hu7
VlnPlot(object = Hu7, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu7.qc.pdf", width=20)
FragmentHistogram(object = Hu7, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu7.frag_by_ns.pdf", width=20)
TSSPlot(Hu7, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu7.tss_enrichment.pdf", width=20)
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
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.UMAP.pdf", width=20)

##Hu8
VlnPlot(object = Hu8, features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),pt.size = 0.1,ncol = 4) + NoLegend()
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu8.qc.pdf", width=20)
FragmentHistogram(object = Hu8, group.by = 'nucleosome_group')
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu8.frag_by_ns.pdf", width=20)
TSSPlot(Hu8, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.Hu8.tss_enrichment.pdf", width=20)
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
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0620.UMAP.pdf", width=20)

##merge data
## create a common peak set
combined.peaks <- UnifyPeaks(object.list = list(Hu5_subset,Hu7_subset,Hu8_subset), mode = "reduce")
combined.peaks

## quantify peaks in each dataset, and add peaks assay to object
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
Hu5_subset$dataset <- 'Hu5'
Hu7_subset$dataset <- 'Hu7'
Hu8_subset$dataset <- 'Hu8'

# make sure to change to the assay containing common peaks
DefaultAssay(Hu5_subset) <- "peaks" 
DefaultAssay(Hu7_subset) <- "peaks" 
DefaultAssay(Hu8_subset) <- "peaks"

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(x = Hu5_subset, y = list(Hu7_subset,Hu8_subset), add.cell.ids = c("Hu5", "Hu7", "Hu8"))

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
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.merged.res_0.4.UMAP_by_cluster.pdf", width=20)
DimPlot(combined, split.by = 'dataset', pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.merged.res_0.4.UMAP_split_by_sample.pdf", width=20)

##merged fragments using unix commands, so add to merged object
combined <- SetFragments(combined, "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments_adult/fragments.tsv.gz")

# save individual samples with peaks assay and combined object
saveRDS(Hu5_subset, file = "./signac_analysis/adult/Hu5_combPeaks.rds")
saveRDS(Hu7_subset, file = "./signac_analysis/adult/Hu7_combPeaks.rds")
saveRDS(Hu8_subset, file = "./signac_analysis/adult/Hu8_combPeaks.rds")
saveRDS(combined, file = "./signac_analysis/adult/combined.rds")
saveRDS(combined.peaks, file = "./signac_analysis/adult/combined_peaks.rds")

##load file if needed
#combined <- readRDS(file = "./signac_analysis/adult/combined.rds")

##find markers
# switch back to working with peaks instead of gene activities
DefaultAssay(combined) <- 'peaks'
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(combined.markers, file='./signac_analysis/adult/human_scatac_adult_0520.merged.markers.csv')
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
write.csv(closest_genes, file='./signac_analysis/adult/human_scatac_adult_0520.merged.marker_gene_ids.csv')


##make gene matrix and add gene names instead of peaks

# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
fragment.path = "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments_adult/fragments.tsv.gz"
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
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.merged.known_marker.dotplot.res0.4.pdf", width = 20)



##run harmony

## load data if needed - from merged data
combined <- readRDS(file = "./signac_analysis/adult/combined.rds")

hm.integrated <- RunHarmony(object = combined, group.by.vars = 'dataset',
                            reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)

# re-compute the UMAP using corrected LSI embeddings
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')

##find clusters for 2 different resolutions
hm.integrated <- FindNeighbors(hm.integrated, reduction = 'harmony', dims = 2:30)
##res 0.4
hm.integrated <- FindClusters(hm.integrated, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(object = hm.integrated, pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.harmony_res0.4.UMAP_by_cluster.pdf", width=20)
DimPlot(hm.integrated, split.by = 'dataset', pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.harmony_res0.4.UMAP_split_by_sample.pdf", width=20)

##save/load file
saveRDS(hm.integrated, file = "./signac_analysis/adult/harmony_integrated.rds")
#hm.integrated <- readRDS(file = "./signac_analysis/adult/harmony_integrated.rds")

##make gene matrix and add gene names instead of peaks

# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
fragment.path = "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/merge_fragments_adult/fragments.tsv.gz"
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
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.harmony.feature_plot1.res0.4.pdf", width=20)

##save/load file
saveRDS(hm.integrated, file = "./signac_analysis/adult/harmony_integrated.rds")
hm.integrated <- readRDS(file = "./signac_analysis/adult/harmony_integrated.rds")

##find markers
# switch back to working with peaks instead of gene activities
DefaultAssay(hm.integrated) <- 'peaks'
harmony.markers <- FindAllMarkers(hm.integrated, only.pos = TRUE, min.pct = 0.25, test.use = 'LR',latent.vars = 'peak_region_fragments')
write.csv(harmony.markers, file='./signac_analysis/adult/human_scatac_adult_0520.harmony.marker.res0.4.csv')
head(harmony.markers)

## get genes closest to marker genes, row names are weird in a few places
##temp way to lead results
#harmony.markers <- read.csv('./signac_analysis/adult/human_scatac_adult_0520.harmony.marker.res0.4.csv')
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
write.csv(closest_genes, file='./signac_analysis/adult/human_scatac_adult_0520.harmony.marker_gene_ids.res0.4.csv')

##make dotplots
DefaultAssay(hm.integrated) <- 'RNA'
#Dot plot - the size of the dot = % of cells and color represents the average expression
markers <- c('RCVRN', 'RHO', 'CRX', 'ARR3', 'GNAT2', 'VSX2', 'LHX4', 'TRPM1', 'GRM6', 'SLC1A3', 'RLBP1', 'PAX6', 'LHX1', 'ONECUT2', 'TFAP2B', 'GAD1', 'SLC6A9', 'RBPMS', 'NEFM', 'GFAP', 'CD74', 'P2RY12', 'BEST1', 'RPE65', 'SFRP2')
DotPlot(hm.integrated, features = markers) + RotatedAxis()
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.harmony.known_marker.dotplot.res0.4.pdf", width = 20)





##add cluster names and then extract cells and cluster ids
hm.integrated <- readRDS(file = "./signac_analysis/adult/harmony_integrated.rds")

#Add celltype info
new.cluster.ids <- c('Rods', 'Rods', 'Mullers', 'Rods_weird', 'Bipolars', 'Cones', 'Horizontals', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Amacrines', 'Amacrines', 'Amacrines', 'Amacrines', 'Ganglions')
names(new.cluster.ids) <- levels(hm.integrated)
hm.integrated <- RenameIdents(hm.integrated, new.cluster.ids)
DimPlot(hm.integrated, reduction = "umap", pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./signac_analysis/adult/human_scatac_adult_0520.harmony_res0.4.UMAP_celltypes.pdf", width=20)
WhichCells(object = hm.integrated, idents = c('Rods', 'Mullers', 'Bipolars', 'Cones', 'Horizontals', 'Amacrines', 'Ganglions'))
hm.integrated$celltype <- plyr::mapvalues(
  x = hm.integrated$seurat_clusters, 
  from = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'), 
  to = c('Rods', 'Rods', 'Mullers', 'Rods_weird', 'Bipolars', 'Cones', 'Horizontals', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Amacrines', 'Amacrines', 'Amacrines', 'Amacrines', 'Ganglions')
)
head(hm.integrated$celltype)
colnames(hm.integrated)
cell.info = hm.integrated$celltype
write.csv(cell.info, file='./signac_analysis/adult/human_scatac_adult_0520.harmony.res0.4.cellinfo.csv')

##get counts per sample/cell class
sample.celltype.info = table(hm.integrated$celltype, hm.integrated$dataset)
write.csv(sample.celltype.info, file='./signac_analysis/adult/human_scatac_adult_0520.harmony.res0.4.sample_cell_info.csv')
