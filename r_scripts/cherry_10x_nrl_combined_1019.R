library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(stringr)
setwd('/data/atimms/cherry_10x_nrl_1019')

##need this to use umap
#reticulate::py_install(packages ='umap-learn')

##set for integrating the datasets
options(future.globals.maxSize = 4000 * 1024^2)


##taking from:
## https://satijalab.org/seurat/v3.1/immune_alignment.html
## https://satijalab.org/seurat/v3.0/integration.html
## https://hbctraining.github.io/scRNA-seq/lessons/02_SC_quality_control-setup.html

# Load the 10x dataset for 0319
nrl_0319.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_data/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319 <- CreateSeuratObject(counts = nrl_0319.data, project = "nrl_0319", min.cells = 3, min.features = 200, names.delim = "-",names.field = 2)
nrl_0319
head(nrl_0319@meta.data)
# store mitochondrial percentage in object meta data (mt- is for mouse)
nrl_0319 <- PercentageFeatureSet(nrl_0319, pattern = "^mt-", col.name = "percent.mt")
# Compute percent mito ratio
nrl_0319$mitoRatio <- PercentageFeatureSet(object = nrl_0319, pattern = "^mt-")
nrl_0319$mitoRatio <- nrl_0319@meta.data$mitoRatio / 100
head(nrl_0319@meta.data)
#nrl_0319 <- RenameIdents(object = nrl_0319, '1' = 'coding', '2' = 'e1', '3' = 'e2', '4' = 'e3', '5' = 'intron', '6' = 'noep', '7'= 'ntc', '8' = 'promoter')
# Add number of genes per UMI for each cell to metadata
nrl_0319$log10GenesPerUMI <- log10(nrl_0319$nFeature_RNA) / log10(nrl_0319$nCount_RNA)

##change metadata
# Create metadata dataframe
metadata <- nrl_0319@meta.data
# Rename columns to make more sense
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
tail(metadata)
# Create treatment column and change ids
metadata$treatment <- NA
metadata$treatment[which(str_detect(metadata$orig.ident, "1"))] <- "coding"
metadata$treatment[which(str_detect(metadata$orig.ident, "2"))] <- "e1"
metadata$treatment[which(str_detect(metadata$orig.ident, "3"))] <- "e2"
metadata$treatment[which(str_detect(metadata$orig.ident, "4"))] <- "e3"
metadata$treatment[which(str_detect(metadata$orig.ident, "5"))] <- "intron"
metadata$treatment[which(str_detect(metadata$orig.ident, "6"))] <- "noep"
metadata$treatment[which(str_detect(metadata$orig.ident, "7"))] <- "ntc"
metadata$treatment[which(str_detect(metadata$orig.ident, "8"))] <- "promoter"
# Add metadata back to Seurat object
nrl_0319@meta.data <- metadata
# Visualize the number of cell counts per cell
metadata %>% 
  ggplot(aes(x=treatment, fill=treatment)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave("nrl_10x_0319.qc_ncells.pdf")
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=treatment, x=nUMI, fill= treatment)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("nrl_10x_0319.qc_cell_density.pdf")
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=treatment, x=nGene, fill= treatment)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
ggsave("nrl_10x_0319.qc_cell_density_umi.pdf")
# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=treatment, y=log10(nGene), fill=treatment)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("nrl_10x_0319.qc_gene_boxplot.pdf")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~treatment)
ggsave("nrl_10x_0319.qc_combined.pdf")
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=treatment, x=mitoRatio, fill=treatment)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave("nrl_10x_0319.qc_mt.pdf")
# Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = treatment, fill=treatment)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave("nrl_10x_0319.qc_novelty.pdf")
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_nrl_0319 <- subset(x = nrl_0319, 
                          subset= (nUMI >= 250) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))
##keep only genes which are expressed in 10 or more cells.....
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_nrl_0319, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0L
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- rowSums(as.matrix(nonzero)) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Create a new Seurat object
clean_nrl_0319 <- CreateSeuratObject(filtered_counts, meta.data = filtered_nrl_0319@meta.data)

# run sctransform
clean_nrl_0319 <- SCTransform(clean_nrl_0319, vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
clean_nrl_0319 <- RunPCA(clean_nrl_0319, verbose = FALSE)
clean_nrl_0319 <- RunUMAP(clean_nrl_0319, dims = 1:30, verbose = FALSE)
clean_nrl_0319 <- FindNeighbors(clean_nrl_0319, dims = 1:30, verbose = FALSE)
clean_nrl_0319 <- FindClusters(clean_nrl_0319, verbose = FALSE)
DimPlot(clean_nrl_0319, label = TRUE) + NoLegend()
dev.copy2pdf(file="nrl_10x_0319_clean.umap.pdf", width = 10)
DimPlot(clean_nrl_0319, reduction = "umap", group.by = "treatment")
dev.copy2pdf(file="nrl_10x_0319_clean.umap_treatment.pdf", width = 10)


# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(clean_nrl_0319) <- "RNA"
#Idents(clean_nrl_0319 = clean_nrl_0319) <- clean_nrl_0319@meta.data$seurat_clusters
nrl_0319.markers <- FindAllMarkers(clean_nrl_0319, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(nrl_0319.markers, file="nrl_10x_0319.cluster_markers.csv")

##write metadata
write.csv(clean_nrl_0319@meta.data, file="nrl_10x_0319.metadata.csv")



# Load the 10x dataset for 1019
nrl_1019.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_1019_data/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_1019 <- CreateSeuratObject(counts = nrl_1019.data, project = "nrl_1019", min.cells = 3, min.features = 200, names.delim = "-",names.field = 2)
nrl_1019
head(nrl_1019@meta.data)
# store mitochondrial percentage in object meta data (mt- is for mouse)
nrl_1019 <- PercentageFeatureSet(nrl_1019, pattern = "^mt-", col.name = "percent.mt")
# Compute percent mito ratio
nrl_1019$mitoRatio <- PercentageFeatureSet(object = nrl_1019, pattern = "^mt-")
nrl_1019$mitoRatio <- nrl_1019@meta.data$mitoRatio / 100
head(nrl_1019@meta.data)
# Add number of genes per UMI for each cell to metadata
nrl_1019$log10GenesPerUMI <- log10(nrl_1019$nFeature_RNA) / log10(nrl_1019$nCount_RNA)

##change metadata
# Create metadata dataframe
metadata <- nrl_1019@meta.data
# Rename columns to make more sense
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
head(metadata)
# Create treatment column and change ids
metadata$treatment <- NA
metadata$treatment[which(str_detect(metadata$orig.ident, "1"))] <- "coding"
metadata$treatment[which(str_detect(metadata$orig.ident, "2"))] <- "e2"
metadata$treatment[which(str_detect(metadata$orig.ident, "3"))] <- "e3"
metadata$treatment[which(str_detect(metadata$orig.ident, "4"))] <- "intron"
metadata$treatment[which(str_detect(metadata$orig.ident, "5"))] <- "noep"
metadata$treatment[which(str_detect(metadata$orig.ident, "6"))] <- "ntc"
metadata$treatment[which(str_detect(metadata$orig.ident, "7"))] <- "promoter"
# Add metadata back to Seurat object
nrl_1019@meta.data <- metadata
# Visualize the number of cell counts per cell
metadata %>% 
  ggplot(aes(x=treatment, fill=treatment)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave("nrl_10x_1019.qc_ncells.pdf")
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=treatment, x=nUMI, fill= treatment)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("nrl_10x_1019.qc_cell_density.pdf")
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=treatment, x=nGene, fill= treatment)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
ggsave("nrl_10x_1019.qc_cell_density_umi.pdf")
# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=treatment, y=log10(nGene), fill=treatment)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("nrl_10x_1019.qc_gene_boxplot.pdf")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~treatment)
ggsave("nrl_10x_1019.qc_combined.pdf")
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=treatment, x=mitoRatio, fill=treatment)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave("nrl_10x_1019.qc_mt.pdf")
# Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = treatment, fill=treatment)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave("nrl_10x_1019.qc_novelty.pdf")
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_nrl_1019 <- subset(x = nrl_1019, 
                            subset= (nUMI >= 250) & 
                              (nGene >= 250) & 
                              (log10GenesPerUMI > 0.80) & 
                              (mitoRatio < 0.20))
##keep only genes which are expressed in 10 or more cells.....
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_nrl_1019, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0L
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- rowSums(as.matrix(nonzero)) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Create a new Seurat object
clean_nrl_1019 <- CreateSeuratObject(filtered_counts, meta.data = filtered_nrl_1019@meta.data)

# run sctransform
clean_nrl_1019 <- SCTransform(clean_nrl_1019, vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
clean_nrl_1019 <- RunPCA(clean_nrl_1019, verbose = FALSE)
clean_nrl_1019 <- RunUMAP(clean_nrl_1019, dims = 1:30, verbose = FALSE)
clean_nrl_1019 <- FindNeighbors(clean_nrl_1019, dims = 1:30, verbose = FALSE)
clean_nrl_1019 <- FindClusters(clean_nrl_1019, verbose = FALSE)
DimPlot(clean_nrl_1019, label = TRUE) + NoLegend()
dev.copy2pdf(file="nrl_10x_1019_clean.umap.pdf", width = 10)
DimPlot(clean_nrl_1019, reduction = "umap", group.by = "treatment")
dev.copy2pdf(file="nrl_10x_1019_clean.umap_treatment.pdf", width = 10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(clean_nrl_1019) <- "RNA"
#Idents(clean_nrl_1019 = clean_nrl_1019) <- clean_nrl_1019@meta.data$seurat_clusters
nrl_1019.markers <- FindAllMarkers(clean_nrl_1019, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(nrl_1019.markers, file="nrl_10x_1019.cluster_markers.csv")

##write metadata
write.csv(clean_nrl_1019@meta.data, file="nrl_10x_1019.metadata.csv")


#select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
nrl.list <- list(clean_nrl_0319, clean_nrl_1019)
nrl.features <- SelectIntegrationFeatures(object.list = nrl.list, nfeatures = 3000)
nrl.list <- PrepSCTIntegration(object.list = nrl.list, anchor.features = nrl.features, 
                               verbose = FALSE)
##identify anchors and integrate the datasets.
nrl.anchors <- FindIntegrationAnchors(object.list = nrl.list, normalization.method = "SCT", 
                                      anchor.features = nrl.features, verbose = FALSE)
nrl.integrated <- IntegrateData(anchorset = nrl.anchors, normalization.method = "SCT", 
                                verbose = FALSE)
##Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
nrl.integrated <- RunPCA(nrl.integrated, verbose = FALSE)
nrl.integrated <- RunUMAP(nrl.integrated, dims = 1:30)
nrl.integrated <- FindNeighbors(nrl.integrated, reduction = "pca", dims = 1:30)
##which is best resolution, 0.5 is default
nrl.integrated <- FindClusters(nrl.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(nrl.integrated, reduction = "umap", group.by = "treatment")
p1
dev.copy2pdf(file="nrl_10x_integrated.umap_treatment.pdf", width = 10)
p2 <- DimPlot(nrl.integrated, reduction = "umap", label = TRUE)
p2
dev.copy2pdf(file="nrl_10x_integrated.umap_cluster.pdf", width = 10)
head(nrl.integrated@meta.data)
plot_grid(p1, p2)
DimPlot(nrl.integrated, reduction = "umap", split.by = "treatment")
dev.copy2pdf(file="nrl_10x_integrated.umap_cluster_treatment.pdf", width = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(nrl.integrated) <- "RNA"
nrl.integrated.markers <- FindAllMarkers(nrl.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(nrl.integrated.markers, file="nrl_10x_integrated.cluster_markers.csv")

##write metadata
write.csv(nrl.integrated@meta.data, file="nrl_10x_integrated.metadata.csv")


