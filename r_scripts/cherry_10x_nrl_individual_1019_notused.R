library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
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
nrl_0319.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_data/raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319 <- CreateSeuratObject(counts = nrl_0319.data, project = "nrl_0319", min.cells = 3, min.features = 200)
nrl_0319
# store mitochondrial percentage in object meta data
nrl_0319 <- PercentageFeatureSet(nrl_0319, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319 <- SCTransform(nrl_0319, vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
nrl_0319 <- RunPCA(nrl_0319, verbose = FALSE)
nrl_0319 <- RunUMAP(nrl_0319, dims = 1:30, verbose = FALSE)
nrl_0319 <- FindNeighbors(nrl_0319, dims = 1:30, verbose = FALSE)
nrl_0319 <- FindClusters(nrl_0319, verbose = FALSE)
DimPlot(nrl_0319, label = TRUE) + NoLegend()
dev.copy2pdf(file="nrl_10x_aggr_0319.umap.pdf", width = 10)

# Load the 10x dataset for 1019
nrl_1019.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_1019_data/raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_1019 <- CreateSeuratObject(counts = nrl_1019.data, project = "nrl_1019", min.cells = 3, min.features = 200)
nrl_1019
# store mitochondrial percentage in object meta data
nrl_1019 <- PercentageFeatureSet(nrl_1019, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_1019 <- SCTransform(nrl_1019, vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
nrl_1019 <- RunPCA(nrl_1019, verbose = FALSE)
nrl_1019 <- RunUMAP(nrl_1019, dims = 1:30, verbose = FALSE)
nrl_1019 <- FindNeighbors(nrl_1019, dims = 1:30, verbose = FALSE)
nrl_1019 <- FindClusters(nrl_1019, verbose = FALSE)
DimPlot(nrl_1019, label = TRUE) + NoLegend()
dev.copy2pdf(file="nrl_10x_aggr_1019.umap.pdf", width = 10)

##lets get the data sample by sample
nrl_list = c("nrl_0319_coding_raw_feature_bc_matrix", "nrl_0319_enhancer1_raw_feature_bc_matrix", "nrl_0319_enhancer2_raw_feature_bc_matrix", 
             "nrl_0319_enhancer3_raw_feature_bc_matrix", "nrl_0319_introniccontrol_raw_feature_bc_matrix", "nrl_0319_noelectroporation_raw_feature_bc_matrix", 
             "nrl_0319_nontargetingcontrol_raw_feature_bc_matrix", "nrl_0319_promoter_raw_feature_bc_matrix", "nrl_1019_coding_raw_feature_bc_matrix", 
             "nrl_1019_E2_raw_feature_bc_matrix", "nrl_1019_E3_raw_feature_bc_matrix", "nrl_1019_Intron_raw_feature_bc_matrix", "nrl_1019_NoEP_raw_feature_bc_matrix", 
             "nrl_1019_NTC_raw_feature_bc_matrix", "nrl_1019_Promoter_raw_feature_bc_matrix")

##nrl_0319_coding
nrl_0319_coding.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_coding_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_coding <- CreateSeuratObject(counts = nrl_0319_coding.data, project = "nrl_0319_coding", min.cells = 3, min.features = 200)
nrl_0319_coding
##add info
nrl_0319_coding$treatment <- "coding"
# store mitochondrial percentage in object meta data
nrl_0319_coding <- PercentageFeatureSet(nrl_0319_coding, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_coding <- SCTransform(nrl_0319_coding, vars.to.regress = "percent.mt", verbose = FALSE)

##nrl_0319_enhancer1
nrl_0319_enhancer1.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_enhancer1_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_enhancer1 <- CreateSeuratObject(counts = nrl_0319_enhancer1.data, project = "nrl_0319_enhancer1", min.cells = 3, min.features = 200)
nrl_0319_enhancer1
##add info
nrl_0319_enhancer1$treatment <- "enhancer1"
# store mitochondrial percentage in object meta data
nrl_0319_enhancer1 <- PercentageFeatureSet(nrl_0319_enhancer1, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_enhancer1 <- SCTransform(nrl_0319_enhancer1, vars.to.regress = "percent.mt", verbose = FALSE)

##nrl_0319_enhancer2
nrl_0319_enhancer2.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_enhancer2_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_enhancer2 <- CreateSeuratObject(counts = nrl_0319_enhancer2.data, project = "nrl_0319_enhancer2", min.cells = 3, min.features = 200)
nrl_0319_enhancer2
##add info
nrl_0319_enhancer2$treatment <- "enhancer2"
# store mitochondrial percentage in object meta data
nrl_0319_enhancer2 <- PercentageFeatureSet(nrl_0319_enhancer2, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_enhancer2 <- SCTransform(nrl_0319_enhancer2, vars.to.regress = "percent.mt", verbose = FALSE)

##nrl_0319_enhancer3
nrl_0319_enhancer3.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_enhancer3_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_enhancer3 <- CreateSeuratObject(counts = nrl_0319_enhancer3.data, project = "nrl_0319_enhancer3", min.cells = 3, min.features = 200)
nrl_0319_enhancer3
##add info
nrl_0319_enhancer3$treatment <- "enhancer3"
# store mitochondrial percentage in object meta data
nrl_0319_enhancer3 <- PercentageFeatureSet(nrl_0319_enhancer3, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_enhancer3 <- SCTransform(nrl_0319_enhancer3, vars.to.regress = "percent.mt", verbose = FALSE)

##nrl_0319_introniccontrol
nrl_0319_introniccontrol.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_introniccontrol_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_introniccontrol <- CreateSeuratObject(counts = nrl_0319_introniccontrol.data, project = "nrl_0319_introniccontrol", min.cells = 3, min.features = 200)
nrl_0319_introniccontrol
nrl_0319_introniccontrol$treatment <- "introniccontrol"
# store mitochondrial percentage in object meta data
nrl_0319_introniccontrol <- PercentageFeatureSet(nrl_0319_introniccontrol, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_introniccontrol <- SCTransform(nrl_0319_introniccontrol, vars.to.regress = "percent.mt", verbose = FALSE)

##nrl_0319_noelectroporation
nrl_0319_noelectroporation.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_noelectroporation_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_noelectroporation$treatment <- "noelectroporation"
nrl_0319_noelectroporation <- CreateSeuratObject(counts = nrl_0319_noelectroporation.data, project = "nrl_0319_noelectroporation", min.cells = 3, min.features = 200)
nrl_0319_noelectroporation
##add info
nrl_0319_noelectroporation$treatment <- "noelectroporation"
# store mitochondrial percentage in object meta data
nrl_0319_noelectroporation <- PercentageFeatureSet(nrl_0319_noelectroporation, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_noelectroporation <- SCTransform(nrl_0319_noelectroporation, vars.to.regress = "percent.mt", verbose = FALSE)

##nrl_0319_nontargetingcontrol
nrl_0319_nontargetingcontrol.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_nontargetingcontrol_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_nontargetingcontrol <- CreateSeuratObject(counts = nrl_0319_nontargetingcontrol.data, project = "nrl_0319_nontargetingcontrol", min.cells = 3, min.features = 200)
nrl_0319_nontargetingcontrol
##add info
nrl_0319_nontargetingcontrol$treatment <- "nontargetingcontrol"
# store mitochondrial percentage in object meta data
nrl_0319_nontargetingcontrol <- PercentageFeatureSet(nrl_0319_nontargetingcontrol, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_nontargetingcontrol <- SCTransform(nrl_0319_nontargetingcontrol, vars.to.regress = "percent.mt", verbose = FALSE)

##nrl_0319_promoter
nrl_0319_promoter.data <- Read10X(data.dir = "/data/atimms/cherry_10x_nrl_1019/nrl_0319_promoter_raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
nrl_0319_promoter <- CreateSeuratObject(counts = nrl_0319_promoter.data, project = "nrl_0319_promoter", min.cells = 3, min.features = 200)
nrl_0319_promoter
##add info
nrl_0319_promoter$treatment <- "promoter"
# store mitochondrial percentage in object meta data
nrl_0319_promoter <- PercentageFeatureSet(nrl_0319_promoter, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
nrl_0319_promoter <- SCTransform(nrl_0319_promoter, vars.to.regress = "percent.mt", verbose = FALSE)

#select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
nrl.list <- list(nrl_0319_coding, nrl_0319_enhancer1, nrl_0319_enhancer2, nrl_0319_enhancer3, nrl_0319_introniccontrol, nrl_0319_noelectroporation, 
                 nrl_0319_nontargetingcontrol, nrl_0319_promoter)
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
dev.copy2pdf(file="nrl_10x_0319_integrated.umap_treatment.pdf", width = 10)
p2 <- DimPlot(nrl.integrated, reduction = "umap", label = TRUE)
p2
dev.copy2pdf(file="nrl_10x_0319_integrated.umap_cluster.pdf", width = 10)
plot_grid(p1, p2)
DimPlot(nrl.integrated, reduction = "umap", split.by = "treatment")
dev.copy2pdf(file="nrl_10x_0319_integrated.umap_cluster_treatment.pdf", width = 10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
retina.integrated.markers <- FindAllMarkers(retina.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(retina.integrated.markers, file="cherry_human_mouse_seurat_0919.cluster_markers.csv")
