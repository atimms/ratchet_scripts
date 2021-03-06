# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA')

# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)

#Import 10x data and convert each dataset to Seurat object (take from individual sample outputs, not cellranger's aggr output). Define sample with project= "samplename"
#embyonic data
d53.data = Read10X(data.dir = './d53/outs/filtered_feature_bc_matrix/')
d53= CreateSeuratObject(counts = d53.data, project = "d53", min.cells = 3, min.features = 200)
d74.data = Read10X(data.dir = './d74/outs/filtered_feature_bc_matrix/')
d74= CreateSeuratObject(counts = d74.data, project = "d74", min.cells = 3, min.features = 200)
d78.data = Read10X(data.dir = './d78/outs/filtered_feature_bc_matrix/')
d78= CreateSeuratObject(counts = d78.data, project = "d78", min.cells = 3, min.features = 200)

#adult Data
Hu37.data = Read10X(data.dir = './H37/outs/filtered_feature_bc_matrix/')
Hu37= CreateSeuratObject(counts = Hu37.data, project = "Hu37", min.cells = 3, min.features = 200)
Hu5.data = Read10X(data.dir = './Hu5/outs/filtered_feature_bc_matrix/')
Hu5= CreateSeuratObject(counts = Hu5.data, project = "Hu5", min.cells = 3, min.features = 200)
Hu7.data = Read10X(data.dir = './Hu7/outs/filtered_feature_bc_matrix/')
Hu7= CreateSeuratObject(counts = Hu7.data, project = "Hu7", min.cells = 3, min.features = 200)

# Merge into one single Seurat object
human=merge(d53, y=c(d74,d78,Hu37,Hu5,Hu7))
#Validate the merge by checking number of cells per group
table(human$orig.ident)

#Store mitochondrial percentage in the Seurat object metadata
human[["percent.mt"]] <- PercentageFeatureSet(human, pattern = "^MT-")

#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.qc.pdf", width=20)

#Add sample and condition information explicitly into the metadata (as oppsoed to storing in 'orig.ident') for future downstream analysis
#Add sample info (simply copying 'orig.ident')
human$sample <- plyr::mapvalues(
  x = human$orig.ident, 
  from = c('d53', 'd74', 'd78', 'Hu37', 'Hu5', 'Hu7'), 
  to = c('d53', 'd74', 'd78', 'Hu37', 'Hu5', 'Hu7')
)
#Add phenotype info (i.e. WT vs KO; allows for merging of different bioreps)
human$time <- plyr::mapvalues(
  x = human$orig.ident, 
  from = c('d53', 'd74', 'd78', 'Hu37', 'Hu5', 'Hu7'),  
  to = c('e50s', 'e70s', 'e70s', 'adult', 'adult', 'adult')
)

#Validate new metadata columns by checking that number of cells per sample/phenotype adds up
table(human$sample)
table(human$time)

#Run same QC metrics by new metadata columns to ensure it is the same as original QC metrics
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by='sample', ncol = 3, pt.size=0.1)
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by='time', ncol = 3, pt.size=0.1)

#Filter the data
human <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human_precheck = human
human_precheck=NormalizeData(human_precheck)
human_precheck <- FindVariableFeatures(human_precheck, selection.method = "vst", nfeatures = 2000)
human_precheck=ScaleData(human_precheck)
human_precheck = RunPCA(human_precheck)
DimPlot(human_precheck, reduction='pca', group.by='sample')
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.pca.pdf", width=20)
ElbowPlot(human_precheck, ndims=30)
dev.copy2pdf(file="./seurat_analysis/human_scrnaseq_0420.elbow_plot.pdf", width=20)

#If samples do not cluster together in PC space (which is the case here), then there is no need to run harmony (and likely no need to use the integrated method either); the merged analysis should do.
#Determine the number of dimensions to use in downstream analyses based on the point at which the Elbow Plot becomes flat (ok to be conservative)

#Save merged and filtered dataset as an R object to be able to re-load it for various future analyses without needing to perform the previous computations
saveRDS(human, file = "./seurat_analysis/human.rds")
saveRDS(human_precheck, file = "./seurat_analysis/human_precheck.rds")
##and load if needed
human <- readRDS(file = "./seurat_analysis/human.rds")
human_precheck <- readRDS(file = "./seurat_analysis/human_precheck.rds")


