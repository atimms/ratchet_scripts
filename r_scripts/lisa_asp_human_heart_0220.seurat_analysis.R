library(Seurat)
library(ggplot2)
library(reticulate)
library(cowplot)
library(devtools)
install_github("ririzarr/rafalib")
library(rafalib)
setwd('/data/atimms/lisa_asp_human_heart_0220')


asp.data <- read.table("all_cells_count_matrix_filtered.tsv",header=T,row.names = 1)
asp <- CreateSeuratObject(counts = asp.data, project = "asp", min.cells = 3, min.features = 200)


# store mitochondrial percentage in object meta data (mt- is for mouse)
asp[["percent.mt"]] <- PercentageFeatureSet(asp, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(asp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.copy2pdf(file="asp_serat_0220.qc.pdf", width = 10)

# filter the data
asp <- subset(asp, subset = nFeature_RNA > 400 & nFeature_RNA < 7500 & percent.mt < 10)

# run sctransform
asp <- SCTransform(asp, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
asp <- RunPCA(asp, verbose = FALSE)
asp <- RunUMAP(asp, dims = 1:30, verbose = FALSE)
asp <- FindNeighbors(asp, dims = 1:30, verbose = FALSE)
asp <- FindClusters(asp, verbose = FALSE, resolution = 0.4)
DimPlot(asp, label = TRUE) 
dev.copy2pdf(file="asp_serat_0220.umap_cluster.pdf")

# find markers for every cluster compared to all remaining cells, report only the positive ones
asp.integrated.markers <- FindAllMarkers(asp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(asp.integrated.markers, file="asp_serat_0220.cluster_markers.csv")

##cluster expression
cluster.averages <- AverageExpression(asp)
head(cluster.averages[["RNA"]][, 1:5])
write.csv(cluster.averages, file="asp_serat_0220.cluster_gene_expression.csv")

##look at clusters etc....
cl_data <- read.table('asp_cluster_data.txt', header=T, row.names=1)
head(cl_data)
##convert data to percentages
cl_perc = cl_data/rowSums(cl_data[1:18])
head(cl_data[1:18])
head(cl_perc)
##Hierarchical clustering
sampleDists <- dist( cl_perc )
sampleDistMatrix <- as.matrix( sampleDists )
hc <- hclust(sampleDists)
hc
##make graphs
myplclust(hc, labels=rownames(cl_data), lab.col=as.fumeric(as.vector(cl_data$in_genelist)), cex=0.5)
abline(h=0.3)
dev.copy2pdf(file='asp_cluster_data_0220.hc.perc.pdf', width = 20, height = 10)
hclusters <- cutree(hc, h=0.3)
cluster_table = table(true=cl_data$in_genelist, cluster=hclusters)
cluster_table
write.table(cluster_table, file="asp_cluster_data_0220.hc.perc.counts.txt")
##add cluster id to data
data_with_cid <- cbind(cl_data, clusterID=hclusters)
write.csv(data_with_cid, file="asp_cluster_data_0220.hc.cluster_id.csv")
data_pc_with_cid <- cbind(cl_perc, clusterID=hclusters)
write.csv(data_pc_with_cid, file="asp_cluster_data_0220.hc.cluster_id.pc.csv")
##heatmap on cluster.... not done...
c11_data <- read.table('hpf_24_set1_cl5.txt', header=T, row.names=1)
mat <- c11_data[1:72]
newdf <- c(c11_data$gene, c11_data$in_genelist)
pheatmap(mat, fontsize_row=2, fontsize_col=10, cluster_cols = F)
dev.copy2pdf(file='neural_tube_cluster_data_1_cluster1.pdf', width = 7, height = 5)