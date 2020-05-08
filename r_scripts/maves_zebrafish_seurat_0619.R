library(Seurat)
library(ggplot2)
library(reticulate)
library(cowplot)
use_python("~/miniconda3/bin/python")

workingDir = "/data/atimms/maves_zf_single_cell_0619";
setwd(workingDir);


#Import Cluster names
Cluster_names <- read.csv(file = "GSE112294_ClusterNames.csv")
#write.csv(Cluster_names, file = "/archive/dobyns_w/sc_data_1118/GSE112294_ClusterNames.csv")

#Read in hpf04 data
hpf04_clustID <- read.table(file = "GSM3067189_04hpf_clustID.txt")
head(hpf04_clustID)
hpf04 <- read.csv(file = "GSM3067189_04hpf.csv", row.names = 1)
hpf04 <- as.data.frame(hpf04)
hpf04[1:5,1:5]
#Seurat Pipeline
hpf04_seurat <- CreateSeuratObject(counts = hpf04, project = "hpf04_zebrafish", meta.data = hpf04_clustID)
hpf04_seurat <- NormalizeData(hpf04_seurat)
hpf04_seurat <- FindVariableFeatures(hpf04_seurat, selection.method = "vst", nfeatures = 2000)
hpf04_seurat <- ScaleData(hpf04_seurat)
hpf04_seurat <- RunPCA(hpf04_seurat, npcs = 50)
ElbowPlot(hpf04_seurat, ndims = 50)
hpf04_seurat <- FindNeighbors(hpf04_seurat, dims = 1:20)
hpf04_seurat <- FindClusters(hpf04_seurat, resolution = 1)
hpf04_seurat <- RunTSNE(hpf04_seurat, dims = 1:20)
hpf04_seurat <- RunUMAP(hpf04_seurat, reduction = "pca", dims = 1:20)
DimPlot(hpf04_seurat, reduction = "tsne", label = T) + NoLegend()
DimPlot(hpf04_seurat, reduction = "umap", label = T) + NoLegend()
#Add clusterID
hpf04_seurat$clustID <- hpf04_clustID$V1
DimPlot(hpf04_seurat, reduction = "tsne", label = F, group.by = "clustID")
#DimPlot(hpf04_seurat, reduction = "umap", label = F, group.by = "clustID")
#table(hpf04_seurat@meta.data$clustID)
#Save hpf04 data set
saveRDS(hpf04_seurat, file = "hpf04_seurat.rds")
readRDS(file = "hpf04_seurat.rds")

#Read in hpf06 data
hpf06_clustID <- read.table(file = "GSM3067190_06hpf_clustID.txt")
head(hpf06_clustID)
hpf06 <- read.csv(file = "GSM3067190_06hpf.csv", row.names = 1)
hpf06 <- as.data.frame(hpf06)
hpf06[1:5,1:5]
#Seurat Pipeline
hpf06_seurat <- CreateSeuratObject(counts = hpf06, project = "hpf06_zebrafish", meta.data = hpf06_clustID)
hpf06_seurat <- NormalizeData(hpf06_seurat)
hpf06_seurat <- FindVariableFeatures(hpf06_seurat, selection.method = "vst", nfeatures = 2000)
#Add clusterID
hpf06_seurat$clustID <- hpf06_clustID$V1
#Save hpf06 data set
saveRDS(hpf06_seurat, file = "hpf06_seurat.rds")
readRDS(file = "hpf06_seurat.rds")

#Read in hpf08 data
hpf08_clustID <- read.table(file = "GSM3067191_08hpf_clustID.txt")
head(hpf08_clustID)
hpf08 <- read.csv(file = "GSM3067191_08hpf.csv", row.names = 1)
hpf08 <- as.data.frame(hpf08)
hpf08[1:5,1:5]
#Seurat Pipeline
hpf08_seurat <- CreateSeuratObject(counts = hpf08, project = "hpf08_zebrafish", meta.data = hpf08_clustID)
hpf08_seurat <- NormalizeData(hpf08_seurat)
hpf08_seurat <- FindVariableFeatures(hpf08_seurat, selection.method = "vst", nfeatures = 2000)
#Add clusterID
hpf08_seurat$clustID <- hpf08_clustID$V1
#Save hpf08 data set
saveRDS(hpf08_seurat, file = "hpf08_seurat.rds")
readRDS(file = "hpf08_seurat.rds")

#Read in hpf10 data
hpf10_clustID <- read.table(file = "GSM3067192_10hpf_clustID.txt")
head(hpf10_clustID)
hpf10 <- read.csv(file = "GSM3067192_10hpf.csv", row.names = 1)
hpf10 <- as.data.frame(hpf10)
hpf10[1:5,1:5]
#Seurat Pipeline
hpf10_seurat <- CreateSeuratObject(counts = hpf10, project = "hpf10_zebrafish", meta.data = hpf10_clustID)
hpf10_seurat <- NormalizeData(hpf10_seurat)
hpf10_seurat <- FindVariableFeatures(hpf10_seurat, selection.method = "vst", nfeatures = 2000)
#Add clusterID
hpf10_seurat$clustID <- hpf10_clustID$V1
#Save hpf10 data set
saveRDS(hpf10_seurat, file = "hpf10_seurat.rds")
readRDS(file = "hpf10_seurat.rds")


#Read in hpf14 data
hpf14_clustID <- read.table(file = "GSM3067193_14hpf_clustID.txt")
head(hpf14_clustID)
hpf14 <- read.csv(file = "GSM3067193_14hpf.csv", row.names = 1)
hpf14 <- as.data.frame(hpf14)
hpf14[1:5,1:5]
#Seurat Pipeline
hpf14_seurat <- CreateSeuratObject(counts = hpf14, project = "hpf14_zebrafish", meta.data = hpf14_clustID)
hpf14_seurat <- NormalizeData(hpf14_seurat)
hpf14_seurat <- FindVariableFeatures(hpf14_seurat, selection.method = "vst", nfeatures = 2000)
#Add clusterID
hpf14_seurat$clustID <- hpf14_clustID$V1
#Save hpf14 data set
saveRDS(hpf14_seurat, file = "hpf14_seurat.rds")
readRDS(file = "hpf14_seurat.rds")

#Read in hpf18 data
hpf18_clustID <- read.table(file = "GSM3067194_18hpf_clustID.txt")
head(hpf18_clustID)
hpf18 <- read.csv(file = "GSM3067194_18hpf.csv", row.names = 1)
hpf18 <- as.data.frame(hpf18)
hpf18[1:5,1:5]
#Seurat Pipeline
hpf18_seurat <- CreateSeuratObject(counts = hpf18, project = "hpf18_zebrafish", meta.data = hpf18_clustID)
hpf18_seurat <- NormalizeData(hpf18_seurat)
hpf18_seurat <- FindVariableFeatures(hpf18_seurat, selection.method = "vst", nfeatures = 2000)
#Add clusterID
hpf18_seurat$clustID <- hpf18_clustID$V1
#Save hpf18 data set
saveRDS(hpf18_seurat, file = "hpf18_seurat.rds")
readRDS(file = "hpf18_seurat.rds")

#Read in hpf24 data
hpf24_clustID <- read.table(file = "GSM3067195_24hpf_clustID.txt")
head(hpf24_clustID)
hpf24 <- read.csv(file = "GSM3067195_24hpf.csv", row.names = 1)
hpf24 <- as.data.frame(hpf24)
hpf24[1:5,1:5]
#Seurat Pipeline
hpf24_seurat <- CreateSeuratObject(counts = hpf24, project = "hpf24_zebrafish", meta.data = hpf24_clustID)
hpf24_seurat <- NormalizeData(hpf24_seurat)
hpf24_seurat <- FindVariableFeatures(hpf24_seurat, selection.method = "vst", nfeatures = 2000)
#Add clusterID
hpf24_seurat$clustID <- hpf24_clustID$V1
hpf24_seurat <- ScaleData(hpf24_seurat)
hpf24_seurat <- RunPCA(hpf24_seurat, npcs = 50)
ElbowPlot(hpf24_seurat, ndims = 50)
hpf24_seurat <- FindNeighbors(hpf24_seurat, dims = 1:20)
hpf24_seurat <- FindClusters(hpf24_seurat, resolution = 1)
hpf24_seurat <- RunTSNE(hpf24_seurat, dims = 1:20)
DimPlot(hpf24_seurat, reduction = "tsne", label = T) + NoLegend()
dev.copy2pdf(file='maves_0619_hpf24.seurat_clusters.pdf', width = 7, height = 5)
#Add clusterID
hpf24_seurat$clustID <- hpf24_clustID$V1
DimPlot(hpf24_seurat, reduction = "tsne", label = F, group.by = "clustID")
dev.copy2pdf(file='maves_0619_hpf24.original_clusters.pdf', width = 7, height = 5)
table(hpf24_seurat@meta.data$clustID)
##get marker genes for seurat clusters
hpf24.markers <- FindAllMarkers(hpf24_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hpf24.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#Save hpf24 data set
saveRDS(hpf24_seurat, file = "hpf24_seurat.rds")

##load hpf24 data 
hpf24_seurat <- readRDS(file = "hpf24_seurat.rds")

##cluster expression
cluster.averages <- AverageExpression(hpf24_seurat)
head(cluster.averages[["RNA"]][, 1:5])
write.csv(cluster.averages, file="maves_0619_hpf24_seurat.cluster_gene_expression.csv")

##get marker genes for seurat clusters using just genes from the gene list
gene_list = c('actc1a', 'actc1b', 'actn1', 'alyref', 'angpt1', 'ank3b', 'anp32e', 'arhgap31', 'arhgap44', 'arhgdia', 'arid3b', 'atf3', 'atp1a1a.1', 'atp1a1a.2', 'atp1a1a.3', 'atp1a1a.4', 'atp1a1a.5', 'atp1a1b', 'atp1b1a', 'atp1b1b', 'atp2a2a', 'atp2a2b', 'atp7a', 'bag6', 'bend4', 'bmp10', 'bmp2b', 'bmp4', 'bmp7a', 'bmp7b', 'cacna1da', 'cacna1db', 'cactin', 'cacybp', 'calcrla', 'calcrlb', 'camk2a', 'casz1', 'ccna2', 'ccnd2a', 'ccnd2b', 'ccnj', 'cct5', 'cct6a', 'cdh11', 'cdip1', 'cited2', 'cmtm8b', 'cnn1a', 'col2a1a', 'col2a1b', 'colec12', 'cox4i1', 'cpeb2', 'csnk2b', 'cyc1', 'cycsa', 'cycsb', 'dag1', 'degs1', 'dkc1', 'dpysl3', 'dusp3a', 'dusp3b', 'dynll1', 'ednraa', 'ednrab', 'eef1a2', 'efna5a', 'efna5b', 'efnb2a', 'efnb2b', 'egln1a', 'egln1b', 'eif3ba', 'eif3bb', 'eif4a1a', 'eif4a1b', 'eif6', 'emd', 'enc1', 'eng', 'epha7', 'ets2', 'fasn', 'fbxo32', 'fhl1a', 'fhl1b', 'flt4', 'foxf1', 'fzr1a', 'gabarapa', 'gadd45ga', 'gadd45gb.1', 'gar1', 'gata2a', 'gata2b', 'gata3', 'gata4', 'git1', 'gli3', 'gnb2', 'gnl1', 'golt1ba', 'grpel1', 'hccsa.1', 'hccsb', 'heatr1', 'hif1aa', 'hif1ab', 'hipk2', 'hmgb2a', 'hmgb2b', 'hnrnpaba', 'hnrnpabb', 'hnrnpl2', 'hsd17b10', 'hsp90ab1', 'hspa8', 'hspb7', 'hspe1', 'id2a', 'id2b', 'idh3g', 'igfbp5a', 'igfbp5b', 'jag1a', 'jag1b', 'kif26ba', 'kita', 'kitb', 'kpnb1', 'lpar3', 'map1b', 'mcm6', 'mdh1aa', 'mdh1ab', 'mdh1b', 'mecom', 'med19a', 'med19b', 'meis1a', 'meis1b', 'mepce', 'mfn2', 'mid2', 'mitfa', 'mitfb', 'mlf2', 'mpped2', 'mpped2a', 'mrto4', 'naf1', 'ncoa4', 'notch1a', 'notch1b', 'npm1a', 'npm1b', 'nr2f2', 'nrbp1', 'odc1', 'pa2g4a', 'pa2g4b', 'paip1', 'pard6b', 'park7', 'pcgf5a', 'pcgf5b', 'pdha1a', 'pdha1b', 'pdlim4', 'pdzd11', 'pfn1', 'phf13', 'phf5a', 'pitx1', 'pkma', 'pkmb', 'plk2a', 'plk2b', 'plxnd1', 'pmp22a', 'pmp22b', 'polr2h', 'pomp', 'ppargc1a', 'ppp1caa', 'ppp1cab', 'prkar1aa', 'prkar1ab', 'prkcea', 'prkceb', 'prkg1a', 'prox1a', 'prrx1a', 'prrx1b', 'psma6a', 'psma6b', 'psmd3', 'psmd6', 'ptbp1a', 'ptbp1b', 'ptprda', 'ptprdb', 'rasgrp2', 'rasip1', 'rnd3a', 'rnd3b', 'rnf128a', 'rpl27a', 'rps23', 'rragd', 'rrm2', 'ruvbl2', 'sall4', 'samd4a', 'schip1', 'scube3', 'sf3b5', 'sik1', 'slc25a5', 'slc2a1a', 'slc2a1b', 'sox9a', 'sox9b', 'srfa', 'srfb', 'srsf1a', 'srsf1b', 'stat3', 'tal1', 'tbx3a', 'tbx5a', 'tbx5b', 'tek', 'tgfb2', 'tmem108', 'trps1', 'tsc22d3', 'tshz1', 'tubb5', 'ubald1a', 'ubald1b', 'uqcrfs1', 'vdac1', 'vrtn', 'vsnl1a', 'vsnl1b', 'wnt3', 'wnt5a', 'ykt6', 'ypel5', 'ywhaba', 'ywhabb', 'ywhae1', 'ywhae2', 'ywhaz', 'zbtb20', 'zfhx4', 'znf462', 'znf711')
hpf24.gl_markers <- FindAllMarkers(hpf24_seurat, features = gene_list, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(hpf24.gl_markers, file="maves_0619_hpf24_seurat.gene_list.serurat_clusters_markers.csv")
#change the current cell identities to clustID via and repeat
hpf24_seurat <- SetIdent(hpf24_seurat, value =  "clustID")
hpf24.gl_markers_2 <- FindAllMarkers(hpf24_seurat, features = gene_list, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(hpf24.gl_markers_2, file="maves_0619_hpf24_seurat.gene_list.original_clusters_markers.csv")

#change the current cell identities to clustID via
hpf24_seurat <- SetIdent(hpf24_seurat, value =  "clustID")
cluster.averages <- AverageExpression(hpf24_seurat)
head(cluster.averages[["RNA"]][, 1:5])
write.csv(cluster.averages, file="maves_0619_hpf24_seurat.gene_expression.csv")




##look at clusters etc....
cl_data <- read.table('hpf_24_cluster_exp.txt', header=T, row.names=1)
head(cl_data)
##convert data to percentages
cl_perc = cl_data/rowSums(cl_data[1:72])
head(cl_data[1:72])
head(cl_perc)
##Hierarchical clustering
sampleDists <- dist( cl_perc )
sampleDistMatrix <- as.matrix( sampleDists )
hc <- hclust(sampleDists)
hc
##make graphs
myplclust(hc, labels=rownames(cl_data), lab.col=as.fumeric(as.vector(cl_data$in_genelist)), cex=0.5)
abline(h=0.2)
dev.copy2pdf(file='maves_0619_hpf24.hc.perc.pdf', width = 20, height = 10)
hclusters <- cutree(hc, h=0.2)
cluster_table = table(true=cl_data$in_genelist, cluster=hclusters)
cluster_table
write.table(cluster_table, file="maves_0619_hpf24.hc.perc.counts.txt")
##add cluster id to data
data_with_cid <- cbind(cl_data, clusterID=hclusters)
write.csv(data_with_cid, file="maves_0619_hpf24.hc.cluster_id.csv")
data_pc_with_cid <- cbind(cl_perc, clusterID=hclusters)
write.csv(data_pc_with_cid, file="maves_0619_hpf24.hc.cluster_id.pc.csv")
##heatmap on cluster5 genes.. manually made these
c11_data <- read.table('hpf_24_set1_cl5.txt', header=T, row.names=1)
mat <- c11_data[1:72]
newdf <- c(c11_data$gene, c11_data$in_genelist)
pheatmap(mat, fontsize_row=2, fontsize_col=10, cluster_cols = F)
dev.copy2pdf(file='neural_tube_cluster_data_1_cluster1.pdf', width = 7, height = 5)

##redo just more highly expressed genes
cl_data <- read.table('hpf_24_cluster_exp_2.txt', header=T, row.names=1)
head(cl_data)
##convert data to percentages
cl_perc = cl_data/rowSums(cl_data[1:72])
head(cl_data[1:72])
head(cl_perc)
##Hierarchical clustering
sampleDists <- dist( cl_perc )
sampleDistMatrix <- as.matrix( sampleDists )
hc <- hclust(sampleDists)
hc
##make graphs
myplclust(hc, labels=rownames(cl_data), lab.col=as.fumeric(as.vector(cl_data$in_genelist)), cex=0.5)
abline(h=0.2)
dev.copy2pdf(file='maves_0619_hpf24_set2.hc.perc.pdf', width = 20, height = 10)
hclusters <- cutree(hc, h=0.2)
cluster_table = table(true=cl_data$in_genelist, cluster=hclusters)
cluster_table
write.table(cluster_table, file="maves_0619_hpf24_set2.hc.perc.counts.txt")







##find anchors etc -- not used right now
reference.list <- list(hpf04_seurat, hpf06_seurat, hpf08_seurat, hpf10_seurat, hpf14_seurat, hpf18_seurat, hpf24_seurat)
zf.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
##integrate all timepoint
zf.integrated <- IntegrateData(anchorset = zf.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(zf.integrated) <- "integrated"
## Run the standard workflow for visualization and clustering
zf.integrated <- ScaleData(zf.integrated, verbose = FALSE)
zf.integrated <- RunPCA(zf.integrated, npcs = 30, verbose = FALSE)
zf.integrated <- RunTSNE(zf.integrated, reduction = "pca", dims = 1:30)

DimPlot(zf.integrated, reduction = "tsne", label = F, group.by = "orig.ident")

##merge data using normailzed data -- not used right now
zf.all <- merge(x = hpf04_seurat, y = c(hpf06_seurat, hpf08_seurat, hpf10_seurat, hpf14_seurat, hpf18_seurat, hpf24_seurat), 
                add.cell.ids = c("hpf04", "hpf06", "hpf08", "hpf10", "hpf14", "hpf18", "hpf24"), project = "zf_all", merge.data = TRUE)
zf.all
GetAssayData(pbmc.combined)[1:10, 1:15]
##save
saveRDS(zf.all, file = "zf_all_seurat.rds")
readRDS(file = "zf_all_seurat.rds")
##normalize and then visulize
#zf.all <- NormalizeData(zf.all)
zf.all <- FindVariableFeatures(zf.all, selection.method = "vst", nfeatures = 2000)
zf.all <- ScaleData(zf.all, verbose = FALSE)
zf.all <- RunPCA(zf.all, npcs = 30, verbose = FALSE)
zf.all <- RunTSNE(zf.all, reduction = "pca", dims = 1:30)
DimPlot(zf.all, reduction = "tsne", label = F, group.by = "orig.ident")
DimPlot(zf.all, reduction = "tsne", label = F, group.by = "clustID")
