### launching interactive 
bsub -n 6 -Is -q interactive bash

### loading most recent R
module load stats/R/3.3.1

### start up R
R

# Enter commands in R (or R studio, if installed)
install.packages("devtools")
library(devtools)

install_github("satijalab/seurat")
library(Seurat)

#### Seurat Main Cluster analysis on Orchestra
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

#Inputs
species="h"
projectName = "20161111_Tim400"
base_dir <- "/groups/neuroduo/tim/161109_indrops/"
folder_names <- c("TJC_Hu1_ret_1","TJC_Hu1_ret_2","Hu1_retina_1","Hu1_retina_2","Hu4_retina_1","Hu4_retina_2","Hu4_retina_3","Hu4_retina_4","Hu8_retina_1","Hu8_retina_2")
mapping_folder <-"/groups/neuroduo/tim/161109_indrops/"
pc_num = 20;
res=0.6;
min_genes = 400;
markers = c("RHO", "SAG", "PDE6B", "OPN1MW", "OPN1SW", "OPN1LW", "ARR3", "CABP5", "LHX1", "VSX2", "OTX2", "SCGN", "GRM6", "PRKCA", "ISL1", "APOE", "RLBP1", "PAX6", "GAD1", "SLC6A9", "SLC17A6", "POU4F1", "POU4F2", "POU4F3", "NEFL", "GFAP", "PECAM1", "KCNJ8", "CX3CR1", "ABCA4", "MYO7A", "USH2A", "PCDH15", "RDH12", "RLBP1", "RP1", "RDS", "PRPF8", "PRPF3", "CRB1", "RPE65", "PDE6H", "CRX", "OTX2", "NRL", "RORB", "MEF2D", "PDE6H", "RAB41", "ARR3", "MPP4", "GUCA1C", "COL4A3", "STRADB", "RP11-1069G10.1", "PEX5L", "GNAT2", "CC2D2A", "FSD1L", "RPGRIP1", "TF", "RASGEF1B", "EMB", "GNGT2", "LBH", "MLXIP", "CLU", "PDE6C", "KCNV2", "GNAS", "CHRNA3", "ANO6", "RABEPK", "CCDC136", "RPL31", "FAM98B", "PCAT6", "DDAH1", "ABHD6", "DOLPP1", "DIP2A", "GPR160", "HSP90AA1", "KIAA1549", "RPS19", "KCNK6", "TSPAN6", "LINC01513", "DHRS12", "FSD2", "KDM3A", "SLC25A25")


#Load files
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)
#[1] 25463	1

for (i in folder_names)
{
  counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
  rownames(counts) <- paste(i,rownames(counts),sep = "_")
  counts_t <- as.data.frame(t(counts))
  counts_t <- counts_t[, colSums(counts_t)>min_genes];
  cd <- cbind(cd,counts_t);
  
}

#Annoying bug where last command gets stuck on .counts.tsv.gz files.  Need to be unzipped.  Cound have probably included that above, but did it manually.

cd$test <- NULL

cd_s <- as(as.matrix(cd), "dgCMatrix")
save(cd_s, file = "cd_s.Robj")


#Create a new folder
date = format(Sys.time(),"%Y%m%d")
dir.create(paste(base_dir,date,projectName,"_Seur_R", sep=""))
setwd(paste(base_dir,date,projectName,"_Seur_R", sep=""))

### Seurat
seurat_mat <- new("seurat", raw.data = cd_s)
seurat_mat <- Setup(seurat_mat, min.cells = 3, min.genes = as.numeric(min_genes), do.logNormalize = T, total.expr = 1e4, project = date)

dim(seurat_mat@data.info)
#[1] 4802    3

##Need to modify the matrix to fix a bug of unknown origin
#seurat_mat@data.info <- seurat_mat@data.info[1:length(seurat_mat@ident),]
#rownames(seurat_mat@data.info)=names(seurat_mat@ident)
#seurat_mat@data.info$orig.ident=seurat_mat@ident

### Find mt genes
if (species == "h"){
  mito.genes <- grep("^MT-", rownames(seurat_mat@data), value = T)
} else{
  mito.genes <- grep("^mt.", rownames(seurat_mat@data), value = T)
}

#Calculate mitohondrial genes for each cell
percent.mito <- colSums(expm1(seurat_mat@data[mito.genes, ]))/colSums(expm1(seurat_mat@data))

#Add mito data back
seurat_mat <- AddMetaData(seurat_mat, percent.mito, "percent.mito")

#Plot UMI vs nGene
pdf("GenePlot.pdf")
GenePlot(seurat_mat, "nUMI", "nGene")
dev.off()

#Plot UMI, nGene, mito violins
pdf("VlnPlot.pdf")
VlnPlot(seurat_mat, c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

#Remove cells that have 
seurat_mat <- SubsetData(seurat_mat, subset.name = "percent.mito", accept.high = 0.10)
seurat_mat <- SubsetData(seurat_mat, subset.name = "nUMI", accept.high = 15000)

#Skipped regressing - took too long for large datases
#seurat_mat <- RegressOut(seurat_mat, latent.vars = c("nUMI"))

#ID variable genes
seurat_mat <- MeanVarPlot(seurat_mat ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)

#Print number of variable genes
print(length(seurat_mat@var.genes))

#Do PCA
seurat_mat <- PCA(seurat_mat, pc.genes = seurat_mat@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5,pcs.store=100)
seurat_mat <- ProjectPCA(seurat_mat)

#Plot PCA Heatmap
pdf("PCHeatmap.pdf")
PCHeatmap(seurat_mat, pc.use = 1, cells.use = 100, do.balanced = TRUE)
dev.off()

#seurat_mat <- JackStraw(seurat_mat, num.replicate = 50, do.print = FALSE) 

#pdf("JackStraw.pdf")
#JackStrawPlot(seurat_mat, PCs = 1:30)
#dev.off()

#Plot PC SD
pdf("PCElbowPlot.pdf")
PCElbowPlot(seurat_mat, num.pc = 100)
dev.off()

#head(seurat_mat@pca.x)
seurat_mat <- FindClusters(seurat_mat, pc.use = 1:pc_num, resolution = res, print.output = 1, save.SNN = T, do.sparse = T)
seurat_mat <- RunTSNE(seurat_mat, dims.use = 1:pc_num, do.fast = T)
write.csv(seurat_mat@tsne.rot, file = paste("tSNECoordinates",pc_num,"_",res,".csv",sep=""))

pdf("tSNE.pdf")
TSNEPlot(seurat_mat)
dev.off()

seuratClusters <- data.frame(seurat_mat@ident)
write.csv(seuratClusters, file = "seuratClusters.csv")


## check that the genes we want to plot exist in our data
markers_filt = intersect(rownames(seurat_mat@data),markers)
markers_filt_split <- split(markers_filt, ceiling(seq_along(markers_filt)/4))

for (i in 1:length(markers_filt_split))
{
  pdf(paste("FeaturePlot_",i,".pdf"))
  FeaturePlot(seurat_mat, markers_filt_split[[i]],cols.use = c("grey","blue"))
  dev.off()
}

clusterMarkers <- FindAllMarkers(seurat_mat, only.pos = F, min.pct = 0.1, thresh.use = 0.25)
write.csv(clusterMarkers, file = "clusterMarkers.csv")

save(seurat_mat, file = "Seurat.Robj")

