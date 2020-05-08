library(monocle)
library(Matrix)
library(ggplot2)

setwd('/data/atimms/cherry_single_cell_0717/seurat_de_1217')


##all samples and data --- add new data
HSMM_expr_matrix <- read.table("retina_seurat_0817.raw_martix.txt",header=T,row.names = 1)
HSMM_sample_sheet <- read.table("retina_seurat_0817.sample_sheet.txt",header=T,row.names = 1)
HSMM_gene_annotation <- read.table("retina_seurat_0817.gene_ann.txt",header=T,row.names = 1)
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                       expressionFamily=negbinomial.size())
HSMM = estimateSizeFactors(HSMM)
HSMM = estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))


##for each of the three clusters 0,3 and 6 subset the monocle object and get the genes that are expressed in 20% of the total cells
HSMM_cluster0 <- HSMM[expressed_genes, pData(HSMM)$s_cluster == "sc0"]
HSMM_cluster0 <- detectGenes(HSMM_cluster0, min_expr=0.1)
print(head(fData(HSMM_cluster0)))
##20% of the total cell number for this cluster
twenty_pc = 1243*0.2
expressed_genes_sc0 <- row.names(subset(fData(HSMM_cluster0), num_cells_expressed >= twenty_pc))
write.csv(expressed_genes_sc0, file="retina_seurat_0817.sc0.20pc_genelist.csv")
ten_pc = 1243*0.1
expressed_genes_sc0 <- row.names(subset(fData(HSMM_cluster0), num_cells_expressed >= ten_pc))
write.csv(expressed_genes_sc0, file="retina_seurat_0817.sc0.10pc_genelist.csv")
five_pc = 1243*0.05
expressed_genes_sc0 <- row.names(subset(fData(HSMM_cluster0), num_cells_expressed >= five_pc))
write.csv(expressed_genes_sc0, file="retina_seurat_0817.sc0.05pc_genelist.csv")
expressed_genes_sc0_10 <- row.names(subset(fData(HSMM_cluster0), num_cells_expressed >= 10))
write.csv(expressed_genes_sc0_10, file="retina_seurat_0817.sc0.10cells_genelist.csv")


HSMM_cluster3 <- HSMM[expressed_genes, pData(HSMM)$s_cluster == "sc3"]
HSMM_cluster3 <- detectGenes(HSMM_cluster3, min_expr=0.1)
print(head(fData(HSMM_cluster3)))
##20% of the total cell number for this cluster
twenty_pc = 433*0.2
expressed_genes_sc3 <- row.names(subset(fData(HSMM_cluster3), num_cells_expressed >= twenty_pc))
write.csv(expressed_genes_sc3, file="retina_seurat_0817.sc3.20pc_genelist.csv")
ten_pc = 433*0.1
expressed_genes_sc3 <- row.names(subset(fData(HSMM_cluster3), num_cells_expressed >= ten_pc))
write.csv(expressed_genes_sc3, file="retina_seurat_0817.sc3.10pc_genelist.csv")
five_pc = 433*0.05
expressed_genes_sc3 <- row.names(subset(fData(HSMM_cluster3), num_cells_expressed >= five_pc))
write.csv(expressed_genes_sc3, file="retina_seurat_0817.sc3.05pc_genelist.csv")
expressed_genes_sc3_10 <- row.names(subset(fData(HSMM_cluster3), num_cells_expressed >= 10))
write.csv(expressed_genes_sc3_10, file="retina_seurat_0817.sc3.10cells_genelist.csv")

HSMM_cluster6 <- HSMM[expressed_genes, pData(HSMM)$s_cluster == "sc6"]
HSMM_cluster6 <- detectGenes(HSMM_cluster6, min_expr=0.1)
print(head(fData(HSMM_cluster6)))
##20% of the total cell number for this cluster
twenty_pc = 343*0.2
expressed_genes_sc6 <- row.names(subset(fData(HSMM_cluster6), num_cells_expressed >= twenty_pc))
write.csv(expressed_genes_sc6, file="retina_seurat_0817.sc6.20pc_genelist.csv")
write.csv(fData(HSMM_cluster6), file='temp.csv')
ten_pc = 343*0.1
expressed_genes_sc6 <- row.names(subset(fData(HSMM_cluster6), num_cells_expressed >= ten_pc))
write.csv(expressed_genes_sc6, file="retina_seurat_0817.sc6.10pc_genelist.csv")
five_pc = 343*0.05
expressed_genes_sc6 <- row.names(subset(fData(HSMM_cluster6), num_cells_expressed >= five_pc))
write.csv(expressed_genes_sc6, file="retina_seurat_0817.sc6.05pc_genelist.csv")
expressed_genes_sc6_10 <- row.names(subset(fData(HSMM_cluster6), num_cells_expressed >= 10))
write.csv(expressed_genes_sc6_10, file="retina_seurat_0817.sc6.10cells_genelist.csv")