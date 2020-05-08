library(pheatmap)


workingDir = "/data/atimms/dave_heatmap_0220";
setwd(workingDir);

##first one
dave_data <- read.table('dave_data_0220.txt', header=T, row.names=1)
head(dave_data)
##make heatmap
pheatmap(dave_data, show_rownames = F, fontsize_col=8)
dev.copy2pdf(file='osteo_chondro_subcluster_top_500_s_het.pdf')

##get genename per branch
pl = pheatmap(dave_data, fontsize_row=3, fontsize_col=8)
hc <- pl$tree_row
lbl <- cutree(hc, 2)
gn = which(lbl==2)
write.csv(gn, file="genenames.csv")
##redo heatmap with labels so can get in order manually
pheatmap(dave_data, fontsize_row=1, fontsize_col=8)
dev.copy2pdf(file='temp.pdf', width = 20, height = 10)


##second one
dave_data <- read.table('dave_data_0320.txt', header=T, row.names=1)
head(dave_data)
##make heatmap
pheatmap(dave_data, fontsize_row=2, fontsize_col=4)
dev.copy2pdf(file='decile_1_2_top_20.pdf', width = 40, height = 20)

##third one
dave_data <- read.table('dave_index_genes_0320.txt', header=T, row.names=1)
head(dave_data)
##make heatmap
pheatmap(dave_data, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='index_genes_0320.pdf', width = 20, height = 10)
##log +1 and repeat
dave_data_log = log(dave_data+1)
head(dave_data)
head(dave_data_log)
pheatmap(dave_data_log, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='index_genes_0320.log_plus1.pdf', width = 20, height = 10)
pheatmap(dave_data_log, fontsize_row=4, fontsize_col=8, cluster_cols=F)
dev.copy2pdf(file='index_genes_0320.log_plus1.no_col_cluster.pdf', width = 20, height = 10)

my_heatmap = pheatmap(dave_data_log, fontsize_row=4, fontsize_col=8, cluster_cols=F)
##gets genes number
order = my_heatmap$tree_row$order
genes = row.names(dave_data_log)
og <- rbind(order, genes)
write.csv(og, file="index_genes_0320.gene_order.csv")
