library("ggplot2")


workingDir = "/data/atimms/kim_rnaseq_0717";
setwd(workingDir);

#de_data = read.csv('kim_rnaseq_1016.de.egl_vs_pk.asd176.csv', header=T, row.names=1)
#de_data = read.csv('kim_rnaseq_1016.de.egl_vs_pk.chd.csv', header=T, row.names=1)
#de_data = read.csv('kim_rnaseq_1016.de.egl_vs_pk.epilesy.csv', header=T, row.names=1)
#de_data = read.csv('kim_rnaseq_0517.pcl_egl.de.asd162.csv', header=T, row.names=1)
#de_data = read.csv('kim_rnaseq_0717.pcl_egl.de.asd161.csv', header=T, row.names=1)
#de_data = read.csv('kim_rnaseq_0717.pcl_egl.de.scz156.csv', header=T, row.names=1)
## try sorting the de_data first
de_data = de_data[order(de_data$genelist, decreasing = TRUE),]
cbPalette <- c("black", "gray")
ggplot(de_data, aes(x=log2FoldChange, y=-log10(padj), color = genelist)) + geom_point(size = 3) +
  theme_bw(base_size = 12) + scale_colour_manual(values=cbPalette)
#ggsave('kim_rnaseq_1016.pk_gr_de.176gene.volcano.pdf')
#ggsave('kim_rnaseq_1016.pk_gr_de.chd.volcano.pdf')
#ggsave('kim_rnaseq_0517.pcl_egl.asd162.volcano.pdf')
#ggsave('kim_rnaseq_0717.pcl_egl.asd161.volcano.pdf')
#ggsave('kim_rnaseq_0717.pcl_egl.scz156.volcano.pdf')

workingDir = "/data/atimms/acomy_rnaseq_0618";
setwd(workingDir);
filename <- dir(workingDir, pattern =".volcano.csv")
for(i in 1:length(filename)){
  de_data = read.csv(filename[i], header=T)
  de_data = de_data[order(de_data$genelist, decreasing = TRUE),]
  cbPalette <- c("black", "gray")
  ggplot(de_data, aes(x=log2FoldChange, y=-log10(padj), color = genelist)) + geom_point(size = 3) +
    theme_bw(base_size = 12) + scale_colour_manual(values=cbPalette)
  gene_pdf <- gsub(".csv",".pdf",filename[i]) 
  ggsave(gene_pdf)
}



workingDir = "/data/atimms/milena_rnaseq_0818";
setwd(workingDir);
filename <- dir(workingDir, pattern =".volcano.csv")
for(i in 1:length(filename)){
  de_data = read.csv(filename[i], header=T)
  de_data = de_data[order(de_data$label, decreasing = TRUE),]
  cbPalette <- c("green", "gray", "black")
  ggplot(de_data, aes(x=log2FoldChange, y=-log10(padj), color = label)) + geom_point(size = 3) +
    theme_bw(base_size = 12) + scale_colour_manual(values=cbPalette)
  gene_pdf <- gsub(".csv",".pdf",filename[i]) 
  ggsave(gene_pdf)
}

workingDir = "/data/atimms/cherry_rnaseq_1019";
setwd(workingDir);
filename <- dir(workingDir, pattern =".volcano.csv")
for(i in 1:length(filename)){
  de_data = read.csv(filename[i], header=T)
  de_data = de_data[order(de_data$genelist, decreasing = TRUE),]
  cbPalette <- c("red", "gray")
  ggplot(de_data, aes(x=log2FoldChange, y=-log10(padj), color = genelist)) + geom_point(size = 3) +
    theme_bw(base_size = 12) + scale_colour_manual(values=cbPalette)
  gene_pdf <- gsub(".csv",".pdf",filename[i]) 
  ggsave(gene_pdf)
}

