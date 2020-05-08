library("ggplot2")
library(reshape2)

workingDir = "/data/atimms/kim_rnaseq_0717";
setwd(workingDir);

##file from python script 
de_data = read.table('kim_rnaseq_1016.4_genes_0517.r_boxplot.txt', header=T, row.names=1)
de_data$tissue <- factor(de_data$tissue, levels=c("EGL", "Whole", "PK"))
##graph for the 4 genes
ggplot(de_data, aes(x=tissue, y=ATOH1)) + geom_boxplot() + theme_bw(base_size = 12)
ggsave("kim_rnaseq_1016.ATOH1.boxplot.pdf")
ggplot(de_data, aes(x=tissue, y=PAX6)) + geom_boxplot()+ theme_bw(base_size = 12)
ggsave("kim_rnaseq_1016.PAX6.boxplot.pdf")
ggplot(de_data, aes(x=tissue, y=SKOR2)) + geom_boxplot()+ theme_bw(base_size = 12)
ggsave("kim_rnaseq_1016.SKOR2.boxplot.pdf")
ggplot(de_data, aes(x=tissue, y=CALB1)) + geom_boxplot()+ theme_bw(base_size = 12)
ggsave("kim_rnaseq_1016.CALB1.boxplot.pdf")


##file from python script -combined version
##box plot 
#de_data = read.table('kim_rnaseq_0717.4_genes.norm.r_boxplot_facet.txt', header=T)
de_data = read.table('kim_rnaseq_0717.4_genes.rlog.r_boxplot_facet.txt', header=T)
de_data$region <- factor(de_data$region, levels=c("EGL", "Bulk", "PCL"))
ggplot(data = de_data, aes(x=region, y=count)) + geom_boxplot() + facet_wrap(~gene,ncol = 2, scales="free_y") + theme_bw(base_size = 12) 
#ggsave("kim_rnaseq_0717.norm_counts.boxplot.pdf")
ggsave("kim_rnaseq_0717.rlog_counts.boxplot.pdf")


##line graph
##rlog
de_data = read.table('kim_rnaseq_0717.4_genes.rlog.r_line_facet.txt', header=T)
de_data$region <- factor(de_data$region, levels=c("EGL", "Bulk", "PCL"))
ggplot(data = de_data, aes(x=age_pcw, y=count,color=region)) +geom_line() + facet_wrap(~gene,ncol = 2, scales="free_y") + theme_bw(base_size = 12) 
ggsave("kim_rnaseq_0717.rlog_counts.line_graph.pdf")
##normalized
de_data = read.table('kim_rnaseq_0717.4_genes.norm.r_line_facet.txt', header=T)
de_data$region <- factor(de_data$region, levels=c("EGL", "Bulk", "PCL"))
ggplot(data = de_data, aes(x=age_pcw, y=count,color=region)) +geom_line() + facet_wrap(~gene,ncol = 2, scales="free_y") + theme_bw(base_size = 12) 
ggsave("kim_rnaseq_0717.norm_counts.line_graph.pdf")