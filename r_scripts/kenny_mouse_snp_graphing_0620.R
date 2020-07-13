
##just use common library, so check which one to use and then set that parameter
.libPaths()
.libPaths( .libPaths()[2] )
##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)   


##mapping
##move to working directory
setwd('/archive/beier_d/non_center_data/off_campus_data/kenny_enu_0620/wgs')

hom_data <- read.table('kenny_wgs_0620.chr17_20-30mb100kb.file_for_graphing.txt', header=T)
ggplot(data=hom_data,aes(x=position, y=count)) + geom_bar(stat="identity", colour="black") + facet_grid(type ~ ., scales="free", space="free_x")

##loop through files
filename <- dir('/archive/beier_d/non_center_data/off_campus_data/kenny_enu_0620/wgs', pattern ="file_for_graphing.txt")
for(i in 1:length(filename)){
  ##read in hom data
  hom_data <- read.csv(filename[i], sep ='\t')
  ##bar graph
  ggplot(data=hom_data,aes(x=position, y=count)) + geom_bar(stat="identity", colour="black") + facet_grid(type ~ ., scales="free", space="free_x")
  gene_pdf <- gsub(".file_for_graphing.txt",".pdf",filename[i])
  ggsave(gene_pdf, width=10, height = 5)
  #ggsave(paste(filename[i], 'pdf',sep="."))
}

