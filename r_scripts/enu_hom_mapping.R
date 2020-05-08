##move to working directory

##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)   

##sergei 1119
##homozygosity mapping
setwd('/data/atimms/sergei_1019')
##loop through files
filename <- dir('/data/atimms/sergei_1019', pattern ="combined_hom_mapping.txt")
for(i in 1:length(filename)){
  ##read in hom data
  hom_data <- read.csv(filename[i], sep ='\t')
  ##reorder chromosome, if using 'chr'
  #goodChrOrder <- paste("chr",c(1:25,"X","Y"),sep="")
  #hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  #goodChrOrder <- paste(c(1:25,"X","Y"),sep="")
  #hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  ##bar graph
  ggplot(data=hom_data,aes(x=chromosome, y=value),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(analysis ~ chr, scales="free", space="free_x") + 
    theme_bw() + theme(axis.text.x = element_blank())
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  ggsave(gene_pdf)
  #ggsave(paste(filename[i], 'pdf',sep="."))
}
hom_data <- read.csv('mut_combined_mm10_10000kb_100kb_combined_hom_mapping.txt')
ggplot(data=hom_data,aes(x=chromosome, y=value),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(analysis ~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
