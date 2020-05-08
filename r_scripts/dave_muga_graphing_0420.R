##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)


setwd('/archive/beier_d/dave_data/dave_muga_0420')

##testing on one
hom_data <- read.csv('smo_tens.mm10_2000kb_1000kb.bed', sep = '\t')
goodChrOrder <- paste("chr",c(1:19,"X","Y"),sep="")
hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
hom_data
ggplot(data=hom_data,aes(x=start, y=ten_count),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid( ~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())

##loop through files
filename <- dir('/archive/beier_d/dave_data/dave_muga_0420', pattern ="kb.bed")
for(i in 1:length(filename)){
  ##read in hom data
  hom_data <- read.csv(filename[i], sep ='\t')
  ##reorder chromosome, if using 'chr'
  goodChrOrder <- paste("chr",c(1:19,"X","Y"),sep="")
  hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  ##bar graph
  ggplot(data=hom_data,aes(x=start, y=ten_count),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid( ~ chr, scales="free", space="free_x") + 
    theme_bw() + theme(axis.text.x = element_blank())
  gene_pdf <- gsub(".bed",".pdf",filename[i])
  ggsave(gene_pdf, width = 20)
}
