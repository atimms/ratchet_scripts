library(ggplot2)


workingDir = "/home/atimms/fetal_0816";
setwd(workingDir);

##graph individual aaf results
#dat <- read.csv("fetal_genomes_0816.dbsnp_passed_rpts_snps.combined.for_r.txt", sep ='\t') 

filename <- dir(workingDir, pattern ="for_r.txt")
for(i in 1:length(filename)){
  dat <- read.csv(filename[i], sep ='\t')
  title <- gsub(".aaf_for_r.txt","",filename[i]) 
  ggplot(dat, aes(aaf, coverage)) +geom_point(color="#0072B2", size = 0.5, position=position_jitter(width=0.02,height=5)) + theme_bw() + ggtitle(title)
  gene_pdf <- gsub(".txt",".pdf",filename[i]) 
  ggsave(gene_pdf)
}


dat <- read.csv("fetal_genomes_0816.dbsnp_passed_rpts_snps.Heart.for_r.txt", sep ='\t')

ggplot(dat, aes(aaf, coverage))+geom_point(color="#0072B2", position=position_jitter(width=0.02,height=5), size = 0.5)

ggplot(dat, aes(aaf, coverage))+geom_point(color="#0072B2")
