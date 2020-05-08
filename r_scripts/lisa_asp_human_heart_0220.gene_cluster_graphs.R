library(ggplot2)

workingDir = "/data/atimms/lisa_asp_human_heart_0220";
setwd(workingDir);


##from paper clusters
dat = read.table("asp_de_genes_candidate_counts_0220.txt", header = T)

# Map the time of day to different fill colors
ggplot(data=dat, aes(x=cluster, y=genelist)) + geom_bar(stat="identity",fill="#E69F00") + scale_x_continuous(breaks=seq(0,14))
ggsave('asp_0220.cluster_genelist.hist_genelist_count.pdf')

ggplot(data=dat, aes(x=cluster, y=genelist.total)) + geom_bar(stat="identity",fill="#56B4E9")+ scale_x_continuous(breaks=seq(0,14))
ggsave('asp_0220.cluster_genelist.hist_genelist_total_count.pdf')

##from my clusters
##from paper clusters
dat = read.table("asp_de_genes_candidate_counts_new_0220.txt", header = T)

# Map the time of day to different fill colors
ggplot(data=dat, aes(x=cluster, y=genelist)) + geom_bar(stat="identity",fill="#E69F00") + scale_x_continuous(breaks=seq(0,17))
ggsave('asp_0220.cluster_genelist.hist_genelist_count_new.pdf')

ggplot(data=dat, aes(x=cluster, y=genelist_total)) + geom_bar(stat="identity",fill="#56B4E9")+ scale_x_continuous(breaks=seq(0,17))
ggsave('asp_0220.cluster_genelist.hist_genelist_total_count_new.pdf')

