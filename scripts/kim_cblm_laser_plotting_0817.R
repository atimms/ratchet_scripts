##move to working directory
setwd('/data/atimms/kim_laser_graphing_0818')

##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)

##make graph with hapmap control data and sample
k6_r5_hgdp <- read.csv('cblm_0818.k6_r5.hgdp_all.txt', sep ='\t')
k10_r5_hgdp <- read.csv('cblm_0818.k10_r5.hgdp_all.txt', sep ='\t')
k10_r10_hgdp <- read.csv('cblm_0818.k10_r10.hgdp_all.txt', sep ='\t')
#k4_hgdp <- read.csv('dan_test.hgdp_all.txt', sep ='\t')
ggplot(aes(x = pc1, y = pc2, color = continent), data = k6_r5_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('k6_r5.pc1_pc2.pdf')
ggplot(aes(x = pc3, y = pc4, color = continent), data = k6_r5_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('k6_r5.pc3_pc4.pdf')

ggplot(aes(x = pc1, y = pc2, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('k10_r10.pc1_pc2.pdf')
ggplot(aes(x = pc3, y = pc4, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('k10_r10.pc3_pc4.pdf')
ggplot(aes(x = pc5, y = pc6, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('k10_r10.pc5_pc6.pdf')
ggplot(aes(x = pc7, y = pc8, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('k10_r10.pc7_pc8.pdf')


##2019 version
##make graph with hapmap control data and sample
k6_r5_hgdp <- read.csv('cblm_0219.k6_r5.hgdp_all.txt', sep ='\t')
k10_r5_hgdp <- read.csv('cblm_0219.k10_r5.hgdp_all.txt', sep ='\t')
k10_r10_hgdp <- read.csv('cblm_0219.k10_r10.hgdp_all.txt', sep ='\t')
#k4_hgdp <- read.csv('dan_test.hgdp_all.txt', sep ='\t')
ggplot(aes(x = pc1, y = pc2, color = continent), data = k6_r5_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('cblm_0219.k6_r5.pc1_pc2.pdf')
ggplot(aes(x = pc3, y = pc4, color = continent), data = k6_r5_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('cblm_0219.k6_r5.pc3_pc4.pdf')

ggplot(aes(x = pc1, y = pc2, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('cblm_0219.k10_r10.pc1_pc2.pdf')
ggplot(aes(x = pc3, y = pc4, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('cblm_0219.k10_r10.pc3_pc4.pdf')
ggplot(aes(x = pc5, y = pc6, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('cblm_0219.k10_r10.pc5_pc6.pdf')
ggplot(aes(x = pc7, y = pc8, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
ggsave('cblm_0219.k10_r10.pc7_pc8.pdf')

