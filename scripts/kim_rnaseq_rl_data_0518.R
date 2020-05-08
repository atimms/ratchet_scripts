library("DESeq2")

workingDir = "/data/atimms/kim_rnaseq_0218";
setwd(workingDir);


##read in count and metadata
countData1 <- read.table('kim_rnaseq_0218_rl_0518.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0218_rl_0518.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData1$rnaaccess_batch <- as.factor(colData1$rnaaccess_batch)
colData1$age_pcw <- as.factor(colData1$age_pcw)
colData1$donor <- as.factor(colData1$donor)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ sex + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 5 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##rlog transform and check
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)

#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_rl_0518.pca.tissue.pdf') 
plotPCA(rld, intgroup = "age_pcw")
ggsave('kim_rnaseq_0218_rl_0518.pca.age.pdf') 
plotPCA(rld, intgroup = "donor")
ggsave('kim_rnaseq_0218_rl_0518.pca.donor.pdf') 
plotPCA(rld, intgroup = "rnaaccess_batch")
ggsave('kim_rnaseq_0218_rl_0518.pca.batch.pdf') 

