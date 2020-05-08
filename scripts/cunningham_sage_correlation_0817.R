##just needs to be done once
source("http://depot.sagebase.org/CRAN.R")
pkgInstall("synapseClient")
install.packages("gplots")
install.packages("VennDiagram")
install.packages("githubr")
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
##issue with ConsensusClusterPlus and DGCA
############################
library(synapseClient)
library(gplots)
library(MASS)
library("VennDiagram")
library('githubr')
library('RColorBrewer')
library(biomaRt)
##issue with ConsensusClusterPlus and DGCA
library(ConsensusClusterPlus)
library('DGCA')



synapseLogin(username="atimms", password="26merliN26", rememberMe=TRUE)
setwd("/Users/atimms/Desktop/ngs_data/cunninham_sage_corr_0817")
getwd()

# KKD for Sage Bionetworks
# Jan. 14, 2014
# biomart query convenience functions

library(biomaRt) 
#Hs = useMart("ensembl") # use this one normally
Hs=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)


getByBiotype=function(biotype="protein_coding", gene=TRUE){
  if (gene==TRUE){
    return(getBM(attributes=c("ensembl_gene_id"),filters="biotype",values=biotype, mart=Hs))    
  }
  else {
    return(getBM(attributes=c("ensembl_transcript_id"),filters="transcript_biotype",values=biotype, mart=Hs))
  }
}

getHGNC=function(inENSG,gene=TRUE){
  geneNames = getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),mart=Hs)    
  if (gene==TRUE){
    return(geneNames[match(inENSG,geneNames[,1]),])    
  }
  else { # write this condition
  }
}

getGeneLengths=function(){
  temp = getBM(attributes=c("ensembl_gene_id", "start_position", "end_position"), mart=Hs) 
  temp$length = temp$end_position - temp$start_position
  return(temp) 
}

getGenesForGOOffspring=function(go_acc,go="CC"){
  library('GO.db')
  temp = switch(go,
                CC = as.list(GOCCOFFSPRING),
                BP = as.list(GOBPOFFSPRING),
                MF = as.list(GOMFOFFSPRING))
  
  offspringGOTerms = temp[[which(names(temp) == go_acc)]]
  allENSG = sapply(offspringGOTerms, getGenesForGOTerm)
  return(unique(unlist(allENSG)))
}


getGenesForGOTerm=function(go_acc){
  return(unique(getBM(attributes=c("ensembl_gene_id"), filters="go_id", values = go_acc, mart=Hs)))
}


addBiotype=function(inCounts,gene=TRUE){
  if (gene==TRUE){
    biotypes = getBM(attributes=c("ensembl_gene_id", "gene_biotype"),filters="ensembl_gene_id",values=rownames(inCounts), mart=Hs)
  }
  else {
    biotypes = getBM(attributes=c("ensembl_transcript_id", "transcript_biotype"),filters="ensembl_transcript_id",values=rownames(inCounts), mart=Hs)
  }
  inCounts$biotype = rep("NA", nrow(inCounts))
  inCounts$biotype = biotypes$gene_biotype[match(rownames(inCounts), biotypes$ensembl_gene_id)]
  return(inCounts)
}

addEntrez=function(inCounts,gene=TRUE){
  if (gene==TRUE){
    entrez = getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters="ensembl_gene_id",values=rownames(inCounts), mart=Hs)
    inCounts$entrez = rep("NA", nrow(inCounts))
    inCounts$entrez = entrez$entrezgene[match(rownames(inCounts), entrez$ensembl_gene_id)]
  }
  else {
    entrez = getBM(attributes=c("ensembl_transcript_id", "entrezgene"),filters="ensembl_transcript_id",values=rownames(inCounts), mart=Hs)
    inCounts$entrez = rep("NA", nrow(inCounts))
    inCounts$entrez = entrez$entrezgene[match(rownames(inCounts), entrez$ensembl_transcript_id)]
  }
  return(inCounts)
}

##get data
resid = read.delim(getFileLocation(synGet("syn8555302")),row.names = 1)
resid[1:5,1:6]
##
genesOfInterest = c("RUNX2", "IGF1", "TWIST1", "GSK3b")
geneNames = getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters="hgnc_symbol",values=genesOfInterest,mart=Hs)    
print(geneNames)

temp = data.matrix(resid[which(rownames(resid) %in% geneNames$ensembl_gene_id),])
rownames(temp) = geneNames$hgnc_symbol[match(rownames(temp),geneNames$ensembl_gene_id)]


heatmap.2(temp,trace="none",scale="none", col=bluered,cexRow = 0.9)






