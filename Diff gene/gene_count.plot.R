rm(list = ls())
setwd("/workdir")
#scatterplot
library(ggplot2)
library(gplots)
library(dplyr)
library(ggbeeswarm)
library(ggsignif)
library(readxl)
library(tidyr)
library(RColorBrewer)

Cancerdf<-read_xlsx("data.xlsx",sheet=1)

rmdup_cancerdf<-Cancerdf %>% select(Tumor_Sample_Barcode,Hugo_Symbol ,Pathological_Subtype ) %>% distinct()
subtype_counts<-table(rmdup_cancerdf$Pathological_Subtype)
IAcounts=as.vector(subtype_counts)[1]
MIAcounts=as.vector(subtype_counts)[2]

GeneCounts=as.matrix(table(rmdup_cancerdf[,c("Hugo_Symbol","Pathological_Subtype")]))

chisqout<-apply(GeneCounts,1,function(x){
  res=chisq.test(matrix(c(x[1],x[2],IAcounts-x[1],MIAcounts-x[2]),2,2));
  return(c(x[1],x[2],res$p.value))
})
chisqout<-as.data.frame(t(chisqout))
names(chisqout)[3]="Pvalue"

sigGenedf=chisqout %>% filter(Pvalue<0.05)
sigGenedf
keepgenes<-row.names(sigGenedf)
write.table(sigGenedf, file = "Gene_chisqout_results_wenjiang.diff_test.xls", sep = "\t", quote = F, row.names = T)