rm(list = ls())
setwd("/workdir")
library("org.Hs.eg.db")
library("enrichplot")
library("clusterProfiler")
library("ggplot2")
Filterpvalue=0.05         
Filterqvalue=0.05         
inputfile= "Gene_cor.txt"
colorSel="qvalue"
if(Filterqvalue>0.05){
  colorSel="pvalue"
}
   
inputfile=read.table(inputfile, header=T, sep="\t", check.names=F)    
genes=unique(as.vector(inputfile[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
inputfile=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       

kkenrich=enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kkenrich)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(inputfile$genes[match(strsplit(x,"/")[[1]],as.character(yifeng$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<Filterpvalue & KEGG$qvalue<Filterqvalue),]
write.table(KEGG, file="KEGGEnrichment_gene_cor.txt", sep="\t", quote=F, row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="keggbubble_gene_cor.pdf", width = 9, height = 7)
dotplot(kkenrich, showCategory=showNum, orderBy="GeneRatio", color=colorSel, label_format = 90)
dev.off()

pdf(file="keggbarplot_gene_cor.pdf", width=9, height=7)
barplot(kkenrich, drop=TRUE, showCategory=showNum, color=colorSel, label_format = 90)
dev.off()