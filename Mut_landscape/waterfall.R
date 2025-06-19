rm(list = ls())
setwd("/workdir")
#scatterplot
library(ggplot2)
library(readxl)
library(dplyr)
library(maftools)

laml.maf = "Mut_data.xlsx"

mafdf=read_xlsx(laml.maf,sheet="SNVINDEL")

clindf<-read_xlsx(laml.maf,sheet="TMB&MSI")

laml = read.maf(maf = mafdf, clinicalData = clindf)

#Shows sample summry.
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#Shows all fields in MAF
getFields(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')
#plotmafSummary
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#show top20 genes.
oncoplot(maf = laml, top = 20)
################plot all########################
#choose RColorBrewer
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
                   'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

annocolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')[c(2,8)]
names(annocolors) = c("MIA", "IA")
annocolors = list(Pathological_subtype = annocolors)
sample_order=clindf %>% arrange(Pathological_subtype) %>% select(Tumor_Sample_Barcode) %>% unlist()
top30genes=head(laml@gene.summary$Hugo_Symbol,n=30)
oncoplot(maf = laml, 
         genes=top30genes,
         colors = vc_cols, top = 30,
         fill=TRUE,
         clinicalFeatures = 'Pathological_subtype',
         sortByAnnotation=TRUE,
         groupAnnotationBySize=TRUE,
         showTumorSampleBarcodes=FALSE,
         annotationColor = annocolors,
         keepGeneOrder=TRUE,
         anno_height=0.3)


##############################plot IA group###################
clindf_IA<-clindf %>% filter(Pathological_subtype=="IA")
mafdf_IA<- mafdf %>% filter(Tumor_Sample_Barcode %in% clindf_IA$Tumor_Sample_Barcode)


clindf_MIA<-clindf %>% filter(Pathological_subtype=="MIA")
mafdf_MIA<- mafdf %>% filter(Tumor_Sample_Barcode %in% clindf_MIA$Tumor_Sample_Barcode)

####plot only IA
lamlIA = read.maf(maf = mafdf_IA, clinicalData = clindf_IA)
plotmafSummary(maf = lamlIA, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
genes30=lamlIA@gene.summary$Hugo_Symbol
o=c(1,2,3,7,4,6,5,53,8,10,85,9,11,15,12,14,13,18,19,20,17,24,16,28,21,290,51,22,31,33)
oncoplot(maf = lamlIA, 
         genes=genes30[o], 
         colors = vc_cols, 
         top = 30,
         fill=TRUE,
         showTumorSampleBarcodes=FALSE,
         annotationColor = annocolors,
         keepGeneOrder=TRUE
         )


##############################plot MIA group###################
lamlMIA = read.maf(maf = mafdf_MIA, clinicalData = clindf_MIA)
plotmafSummary(maf = lamlMIA, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
S=c(1,10,4,2,8,5,12,3,23,6,67,26,69,37,97,53,35,61,70,19,16,102,7,9,160,21,22)
oncoplot(maf = lamlMIA, 
         genes=genes27_MIA[S], 
         colors = vc_cols, 
         top = 27,
         fill=TRUE,
         showTumorSampleBarcodes=FALSE,
         annotationColor = annocolors,
         #sampleOrder=sample_order
         keepGeneOrder=TRUE
)
