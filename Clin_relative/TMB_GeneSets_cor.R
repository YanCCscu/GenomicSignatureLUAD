######Visualization of Cancer data
rm(list = ls())
setwd("/workdir")
#scatterplot
library(ggplot2)
library(readxl)
library(dplyr)
#scatterplot
library(ggplot2)
library(gplots)
library(dplyr)
library(ggbeeswarm)
library(ggsignif)
library(readxl)
library(tidyr)
library(RColorBrewer)
laml.maf = "data1.xlsx"
###################GeneSets_group with TMB#################################
TBMdf<-read_xlsx(laml.maf,sheet="TMB&MSI")
p1<-ggplot(TBMdf,aes(x=factor(GeneSets_group,levels=c("LowRisk","HighRisk")),y=TMB_nums)) 
Tdiff=t.test(TMB_nums~GeneSets_group,TBMdf)
p1<-p1+
  geom_jitter(aes(color=GeneSets_group,shape=GeneSets_group),width = 0.25, size=4, height = 0.1)+
  scale_color_manual(values = c("#ca1212","#0d9ccb"))+ 
  geom_signif(test = "wilcox.test",
              comparisons = list(c("LowRisk","HighRisk")), 
              y_position = c(max(TBMdf$TMB_nums)),
              map_signif_level = FALSE )+
  ylab("TMB(muts/Mb)")+
  xlab("")+
  theme_classic() +
  theme( 
    #legend.position="none",
    axis.line.x=element_line(linetype=1,color="black",size=1), 
    axis.line.y=element_line(linetype=1,color="black",size=1),
    axis.ticks.x=element_line(color="black",size=1,lineend = 1),
    axis.ticks.y=element_line(color="black",size=1,lineend = 1),
    text=element_text(family = "sans",size=14),
    axis.text=element_text(family = "sans",size=14),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

qfun<-function(x){
  q=quantile(x, c(0.5,0.25,0.75))
  qdf=data.frame(y=q[1],ymin=q[2],ymax=q[3])
  return(qdf)
}

p1<- 
p1+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = 'crossbar', width = 0.5, linewidth = 0.25, color = 'grey10') +
  stat_summary(fun.data =qfun, #median_hilow(x, conf.int=.5, na.rm=TRUE) ,
               geom = 'errorbar', width = 0.25, color = 'grey10')
ggsave(p1,file='TMB_GeneSets_StageI.pdf')


##########################################################################
###########################boxplot with CTs###############################
##########################################################################
Cancerdf<-read_xlsx("data1.xlsx",sheet=1)
genedivdf <- Cancerdf %>% 
  dplyr::select(Tumor_Sample_Barcode,Hugo_Symbol,GeneSets_group,CTs,GeneSets_group,TMB_MB,) %>% filter(CTs !="*") %>% distinct()

genedivdf_V1 <- Cancerdf %>% 
  dplyr::select(Tumor_Sample_Barcode,GeneSets_group,CTs,GeneSets_group,TMB_MB) %>% filter(CTs !="*") %>% distinct()
if(length(unique(genedivdf_V1$GeneSets_group)) >1 && min(table(genedivdf_V1$GeneSets_group)) >1){
  Tdiff=aov(TMB_MB~CTs,genedivdf_V1)
  res=summary(Tdiff)
  maxy=max(genedivdf_V1$TMB_MB)
  genedivdf_V1$CTs = factor(genedivdf_V1$CTs,levels=c("pGGO","mGGO","Solid"))
  genedivdf_V1$GeneSets_group = factor(genedivdf_V1$GeneSets_group,levels=c("LowRisk","HighRisk"))
  p2<-ggplot(data=genedivdf_V1,
             aes(x=GeneSets_group,y=TMB_MB))
  p2<-p2+
    geom_boxplot(aes(color=GeneSets_group,fill=GeneSets_group),width=0.6) +
    geom_jitter(aes(color=GeneSets_group),width = 0.25, height = 0.1,alpha=1)+
    scale_color_manual(values = c("#258eda","#731a6a"))+
    scale_fill_manual(values = c("#87C1EB","#B281AD"))+
    geom_signif(
      comparisons = list(c("LowRisk","HighRisk")), 
      y_position = c(maxy,maxy+0.8,maxy+1.5),
      map_signif_level = FALSE )+
    scale_x_discrete(sprintf("subtype pvalue:%s",signif(res[[1]]$`Pr(>F)`[1],4)), 
                     labels = c("LowRisk" = "L","HighRisk" = "H"))+
    labs(y="TMB(muts/Mb)",
         color="GeneSets_group"
    )+
    theme_classic() +
    theme(
      #legend.position="none",
      text=element_text(family = "sans",size=14),
      axis.text=element_text(family = "sans",size=14),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    facet_wrap(CTs~.,strip.position="top")+
    theme(strip.background = element_blank(), strip.text = element_text(family = "sans",size=18),
          strip.placement = "outside")
}


##########################################################################
###########################boxplot with Grade###############################
##########################################################################
genedivdf <- Cancerdf %>% 
  dplyr::select(Tumor_Sample_Barcode,Hugo_Symbol,GeneSets_group,Grade,GeneSets_group,TMB_MB,) %>% filter(Grade !="*") %>% distinct()

genedivdf_V1 <- Cancerdf %>% 
  dplyr::select(Tumor_Sample_Barcode,GeneSets_group,Grade,GeneSets_group,TMB_MB) %>% filter(Grade !="*") %>% distinct()
if(length(unique(genedivdf_V1$GeneSets_group)) >1 && min(table(genedivdf_V1$GeneSets_group)) >1){
  Tdiff=aov(TMB_MB~Grade,genedivdf_V1)
  res=summary(Tdiff)
  maxy=max(genedivdf_V1$TMB_MB)
  genedivdf_V1$Grade = factor(genedivdf_V1$Grade,levels=c("Grade1","Grade2","Grade3"))
  genedivdf_V1$GeneSets_group = factor(genedivdf_V1$GeneSets_group,levels=c("LowRisk","HighRisk"))
  p2<-ggplot(data=genedivdf_V1,
             aes(x=GeneSets_group,y=TMB_MB))
  p2<-p2+
    geom_boxplot(aes(color=GeneSets_group,fill=GeneSets_group),width=0.6) +
    geom_jitter(aes(color=GeneSets_group),width = 0.25, height = 0.1,alpha=1)+
    scale_color_manual(values = c("#258eda","#731a6a"))+
    scale_fill_manual(values = c("#87C1EB","#B281AD"))+
    geom_signif(
      comparisons = list(c("LowRisk","HighRisk")), 
      y_position = c(maxy,maxy+0.8,maxy+1.5),
      map_signif_level = FALSE )+
    scale_x_discrete(sprintf("subtype pvalue:%s",signif(res[[1]]$`Pr(>F)`[1],4)), 
                     labels = c("LowRisk" = "L","HighRisk" = "H"))+
    labs(y="TMB(muts/Mb)",
         color="GeneSets_group"
    )+
    theme_classic() +
    theme(
      #legend.position="none",
      text=element_text(family = "sans",size=14),
      axis.text=element_text(family = "sans",size=14),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    facet_wrap(Grade~.,strip.position="top")+
    theme(strip.background = element_blank(), strip.text = element_text(family = "sans",size=18),
          strip.placement = "outside")
  ggsave(p2,file='TMB_GeneSets_Grade.pdf')
  p2
}


##########################################################################
###########################boxplot with Tumor_size###############################
##########################################################################
genedivdf <- Cancerdf %>% 
  dplyr::select(Tumor_Sample_Barcode,Hugo_Symbol,GeneSets_group,Tumor_size,GeneSets_group,TMB_MB,) %>% filter(Tumor_size !="*") %>% distinct()

genedivdf_V1 <- Cancerdf %>% 
  dplyr::select(Tumor_Sample_Barcode,GeneSets_group,Tumor_size,GeneSets_group,TMB_MB) %>% filter(Tumor_size !="*") %>% distinct()
if(length(unique(genedivdf_V1$GeneSets_group)) >1 && min(table(genedivdf_V1$GeneSets_group)) >1){
  Tdiff=aov(TMB_MB~Tumor_size,genedivdf_V1)
  res=summary(Tdiff)
  maxy=max(genedivdf_V1$TMB_MB)
  genedivdf_V1$Tumor_size = factor(genedivdf_V1$Tumor_size,levels=c("<=1cm","1-2cm",">2cm"))
  genedivdf_V1$GeneSets_group = factor(genedivdf_V1$GeneSets_group,levels=c("LowRisk","HighRisk"))
  p2<-ggplot(data=genedivdf_V1,
             aes(x=GeneSets_group,y=TMB_MB))
  p2<-p2+
    geom_boxplot(aes(color=GeneSets_group,fill=GeneSets_group),width=0.6) +
    geom_jitter(aes(color=GeneSets_group),width = 0.25, height = 0.1,alpha=1)+
    scale_color_manual(values = c("#258eda","#731a6a"))+
    scale_fill_manual(values = c("#87C1EB","#B281AD"))+
    geom_signif(
      comparisons = list(c("LowRisk","HighRisk")), 
      y_position = c(maxy,maxy+0.8,maxy+1.5),
      map_signif_level = FALSE )+
    scale_x_discrete(sprintf("subtype pvalue:%s",signif(res[[1]]$`Pr(>F)`[1],4)), 
                     labels = c("LowRisk" = "L","HighRisk" = "H"))+
    labs(y="TMB(muts/Mb)",
         color="GeneSets_group"
    )+
    theme_classic() +
    theme(
      #legend.position="none",
      text=element_text(family = "sans",size=14),
      axis.text=element_text(family = "sans",size=14),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    facet_wrap(Tumor_size~.,strip.position="top")+
    theme(strip.background = element_blank(), strip.text = element_text(family = "sans",size=18),
          strip.placement = "outside")
  ggsave(p2,file='TMB_GeneSets_Tumor_size.pdf')
}
