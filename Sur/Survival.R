rm(list = ls())
setwd("/workdir")

library(survival)
library(survminer)
inputfile="data.txt"

biodateSurvival=function(biodateinputFile=null, myoutFile=null){
  datas=read.table(biodateinputFile, header=T, sep="\t", check.names=F)
  diff=survdiff(Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ DFS_Groups, data=datas)
  pValue=1-pchisq(diff$chisq, df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ DFS_Groups, data = datas)
  
  suryifengPlot=ggsurvplot(fit, 
                           data=datas,
                           conf.int=T,
                           pval=pValue,
                           pval.size=6,
                           #surv.median.line = "hv",
                           legend.title="GeneSets11",
                           legend.labs=c("HighRisk", "LowRisk"),
                           xlab="Time(Months)",
                           break.time.by = 20,
                           palette=c("red", "blue"),
                           risk.table=TRUE,
                           risk.table.title="",
                           risk.table.col = "strata",
                           risk.table.height=.35)
  pdf(file=myoutFile, onefile=FALSE, width=7, height=7)
  print(suryifengPlot)
  dev.off()
}

biodateSurvival(inputfile, myoutFile="survival.pdf")
