##############################################################
rm(list = ls())
setwd("/workdir")
library(glmnet)
library(ggplot2)
library(readxl)

data=read_xlsx("data.xlsx",sheet=1)
x <- as.matrix(data[, 3:21])
#y <- data$Pathology
y <- ifelse(data$Pathology == "MIA", 0,1)
######################################
####LASSO####################
######################################
fit <- glmnet(x, y, family = "binomial", alpha=1, maxit = 1000) 
plot(fit,xvar = "lambda", label= TRUE)

set.seed(10)
cvfit <- cv.glmnet(as.matrix(x), as.matrix(y), alpha=1,family="binomial", nfolds=5, maxit = 1000)
cvfit
plot(cvfit)

cvfit$lambda.lse

myCoefs <- coef(cvfit, s=cvfit$lambda.1se)

index = which(myCoefs !=0)
actCofe = myCoefs[index]
lassoGene = row.names(myCoefs)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCofe)
geneCoef
write.csv(geneCoef, file = "geneCoef.csv",row.names = T)

myCoefs2 <- data.frame(Label = rownames(myCoefs), Coefficient = myCoefs[,1]) %>% filter(Coefficient !="0")
myCoefs1 <- myCoefs2[-1,]
coef_plot <- ggplot(myCoefs1, aes(x = Label, y = Coefficient, color=Label)) +
  geom_bar(stat="identity",aes(fill=Label)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Coefficients at Best Lambda") +
  theme_classic() +
  theme( 
    #legend.position="none",
    text=element_text(family = "sans",size=10),
    axis.text=element_text(family = "sans",size=10),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
coef_plot
