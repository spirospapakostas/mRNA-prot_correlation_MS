#to run sequentially with mRNA (GE variable in dataset) and protein data (prot variable)
library(pbkrtest)
library(lme4)
data_g4=data2
#Empty vectors
coef=c()
pval2=c()
IDs=c()

for(i in 1:length(levels(data_g4$ID))){#loop for each gene
  data_g4b=data_g4[which(data_g4$ID==levels(data_g4$ID)[i]),]#select one gene ID
  data_g4b[,c(3,5)]=scale(data_g4b[,c(3,5)])#scale
  
  lmer1=lmer(GE~treatment+(1|pop:treatment),data=data_g4b, REML=FALSE)#for anova, REML=FALSE
  lmer2=lmer(GE~1+(1|pop:treatment),data=data_g4b,REML=FALSE)
  
  coef=c(coef,coefficients(lmer1)$`pop:treatment`[1,2])
  pval2=c(pval2, KRmodcomp(lmer1,lmer2)$test$p.value[1])#, KRmodcomp refit the model or not depending on the REML option, extract pvalue
  IDs=c(IDs,levels(data_g4$ID)[i])#take gene ID
}

#see: http://glmm.wikidot.com/faq for pvalues

hist(pval2)

#select significant fdr adjusted pvalue
which(p.adjust(pval2,method="fdr")<0.05)


