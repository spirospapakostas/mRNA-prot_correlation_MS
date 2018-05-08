library(MCMCglmm)
load("datag4")

#Variables for the loop accross genes
gel1=c()#Gelman diagnostic, 1 model with quasi=T
gel2=c()#Gelman diagnostic for random terms
ID=c()
autocorr=c()#autocorrelation, but we can also check the final effective sampling sizes
autocorr_RE=c()#random terms autocorrelation
prot_treat=c()#treatment effect on proteins
GE_treat=c()#treatment effect on mRNA
trait_corr=c()#correlation mRNA/proteins
sd1=c()#sd for eventual post-analyses

prior = list(R = list(V = diag(2)/1000000, nu = 3),#inverse gamma 0.01...
             G = list(G1 = list(V = diag(2), nu = 1.002)))

n1=1#first model in the loop
n2=654#last model in the loop


for(i in n1:n2){
  print(i)
  data3=data2[which(data2$ID==levels(data2$ID)[i]),]#select the gene
  data3[,c(3,5)]=scale(data3[,c(3,5)])#scale
  
  #two models for convergence
  mod1=MCMCglmm(cbind(prot,GE)~trait*treatment-1,random=~us(trait):pop:treatment, 
                rcov = ~ us(trait):units, data=data3,
                prior=prior,
                family = rep("gaussian", 2), nitt = 600000, burnin = 100000, thin=4)
  
  mod2=MCMCglmm(cbind(prot,GE)~trait*treatment-1,random=~us(trait):pop:treatment, 
                rcov = ~ us(trait):units, data=data3,
                prior=prior,start=list(QUASI=FALSE),
                family = rep("gaussian", 2), nitt = 600000, burnin = 100000, thin=4)
  
  #some diagnostics
  #Gelman convergence
  gel1=rbind(gel1,sapply(gelman.diag(mcmc.list(mod1$Sol, mod2$Sol)),max))#try to keep it under 1.1
  gel2=rbind(gel2,max(gelman.diag(mcmc.list(mod1$VCV, mod2$VCV),multivariate = F)$psrf[,1]))#try to keep it under 1.1
  
  #autocorrelation
  autocorr=c(autocorr,sapply(autocorr(mcmc.list(mod1$Sol),lags=1),max))
  autocorr_RE=c(autocorr_RE,sapply(autocorr(mcmc.list(mod1$VCV),lags=1),max))
  
  #now we extract interesting features
  #trait correlation
  
  #pMCMC, correlation values
  pval=min(prop.table(table(posterior.cor(mod1$VCV[,5:8])[,2]>0)))*2
  post.mean=mean(posterior.cor(mod1$VCV[,5:8])[,2])#Transforms posterior distribution of covariances into correlations
  lCI=summary(posterior.cor(mod1$VCV[,5:8])[,2], quantiles=c(0.025,0.975))[[2]][1]#lower confidence
  uCI=summary(posterior.cor(mod1$VCV[,5:8])[,2], quantiles=c(0.025,0.975))[[2]][2]#upper
  
  trait_corr=rbind(trait_corr, cbind(post.mean,lCI,uCI,eff_size=NA,pval))#combine the results
  
  sd1=c(sd1,sd(posterior.cor(mod1$VCV[,5:8])[,2]))
  
  #treatment effect
  prot_treat=rbind(prot_treat,summary(mod1)$solution[3,])
  #calculate it for GE/mRNA
  interest1="treatment6"
  interest2="traitGE:treatment6"
  iter=mod1$Sol#all iterations
  values=iter[,pmatch(interest1,colnames(iter))]+iter[,pmatch(interest2,colnames(iter))]#calculate mRNA treatment at each iteration from protein treatment effect and differences
  
  #calculate pvalues and confidence interval by ourself: 
  #GE_quantiles=summary(values,quantiles=c(0.025,0.975))[[2]]
  #GE_mean=summary(values,quantiles=c(0.025,0.975))[[1]][1]
  #pval2=min(prop.table(table(values>0)))*2
  #GE_treat=rbind(GE_treat,cbind(GE_mean,GE_quantiles[1],GE_quantiles[2],effective_sample=NA,p=pval2))

  #or
  #replace it in the model so that it is calculated exactly in the same way (effective size accounted for)
  mod1$Sol[,3]=values
  GE_treat=rbind(GE_treat,summary(mod1)$solution[3,])
  
  ID=c(ID,levels(data2$ID)[i])
}

#save
name="bivar_g4_6"

obj=list(autocorr=autocorr,autocorr_RE=autocorr_RE,sd1=sd1, ID=ID, gel1=gel1,gel2=gel2,
         prot_treat=prot_treat,GE_treat=GE_treat, trait_corr=trait_corr)

save(obj,file=name)

#save(gel,file=paste("gel_",name,sep=""))
#save(autocorr,file=paste("autocorr_",name,sep=""))
#save(autocorr_RE,file=paste("autocorr_RE_",name,sep=""))
#save(models,file=paste("models_",name,sep=""))
#save(prot_treat,file=paste("prot_treat_",name,sep=""))
#save(GE_treat,file=paste("GE_treat_",name,sep=""))
#save(trait_corr,file=paste("trait_corr_",name,sep=""))
