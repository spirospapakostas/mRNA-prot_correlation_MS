#Read file
data1=read.table("4_G_express_LFQimputed_RNAseq_654genes.txt",h=T)
str(data1)
dim(data1)
names(data1)

#####convert data in order to have individuals in row####
IDs=data1$ID#take the ID
names(data1[,2:36])#prot
names(data1[,c(37:71)])#mrna

data_prot=t(data1[,2:36])
data_mRNA=t(data1[,c(37:71)])
head(data_mRNA)
data_mRNA=as.data.frame(data_mRNA)
data_prot=as.data.frame(data_prot)

names(data_mRNA)=IDs
names(data_prot)=IDs

ind=row.names(data_mRNA)
ind2=row.names(data_prot)

ind= substr(ind,2,nchar(ind))#remove p or t
ind2= substr(ind2,2,nchar(ind2))#remove p or t
ind==ind2#same individuals in row

#replace row.names
row.names(data_mRNA)=ind
row.names(data_prot)=ind2

#check class
sapply(data_mRNA,class)
sapply(data_prot,class)

#data_mRNA = apply(data_mRNA, 2, function(x) as.numeric(as.character(x)))
#data_prot = apply(data_prot, 2, function(x) as.numeric(as.character(x)))

#####treatment and pop####
pop3=ind
treatment=ind
pos=gregexpr("\\.[^\\.]*$", treatment)[[1]][1]
treatment=substr(treatment,pos-1,pos-1)
treatment=as.factor(treatment)
levels(treatment)[1]="10"

pop=substr(pop3,1,pos-2)
pop=as.factor(pop)
table(pop)

#####Start analyses of co-inertia####
library(ade4)

pca1=dudi.pca(data_mRNA)#acp sur mRNA, 2 axes, center = TRUE, scale = TRUE automatically
pca2=dudi.pca(data_prot)#acp sur prot, 3 axes

#some plots
par(mfrow=c(1,2))
s.class(pca1$li, fac = treatment, xax = 1, yax = 2)#treatment
s.class(pca2$li, fac = treatment, xax = 1, yax = 2)

par(mfrow=c(1,2))
s.class(pca1$li, fac = pop, xax = 1, yax = 2)#population
s.class(pca2$li, fac = pop, xax = 1, yax = 2)

par(mfrow=c(1,2))
s.class(pca1$li, fac = pop:treatment, xax = 1, yax = 2)#interaction
s.class(pca2$li, fac = pop:treatment, xax = 1, yax = 2)

inertia.dudi(pca1)#55% for 2 axes
inertia.dudi(pca2)#58% for 3 axes

pca1$cw#uniform weights for co-inertia
pca2$cw

#co-inertia
coiner1=coinertia(pca1,pca2)#2
summary(coiner1)
inertia.dudi(coiner1)#globally, 77% with 2 axes

#tests with permutation
set.seed(110489)
rvtest2100k=randtest(coiner1,nrepet = 100000)#observed some var in pvalues because close to 0.05, 10k repetitions would be better?
hist(rvtest2100k$sim)
str(rvtest2100k)
dat=rvtest2100k$sim
dat=as.data.frame(dat)
dat$repetitions=1:100000

qplot(dat$dat, geom="histogram", binwidth = 0.005) +
  geom_segment(aes(x=0.2203426,y=0,xend=0.2203426,yend=5000))+
  geom_point(aes(x=0.2203426,y=5000))+ylab("Count")+xlab("RV")+
  geom_text(aes(y=5200,x=0.2203426, label="Observed RV"),size=6)+
  theme_bw()+
  theme(plot.margin=unit(c(0.5,2.6,0.5,0.3),"cm"),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(10,0,0, 0)),
        axis.text = element_text(size = rel(3)),
        axis.title = element_text(size = rel(3)))

#dev.copy(pdf,"RV_permutations.pdf", width=8, height=8)
#dev.off()

# 
#two sides, exactly same in fact (greater, not 2 sides)
set.seed(110489)
rvtest3=randtest(coiner1,nrepet = 10000,alter="two-sided")#observed some var in pvalues because close to 0.05, 10k repetitions would be better?
plot(rvtest3)
#save(rvtest3,file="rvtest3")

#check things
names(coiner1)

coiner1$eig[1]/sum(coiner1$eig)#var explained by axe 1
coiner1$eig[2]/sum(coiner1$eig)
inertia.dudi(coiner1)#globally, 77% with 2 axes

#plots
s.corcircle(coiner1$aX)#axes acp1
s.corcircle(coiner1$aY)

#Adonis####
#package data
library(vegan)

#test for multivariate homogeneity of group dispersions before ADONIS
#mRNA
dist11<-vegdist(data_mRNA, method="euclidean")#euclidian distance
anova(betadisper(dist11 , treatment))
anova(betadisper(dist11 , treatment:pop))

#with permutations
permutest(betadisper(dist11 , treatment:pop))

#protein
dist11<-vegdist(data_prot, method="euclidean")
anova(betadisper(dist11 , treatment))
anova(betadisper(dist11 , treatment:pop))

#with permutations
permutest(betadisper(dist11 , treatment:pop))

#Multivariate analysis of variance using distance matrix, mRNA
adonis(data_mRNA ~ pop*treatment, permutations=10000)
adonis(data_mRNA ~ pop+treatment, permutations=10000)
adonis(data_mRNA ~ treatment+pop, permutations=10000)#inverse to see (sequential testing), almost the same because balanced design?

#Multivariate analysis of variance using distance matrix, protein
adonis(data_prot ~ pop*treatment,  permutations=10000)
adonis(data_prot ~ pop+treatment,  permutations=10000)
adonis(data_prot~ treatment+pop, permutations=10000)#inverse to see, ok because balanced design?
