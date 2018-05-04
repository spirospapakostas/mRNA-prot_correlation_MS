#read file
data1=read.table("4_G_express_LFQimputed_RNAseq_isoforms.txt",h=T, sep="\t")
str(data1)
dim(data1)
names(data1)

#reshape in long format
library(reshape2)
aa=melt(data1[,1:36],id=c("ID"))#melt for proteins, dataset 1
plot(data1$pOT6.1,subset(aa,variable=="pOT6.1")$value)#check some values
bb=melt(data1[,c(1,37:71)],id=c("ID"))#melt for mRNA, datset 2
data2=cbind(aa,bb[,-c(1)])#combine 2 datasets

head(data2)
names(data2)[2:5]=c("pop","prot","pop2","GE")
plot(data2$GE,data2$prot)

#create treatment and pop variables
data2$pop3=data2$pop
data2$treatment=data2$pop
pos=gregexpr("\\.[^\\.]*$", data2$treatment[1])[[1]][1]#find the .
data2$treatment=substr(data2$treatment,pos-1,pos-1)#extract treatment, "1" is 10 degrees
data2$treatment=as.factor(data2$treatment)#convert as factor
data2$pop=substr(data2$pop2,2,pos-2)#extract population
data2$pop=as.factor(data2$pop)#convert as factor
levels(data2$pop)#check levels
levels(data2$treatment)
names(data2)
table(data2$treatment,data2$ID)

#other visual checks
summary(data2)
summary(data1)
plot(data1$tOT6.1,subset(data2,pop2=="tOT6.1")$GE)
plot(data1$pVA10.3,subset(data2,pop3=="pVA10.3")$prot)

#the format of the new dataset is appropriated for mixed model analyses
