rm(list=ls())
setwd("/home/ghepard/deeper_spike/14_10_2020/")
library(RColorBrewer)
library(wesanderson)
library(pheatmap)
library(reshape2)
library(ape)
library(gplots)
library(dplyr)
library(plyr)
library(data.table)

dd<-read.csv("new_tabs/MutationPatterns.csv",sep=";",stringsAsFactors = F)
rownames(dd)<-dd$X
dd$X<-NULL
ff<-dd %>% replace(is.na(.),0)
tt<-data.frame(t(ff))
###
sel<-grep("BRASP",rownames(tt),value = TRUE)
sel2<-grep("BAL",rownames(tt),value = TRUE)
list<-c(sel,sel2)
tt$patient<-ifelse(rownames(tt) %in% list,"lower","upper")
aa<-melt(tt)
##### STAT
aa$stat<-ifelse(aa$patient=="lower",0,1)
#wilcox.test(aa$value,aa$stat)

df<-data.frame()
for (i in aa$variable){
  ss<-aa[aa$variable==i,]
  wilc<-wilcox.test(ss$value~as.factor(ss$stat))$p.value
  #print (paste(aa[aa$variable==i,]$variable[1],wilc,sep=" "))
  df[i,1]<-aa[aa$variable==i,]$variable[1]
  df[i,2]<-wilc
  df<-unique(df)
}



##########################
aa$variable<-as.character(aa$variable)
merged<-merge(aa,df,by.x = "variable",by.y = "V1")
avoid<-unique(grep("Del",merged$variable,value = T))

merged$variable<-ifelse(merged$variable %in% avoid,NA,merged$variable)
merged<-merged[!is.na(merged$variable),]

###
#good<-unique(c(merged[merged$V2<0.05,]$variable))
gg<-ggplot(merged,aes(x=variable,y = value,fill=patient))+
   geom_boxplot()+facet_grid(~variable,scale = "free")+
  geom_text(data = merged,aes(label = ifelse(V2<0.05, "*", ""),y = 50),color="red",size=10)+
  theme_light()
  
 
gg
