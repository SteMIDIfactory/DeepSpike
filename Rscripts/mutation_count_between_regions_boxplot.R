rm(list=ls())
setwd("/home/ghepard/deeper_spike/14_10_2020/")
library(RColorBrewer)
library(wesanderson)
library(cluster)
library(ape)
library(gplots)
library(dplyr)
library(plyr)
library(data.table)
dd<-read.csv("new_tabs/MutationCount_Regions_new.csv",sep=";")
w<-data.frame("8076_TNAS",0,0,0,0,0,0,0,0,0,0)
colnames(w)<-c("X","NTD","RBD","SD1","SD2","FP","HR1","HR2","TM","S1","S2")
###

dd<-dd[,c("X","NTD","RBD","SD1","SD2","FP","HR1","HR2","TM","S1","S2")]
colnames(dd)

dd<-rbind.data.frame(dd,w)
rownames(dd)<-dd$X
dd$X<-NULL
sel<-grep("BRASP",rownames(dd),value = TRUE)
sel2<-grep("BAL",rownames(dd),value = TRUE)
list<-c(sel,sel2)
dd$patient<-ifelse(rownames(dd) %in% list,"Lower","Upper")
dd[is.na(dd)]<-0
aa<-melt(dd)

# ggplot(aa,aes(x=patient,y =value,fill=patient))+
#   geom_boxplot()+facet_grid(~variable,scale = "free")+theme_light()
aa$stat<-ifelse(aa$patient=="Lower",0,1)


df<-data.frame()
for (i in aa$variable){
  ss<-aa[aa$variable==i,]
  wilc<-wilcox.test(ss$value~ss$stat)$p.value
  #print (paste(aa[aa$variable==i,]$variable[1],wilc,sep=" "))
  df[i,1]<-aa[aa$variable==i,]$variable[1]
  df[i,2]<-wilc
  df<-unique(df)
}


##########################
aa$variable<-as.character(aa$variable)
merged<-merge(aa,df,by.x = "variable",by.y = "V1")
###
lista<-c("NTD","RBD","SD1","SD2","FP","HR1","HR2","TM","S1","S2")
library(gdata)
merged$variable <- reorder.factor(merged$variable, new.order=lista)

gg<-ggplot(merged,aes(x=variable,y = value,fill=patient))+
  geom_boxplot()+facet_grid(~variable,scale = "free")+
  geom_text(data = merged,aes(label = ifelse(V2<0.05, "*", ""),y = 120),color="red",size=10)+
  theme_light()+theme(axis.title.x = element_blank())+labs(y="Incidence")+
  guides(fill=guide_legend(title="District"))


gg

