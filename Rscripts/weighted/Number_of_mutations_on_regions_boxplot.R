rm(list=ls())
setwd("/home/ghepard/deeper_spike/14_10_2020/")
library(RColorBrewer)
library(wesanderson)
library(pheatmap)
library(reshape2)
library(ape)
library(exactRankTests)
library(ggplot2)
library(dplyr)
library(data.table)

dd_count<-read.csv("new_tabs/MutationCount.csv",sep=";")
w<-data.frame("8076_TNAS",0)
colnames(w)<-c("X","Total.mutations")
dd_count<-rbind.data.frame(dd_count,w)
###
rownames(dd_count)<-dd_count$X
dd_count$X<-NULL
sel<-grep("BRASP",rownames(dd_count),value = TRUE)
sel2<-grep("BAL",rownames(dd_count),value = TRUE)
list<-c(sel,sel2)
# ##########
# dd_synnotsyn[is.na(dd_synnotsyn)==TRUE]<-0
###
#############################################
dd_region<-read.csv("new_tabs/MutationCount_Regions_new.csv",sep=";")
rownames(dd_region)<-dd_region$X
dd_region$X<-NULL
sel<-grep("BRASP",rownames(dd_region),value = TRUE)
sel2<-grep("BAL",rownames(dd_region),value = TRUE)
list<-c(sel,sel2)
dd_region$patient<-ifelse(rownames(dd_region) %in% list,"lower","upper")
dd_region[is.na(dd_region)]<-0
#####################################################
dd_count$id<-rownames(dd_count)
dd_region$id<-rownames(dd_region)
####################################################
merged<-merge(dd_count,dd_region,by="id")

#merged_weighted[,c('patient_type',"NTD","RBD","SD1","SD2","FP","HR1","HR2","TM","S1","S2")
head(merged)
merged_weighted<-(merged[,c(3:12)]/merged$Total.mutations)
###

rownames(merged_weighted)<-merged$id
merged_weighted$patient_type<-ifelse(rownames(merged_weighted) %in% list,0,1)
colnames(merged_weighted)

ss_merged<-merged_weighted[,c('patient_type',"NTD","RBD","SD1","SD2","FP","HR1","HR2","TM","S1","S2")]


aa<-melt(ss_merged,id.vars = "patient_type")

tail(aa$variable)


####################################################

df<-data.frame()
for (i in aa$variable){
  ss<-aa[aa$variable==i,]
  wilc<-wilcox.test(ss$value~ss$patient_type)$p.value
  df[i,1]<-aa[aa$variable==i,]$variable[1]
  df[i,2]<-wilc
  df<-unique(df)
}

###

aa$variable<-as.character(aa$variable)
df$V1<-as.character(df$V1)
final<-merge(aa,df,by.x = "variable",by.y = "V1")
final$patient<-ifelse(final$patient==0,"Lower","Upper")
colnames(final)[1]<-"regions"
###
lista<-c("NTD","RBD","SD1","SD2","FP","HR1","HR2","TM","S1","S2")
library(gdata)
final$regions <- reorder.factor(final$regions, new.order=lista)

ggplot(final,aes(x=regions,y = value,fill=patient))+
  geom_boxplot()+facet_grid(~regions,scale = "free")+
  geom_text(data = final,aes(label = ifelse(V2<0.05, "*", ""),y = 1.1),color="red",size=10)+
  theme_light()+theme(axis.title.x=element_blank(),axis.text.x = element_blank())+
  guides(fill=guide_legend(title="District"))+labs(y="Weighted incidence")


df[df$V2<0.05,]
