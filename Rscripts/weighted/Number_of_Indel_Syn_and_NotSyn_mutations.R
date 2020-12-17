rm(list=ls())
setwd("/home/ghepard/deeper_spike/14_10_2020//")
library(RColorBrewer)
library(wesanderson)
library(pheatmap)
library(reshape2)
library(ape)
library(ggplot2)
library(dplyr)
#library(plyr)
library(data.table)
######################################## MUTATION PATT
dd_synnotsyn<-read.csv("new_tabs/CountSynNotSyn.csv",sep=";")
w<-data.frame("8076_TNAS",0,0,0,0,0)
colnames(w)<-c("X","NotSyn","Syn","Indel","STOP","NA.")
dd_synnotsyn<-rbind.data.frame(dd_synnotsyn,w)
dd_synnotsyn[is.na(dd_synnotsyn)]<-0
dd_synnotsyn<-dd_synnotsyn[-6]
rownames(dd_synnotsyn)<-dd_synnotsyn$X
dd_synnotsyn$X<-NULL
sel<-grep("BRASP",rownames(dd_synnotsyn),value = TRUE)
sel2<-grep("BAL",rownames(dd_synnotsyn),value = TRUE)
list<-c(sel,sel2)
dd_synnotsyn$patient<-ifelse(rownames(dd_synnotsyn) %in% list,"lower","upper")



################################################ MUTATION COUNT
dd_count<-read.csv("new_tabs/MutationCount.csv",sep=";")
w<-data.frame("8076_TNAS",0)
colnames(w)<-c("X","Total.mutations")
dd_count<-rbind.data.frame(dd_count,w)
rownames(dd_count)<-dd_count$X
dd_count$X<-NULL
sel<-grep("BRASP",rownames(dd_count),value = TRUE)
sel2<-grep("BAL",rownames(dd_count),value = TRUE)
list<-c(sel,sel2)
# dd_count$patient<-ifelse(rownames(dd_count) %in% list,"lower","upper")
# ##########
# dd_count$patient_type<-ifelse(dd_count$patient=="lower",0,1)
# dd_synnotsyn[is.na(dd_synnotsyn)==TRUE]<-0
###
####
dd_count$id<-rownames(dd_count)
dd_synnotsyn$id<-rownames(dd_synnotsyn)

merged<-merge(dd_count,dd_synnotsyn,by="id")
merged_weighted<-(merged[,c(3:6)]/merged$Total.mutations)
###
merged_weighted[is.na(merged_weighted)==TRUE]<-0
rownames(merged_weighted)<-merged$id
###
merged_weighted$patient_type<-ifelse(rownames(merged_weighted) %in% list,0,1)
ss_merged<-merged_weighted[,c('patient_type',"NotSyn","Syn","Indel","STOP")]
aa<-melt(ss_merged,id.vars = "patient_type")

#cat(shQuote(colnames(merged[,c(5:27)])),sep = ",")

df<-data.frame()
for (i in aa$variable){
  ss<-aa[aa$variable==i,]
  wilc<-wilcox.test(ss$value~ss$patient_type)$p.value
  
  #print (paste(aa[aa$variable==i,]$variable[1],wilc,sep=" "))
  df[i,1]<-aa[aa$variable==i,]$variable[1]
  df[i,2]<-wilc
  df<-unique(df)
}


###
aa$variable<-as.character(aa$variable)
df$V1<-as.character(df$V1)
final<-merge(aa,df,by.x = "variable",by.y = "V1")
final$patient<-ifelse(final$patient==0,"lower","upper")
###
final<-final[final$variable!="STOP",]
ggplot(final,aes(x=variable,y = value,fill=patient))+
  geom_boxplot()+facet_grid(~variable,scale = "free")+
  geom_text(data = final,aes(label = ifelse(V2<0.05, "*", ""),y = 2),color="red",size=10)+
  theme_light()+theme(axis.title.x=element_blank())

