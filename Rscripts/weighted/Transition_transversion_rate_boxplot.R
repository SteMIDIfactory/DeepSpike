rm(list=ls())
setwd("/home/ghepard/deeper_spike/14_10_2020/")
library(RColorBrewer)
library(exactRankTests)
library(wesanderson)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ape)
library(gplots)
library(dplyr)
library(data.table)
######################################## MUTATION PATT
dd_mut_pat<-read.csv("new_tabs/MutationPatterns.csv",sep=";",stringsAsFactors = F)
rownames(dd_mut_pat)<-dd_mut_pat$X
dd_mut_pat$X<-NULL
ff<-dd_mut_pat %>% replace(is.na(.),0)
tt<-data.frame(t(ff))
###
sel<-grep("BRASP",rownames(tt),value = TRUE)
sel2<-grep("BAL",rownames(tt),value = TRUE)
list<-c(sel,sel2)
tt$patient<-ifelse(rownames(tt) %in% list,"lower","upper")
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
dd_count$patient<-ifelse(rownames(dd_count) %in% list,"lower","upper")
##########
dd_count$patient_type<-ifelse(dd_count$patient=="lower",0,1)
dd_mut_pat[is.na(dd_mut_pat)==TRUE]<-0
###
dd_count$id<-rownames(dd_count)
dd_mut_pat<-as.data.frame(t(dd_mut_pat))
rownames(dd_mut_pat)<-gsub("X","",rownames(dd_mut_pat))
dd_mut_pat$id<-rownames(dd_mut_pat)
# 
#dd_mut_pat<-dd_mut_pat[rownames(dd_count),]
merged<-merge(dd_count,dd_mut_pat,by="id")
merged_weighted<-(merged[,c(5:26)]/merged$Total.mutations)
merged_weighted[is.na(merged_weighted)==TRUE]<-0
###
rownames(merged_weighted)<-merged$id
###
merged_weighted$patient_type<-ifelse(rownames(merged_weighted) %in% list,0,1)
ss_merged<-merged_weighted[,c('patient_type','A_C','A_Del','A_G','A_In','A_T','C_A','C_Del','C_G','C_In','C_T','Del_T','G_A','G_C','G_Del','G_In','G_T','In_G','T_A','T_C','T_Del','T_G','T_In')]
aa<-melt(ss_merged,id.vars = "patient_type")


df<-data.frame()
for (i in aa$variable){
  ss<-aa[aa$variable==i,]
  wilc<-wilcox.test(ss$value~as.factor(ss$patient_type))$p.value
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
avoid1<-unique(grep("Del",final$variable,value = T))
avoid2<-unique(grep("In",final$variable,value = T))
avoid<-c(avoid1,avoid2)

final$variable<-ifelse(final$variable %in% avoid,NA,final$variable)
final<-final[!is.na(final$variable),]


ggplot(final,aes(x=variable,y = value,fill=patient))+
   geom_boxplot()+facet_grid(~variable,scale = "free")+
   geom_text(data = final,aes(label = ifelse(V2<0.05, "*", ""),y = 0.75),color="red",size=10)+
   theme_light()+theme(axis.title.x=element_blank(),axis.text.x = element_blank())+
  guides(fill=guide_legend(title="District"))
                       

