rm(list=ls())
setwd("/home/ghepard/deeper_spike/29_09_2020/")

library(RColorBrewer)
library(wesanderson)
library(pheatmap)
library(gridExtra)
library(sqldf)
library(ape)
library(ggplot2)
library(gplots)
library(dplyr)
library(plyr)
library(ggrepel)
library(data.table)
# ###
dd<-read.csv("new_tabs/MutatedPosition.csv",sep=";")
#dd<-read.csv("new_tabs/MutatedPosition_min.csv",sep=";")
dd$X<-1:nrow(dd)
dd$count<-NULL
rownames(dd)<-dd$X
dd$X<-NULL
ss<-dd
tt<-data.frame(t(ss))
sel<-grep("BRASP",rownames(tt),value = TRUE)
sel2<-grep("BAL",rownames(tt),value = TRUE)
list<-c(sel,sel2)
tt$patient<-ifelse(rownames(tt) %in% list,"lower","upper")
pos<-colnames(tt)
###

df<-data.frame()
###
ll<-c()

c<-0

while (c<length(colnames(tt))){

for (i in pos){
    if (nrow(table(tt$patient,tt[,i]))==2 & ncol(table(tt$patient,tt[,i]))==2){
    f<-as.data.frame(table(tt$patient,tt[,i]))
    f$pos<-gsub("X","",i)
    #ll<-append(ll,f)
    write.csv(f,sprintf("dataframes/df_%i.csv",c),sep='\t')

    c<-c+1}
  }
}


###
file_list <- list.files("dataframes/")
setwd("dataframes/")
file_list<-list.files()
dataset <- ldply(file_list, read.csv, header=TRUE, sep=",")
###
dataset$X<-NULL
results<-dataset
#45 upper
#36 lower
results<-results[results$Var2!="2",]
results$pos<-as.numeric(results$pos)
results<-unique(results)
###

results<-results[results$Var2=="1",]

results$Freq<-ifelse(results$Var1=="lower",-results$Freq,results$Freq)

ggplot(data = results,aes(as.numeric(pos),Freq,group=Var2,fill=Var1))+
  geom_bar(aes(pos),stat = "identity",width = 3)+theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank())+
  scale_x_continuous(breaks = seq(from = 0, to = max(as.numeric(results$pos)), by = 100))



results[results$pos>0 & results$pos<10,]

# results$Var2<-as.character(results$Var2)
# ggplot(data = results,aes(as.numeric(pos),Freq,group=po,fill=Var2))+geom_bar(stat = "identity",width = 1)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   theme(axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),axis.title.x=element_blank())+
#   scale_x_continuous(breaks = seq(from = 0, to = max(as.numeric(results$pos)), by = 100))
