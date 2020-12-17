rm(list=ls())
setwd("/home/ghepard/deeper_spike/1_09_2020/")
library(sqldf)
library(RColorBrewer)
library(wesanderson)
library(pheatmap)
library(ape)
library(gplots)
library(ggplot2)
library(dplyr)
library(plyr)
library(data.table)
###
dd<-read.csv("new_tabs/MutatedCodon.csv",sep=";")
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
for (i in pos){
  if (nrow(table(tt$patient,tt[,i]))==2 & ncol(table(tt$patient,tt[,i]))==2){
    f<-fisher.test(table(tt$patient,tt[,i]))
    df[i,1]<-f$p.value
    #df[i,2]<-f$conf.int[1]
  }else{
    ll<-append(ll,i)     
  }
}

df$log<-abs(log10(df$V1))
df<-df[1:nrow(df)-1,]
ll<-as.data.frame(ll)
ll$V1<-1
rownames(ll)<-ll$ll
ll$ll<-NULL
df$log<-NULL
results<-rbind(df,ll)
rownames(results)<-gsub("X","",rownames(results))
results$log<-abs(log10(results$V1))
results$pos<-as.numeric(rownames(results))

results<-results[order(results$pos),]


##########################
#per la regione mettere etichette al limite, add column con posizione e regione
regions<-read.csv("new_tabs/regions_aa_coord.txt",sep='\t')

myColors <- brewer.pal(nrow(regions),"Set1")
regions$colors<-myColors
results$pos<-as.numeric(results$pos)

merged_res<-sqldf("select * from results f1 left join regions f2 
   on (f1.pos >= f2.start and f1.pos<= f2.stop)")

###
xx<-ifelse(merged_res$log>abs(log10(0.05)),as.numeric(merged_res$pos),NA)
yy<-ifelse(merged_res$log>abs(log10(0.05)),as.numeric(merged_res$log),NA)
###
lab<-ifelse(merged_res$log>abs(log10(0.05)),as.character(merged_res$pos),NA)
###
posizioni<-results[results$log>abs(log10(0.05)),]$pos
##########################
merged_res<-select(merged_res,regions,log,pos,start,stop,colors)
merged_res$pos<-as.character(merged_res$pos)
merged_res$middle<-(merged_res$start+merged_res$stop)/2

#merged_res$pos<-as.character(merged_res$pos)
##########################

df1<-data.frame()
for (p in posizioni){
  p<-paste("X",p,sep = "")
  aaa<-as.data.frame(table(tt$patient, tt[,p]))
  absence<-aaa[aaa$Var2==0,]
  pres<-aaa[aaa$Var2==1,]
  values<-c(abs(diff(absence$Freq)),abs(diff(pres$Freq)))
  if (values[1]<values[2]){
    df1[p,1]<-"presence"
    df1[p,2]<-pres[pres$Freq==max(pres$Freq),]$Var1
  }else{
    
    df1[p,1]<-"absence"
    df1[p,2]<-absence[absence$Freq==max(absence$Freq),]$Var1}
  
}


df1<-df1[,c(2,1)]
colnames(df1)<-c("V1","V2")
###################################
df1$pos<-rownames(df1)
df1$pos<-gsub("X","",df1$pos)
df1$pos<-as.numeric(df1$pos)
#upper absence è lower

df1[df1$V1=="upper" & df1$V2=="absence",]$V1<-"lower"
merged<-merge(df1,merged_res,by="pos",all.y = T)
lower<-df1[df1$V1=="lower",]$pos
upper<-df1[df1$V1=="upper",]$pos
merged$log<-ifelse(merged$pos %in% lower,-merged$log,merged$log)


###################################

merged$filled<-ifelse(merged$V1=="lower","#F8766D","#00BFC4")


ggplot(data = merged,aes(as.numeric(pos),log,fill=V1))+geom_bar(stat = "identity",position="dodge",width = 2)+
        geom_hline(yintercept=abs(log10(0.05)), linetype="dashed",color = "red", size=0.5)+
        geom_hline(yintercept=log10(0.05), linetype="dashed",color = "red", size=0.5)+
        scale_x_continuous(breaks = seq(from = 0, to = max(as.numeric(merged$pos)), by = 50))+
        geom_rect(aes(NULL, NULL, xmin = start, xmax = stop,ymin=6.5,ymax=8, fill = regions))+
        geom_label(aes(x = middle,y=7, label = regions),fill="white",colour="black",label.size = 0.4, label.padding = unit(0.2, "lines"))+
        theme_light()+theme(axis.text.x=element_text(angle = 45, hjust = 1))+
        theme(axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+
        theme(legend.position = "none")

res<-results[which(results$V1<0.05),]

site<-gsub("X","",rownames(res))
#write.csv(site,"codon.tab",row.names = F)
site
