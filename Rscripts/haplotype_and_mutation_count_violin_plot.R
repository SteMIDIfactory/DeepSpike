rm(list=ls())
setwd("/home/ghepard/deeper_spike/14_10_2020/")
library(ggfortify)
library(lmtest)
library(exactRankTests)
library(RColorBrewer)
library(wesanderson)
library(cluster)
library(pheatmap)
library(ape)
library(pROC)
library(MASS)
library(InformationValue)
library(gplots)
library(dplyr)
library(plyr)
library(data.table)

dd<-read.csv("new_tabs/MutationCount.csv",sep=";")
w<-data.frame("8076_TNAS",0)
colnames(w)<-c("X","Total.mutations")
dd<-rbind.data.frame(dd,w)
rownames(dd)<-dd$X
dd$X<-NULL

#tt<-data.frame(t(dd))
sel<-grep("BRASP",rownames(dd),value = TRUE)
sel2<-grep("BAL",rownames(dd),value = TRUE)
list<-c(sel,sel2)
dd$patient<-ifelse(rownames(dd) %in% list,"lower","upper")

##########
dd$stat<-ifelse(dd$patient=="lower",0,1)
###
wilc<-wilcox.test(dd$Total.mutations~dd$stat)$p.value



#########
p<-wilcox.exact(dd$Total.mutations~as.factor(dd$stat))
pval<-format.pval(p$p.value, digits = max(2, getOption("digits") - 2))

gg<-ggplot(dd,aes(x=patient,y = Total.mutations,fill=patient))+
  geom_violin(draw_quantiles = T)+theme_light()+
  geom_boxplot(width=0.1,fill = "white")+labs(x=NULL,y="Total mutations")+
  guides(fill=guide_legend(title="District"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
        


                                              
###################################

dd<-read.csv("new_tabs/Number_of_haplotypes.csv",sep=";")
rownames(dd)<-dd$X
dd$X<-NULL

sel<-grep("BRASP",rownames(dd),value = TRUE)
sel2<-grep("BAL",rownames(dd),value = TRUE)
list<-c(sel,sel2)
dd$patient<-ifelse(rownames(dd) %in% list,"lower","upper")
dd$stat<-ifelse(dd$patient=="lower",0,1)

###
wilc<-wilcox.test(dd$Number.of.haplotypes~as.factor(dd$stat))$p.value

#########
#p<-t.test(dd$Number.of.haplotypes~as.factor(dd$stat))
pval<-format.pval(wilc, digits = max(2, getOption("digits") - 2))


gg1<-ggplot(dd,aes(x=patient,y = Number.of.haplotypes,fill=patient))+
  geom_violin(draw_quantiles = T)+theme_light()+
  geom_boxplot(width=0.1,fill = "white")+
  labs(x=NULL,y="Number of haplotypes")+theme(legend.position = "none",
                                              axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(),
                                                    axis.ticks.x=element_blank())
                                              

library(gridExtra)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(gg)

getwd()
p3 <- grid.arrange(arrangeGrob(gg + theme(legend.position="none",
                                          axis.title.y = element_text(size=15),
                                          axis.text.y = element_text(size=12)),
                               gg1 + theme(legend.position="none",
                              axis.title.y = element_text(size=15),
                              axis.text.y = element_text(size=12)),
                               nrow=2),mylegend, ncol=3,widths=c(2.3, 0.3, 0.2))
                   
