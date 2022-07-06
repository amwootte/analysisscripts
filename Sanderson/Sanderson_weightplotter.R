source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(mapdata)
library(maptools)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(ggplot2)
library(modi)

setwd("/home/woot0002/DS_ind/")
type="ann"

vars=c("tmax","pr")
WUs = c("full","louisiana","new mexico")

GCMweights_all = LOCAweights_all = BMAweightsGCM_all = BMAweightsLOCA_all = NULL


for(v in 1:length(vars)){
  for(w in 1:length(WUs)){
    var = varin = vars[v]
    weightingused = WUs[w]
    
    if(weightingused=="full"){
      load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,".Rdata",sep=""))
      BMAweights_GCM = read.table(paste("best_BMA_combo_",var,".txt",sep=""))
      BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,".txt",sep=""))
      BMAweightsGCM = read.table(paste("posterior_BMA_combo_",var,".txt",sep=""))
      BMAweightsLOCA = read.table(paste("posterior_BMA_combo_LOCA_",var,".txt",sep=""))
    } else {
      load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,"_",weightingused,".Rdata",sep=""))
      BMAweights_GCM = read.table(paste("best_BMA_combo_",var,"_",weightingused,".txt",sep=""))
      BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
      BMAweightsGCM = read.table(paste("posterior_BMA_combo_",var,"_",weightingused,".txt",sep=""))
      BMAweightsLOCA = read.table(paste("posterior_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
    }
    
    GCMweights = GCMhdat
    LOCAweights = LOCAhdat
    
    GCMweights$Wh = (GCMweights$Wuh*GCMweights$Wqh)/sum(GCMweights$Wuh*GCMweights$Wqh)
    GCMweights$Wc = (GCMweights$Wuc*GCMweights$Wqc)/sum(GCMweights$Wuc*GCMweights$Wqc)
    GCMweights$Ws = GCMweights$Wqh/sum(GCMweights$Wqh)
    
    LOCAweights$Wh = (LOCAweights$Wuh*LOCAweights$Wqh)/sum(LOCAweights$Wuh*LOCAweights$Wqh)
    LOCAweights$Wc = (LOCAweights$Wuc*LOCAweights$Wqc)/sum(LOCAweights$Wuc*LOCAweights$Wqc)
    LOCAweights$Ws = LOCAweights$Wqh/sum(LOCAweights$Wqh)
    
    GCMweights$BMA = t(BMAweights_GCM)[,1]
    LOCAweights$BMA = t(BMAweights_LOCA)[,1]
    
    GCMweights$region = weightingused
    LOCAweights$region = weightingused
    
    GCMweights$var = varin
    LOCAweights$var = varin
    
    GCMweights_all = rbind(GCMweights_all,GCMweights)
    LOCAweights_all = rbind(LOCAweights_all,LOCAweights)
    
    GCMbase = GCMweights[,c(1,3,4,5,14,15)]
    LOCAbase = LOCAweights[,c(1,3,4,5,14,15)]
    
    BMAweightsGCM_tmp = NULL
    BMAweightsLOCA_tmp = NULL
    
    for(i in 1:100){
      BMA = t(BMAweightsGCM[i,])
      tmp = data.frame(GCMbase,BMA)
      names(tmp)[ncol(tmp)] = "BMA"
      tmp$BMAit = i
      BMAweightsGCM_tmp = rbind(BMAweightsGCM_tmp,tmp)
      
      BMA = t(BMAweightsLOCA[i,])
      tmp = data.frame(LOCAbase,BMA)
      names(tmp)[ncol(tmp)] = "BMA"
      tmp$BMAit = i
      BMAweightsLOCA_tmp = rbind(BMAweightsLOCA_tmp,tmp)
    }
    
    BMAweightsGCM_all = rbind(BMAweightsGCM_all,BMAweightsGCM_tmp)
    BMAweightsLOCA_all = rbind(BMAweightsLOCA_all,BMAweightsLOCA_tmp)
    
    message("Finished combining data on ",w," / 3 weight domains and ",v," / 2 weight variables")
  }
}

GCMwa1 = GCMweights_all[,c(1:5,14,15,12)]
GCMwa2 = GCMweights_all[,c(1:5,14,15,10)]
GCMwa3 = GCMweights_all[,c(1:5,14,15,11)]
GCMwa4 = GCMweights_all[,c(1:5,14,15,13)]

LOCAwa1 = LOCAweights_all[,c(1:5,14,15,12)]
LOCAwa2 = LOCAweights_all[,c(1:5,14,15,10)]
LOCAwa3 = LOCAweights_all[,c(1:5,14,15,11)]
LOCAwa4 = LOCAweights_all[,c(1:5,14,15,13)]

names(GCMwa1)[ncol(GCMwa1)] = names(GCMwa2)[ncol(GCMwa2)] = names(GCMwa3)[ncol(GCMwa3)] = names(GCMwa4)[ncol(GCMwa4)] = "W"
names(LOCAwa1)[ncol(LOCAwa1)] = names(LOCAwa2)[ncol(LOCAwa2)] = names(LOCAwa3)[ncol(LOCAwa3)] = names(LOCAwa4)[ncol(LOCAwa4)] = "W"

GCMwa1$weightscheme = "Skill"
GCMwa2$weightscheme = "SI-h"
GCMwa3$weightscheme = "SI-c"
GCMwa4$weightscheme = "BMA best"

LOCAwa1$weightscheme = "Skill"
LOCAwa2$weightscheme = "SI-h"
LOCAwa3$weightscheme = "SI-c"
LOCAwa4$weightscheme = "BMA best"

GCMwa = rbind(GCMwa1,GCMwa2)
GCMwa = rbind(GCMwa,GCMwa3)
GCMwa = rbind(GCMwa,GCMwa4)

LOCAwa = rbind(LOCAwa1,LOCAwa2)
LOCAwa = rbind(LOCAwa,LOCAwa3)
LOCAwa = rbind(LOCAwa,LOCAwa4)

names(BMAweightsGCM_all)[7]="W"
names(BMAweightsLOCA_all)[7]="W"


ggplot(BMAweightsGCM_all, aes(x=GCM, y=W)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="grey40", outlier.shape=8,outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=GCMwa, mapping=aes(x=GCM,y=W,colour = factor(weightscheme)),size=2) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ylim(0,1) +
  ggtitle("Weight to each GCM - CMIP5")+xlab("")+ylab("Weight") +facet_grid(cols=vars(var),rows=vars(region))
ggsave("/home/woot0002/DS_ind/GCMweightsv2.pdf",device = "pdf",width=10,height=8)

ggplot(BMAweightsLOCA_all, aes(x=GCM, y=W)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="grey40", outlier.shape=8,outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=LOCAwa, mapping=aes(x=GCM,y=W,colour = factor(weightscheme)),size=2) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ylim(0,1) +
  ggtitle("Weight to each GCM - LOCA")+xlab("")+ylab("Weight") +facet_grid(cols=vars(var),rows=vars(region))
ggsave("/home/woot0002/DS_ind/LOCAweightsv2.pdf",device = "pdf",width=10,height=8)


####
# plot without BMA

GCMwa_1o = subset(GCMwa,weightscheme!="BMA best")
LOCAwa_1o = subset(LOCAwa,weightscheme!="BMA best")

ggplot(GCMwa_1o, aes(x=GCM, y=W)) + 
  geom_point(mapping=aes(x=GCM,y=W,colour = factor(weightscheme)),size=1.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Weight to each GCM - CMIP5")+xlab("")+ylab("Weight") +facet_grid(cols=vars(var),rows=vars(region))

ggplot(LOCAwa_1o, aes(x=GCM, y=W)) + 
  geom_point(mapping=aes(x=GCM,y=W,colour = factor(weightscheme)),size=1.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Weight to each GCM - LOCA")+xlab("")+ylab("Weight") +facet_grid(cols=vars(var),rows=vars(region))

ggplot(GCMwa, aes(x=GCM, y=W)) + 
  geom_point(mapping=aes(x=GCM,y=W,colour = factor(weightscheme)),size=1.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Weight to each GCM - CMIP5")+xlab("")+ylab("Weight") +facet_grid(cols=vars(var),rows=vars(region))

ggplot(LOCAwa, aes(x=GCM, y=W)) + 
  geom_point(mapping=aes(x=GCM,y=W,colour = factor(weightscheme)),size=1.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Weight to each GCM - LOCA")+xlab("")+ylab("Weight") +facet_grid(cols=vars(var),rows=vars(region))


