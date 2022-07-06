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

weighted.var2 <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =na.rm)
}

weighted.var3 <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
}



setwd("/home/woot0002/DS_ind/")

load(file="Sanderson_EnsembleWeights_v4pronly_ann.Rdata")
BMAweights = read.table("best_BMA_combo.txt")
GCMhdat$BMA = t(BMAweights)[,1]

BMAweights = read.table("best_BMA_combo_LOCA.txt")
LOCAhdat$BMA = t(BMAweights)[,1]

GCMweights= GCMhdat
LOCAweights = LOCAhdat

var = varin = "pr"
type="ann"

if(var=="tmin"){
  varin="tasmin"
}
if(var=="tmax"){
  varin="tasmax"
}
if(var=="tmax95"){
  varin="tasmax"
}
if(var=="tmin32"){
  varin="tasmin"
}
if(var=="r1mm"){
  varin="pr"
}
if(var=="pr50"){
  varin="pr"
}
if(var=="rx1day"){
  varin="pr"
}
if(var=="rx5day"){
  varin="pr"
}

GCMfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)

GCMprojfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)
LOCAprojfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)

LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)

load("/home/woot0002/DS_ind/GCMlist.Rdata")

GCM_hfiles = GCM_pfiles = LOCA_hfiles = LOCA_pfiles = c()

for(i in 1:length(GCMlist)){
  
  GCM_hfiles[i] = GCMfiles[grep(paste(GCMlist[i],"_",sep=""),GCMfiles)]
  GCM_pfiles[i] = GCMprojfiles[grep(paste(GCMlist[i],"_",sep=""),GCMprojfiles)]
  
  LOCA_hfiles[i] = LOCAfiles[grep(paste(GCMlist[i],"_",sep=""),LOCAfiles)]
  LOCA_pfiles[i] = LOCAprojfiles[grep(paste(GCMlist[i],"_",sep=""),LOCAprojfiles)]
  
}

###
# create full filelist + metadata table - historical

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_hfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMhdat = filelist2[,c(2,3,4,6)]
names(GCMhdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_hfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCAhdat = filelist2[,c(2,3,4,6)]
names(LOCAhdat) = names(GCMhdat)

#All metadata
GCM = rep(NA,1)
exp = rep(NA,1)
DS = rep(NA,1)
training = "LIVNEH"
obsdat = data.frame(GCM,exp,DS,training)

GCMhdat = rbind(GCMhdat,obsdat)
LOCAhdat= rbind(LOCAhdat,obsdat)

# all files
GCMgroup = c(GCM_hfiles,LIVNEHfile)
LOCAgroup = c(LOCA_hfiles,LIVNEHfile)

###
# create full filelist + metadata table - projected

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_pfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMpdat = filelist2[,c(2,3,4,6)]
names(GCMpdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_pfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCApdat = filelist2[,c(2,3,4,6)]
names(LOCApdat) = names(GCMpdat)

# all files
GCMpgroup = GCM_pfiles
LOCApgroup = LOCA_pfiles

######
# state mask

statelist = c("oklahoma","texas","louisiana","new mexico")
regionmasklist = list()
for(s in 1:length(statelist)){
  test = nc_open(paste("/home/woot0002/DS_ind/",statelist[s],"_mask.nc",sep=""))
  regionmasklist[[s]] = ncvar_get(test,"mask")
  nc_close(test)
}


######
# Gather data

ncvarname = paste(var,"climo",sep="")

### GCM hist + Livneh
GCMhvarmat = matrix(NA,nrow=length(GCMgroup),ncol=length(statelist)+1)
for(i in 1:length(GCMgroup)){
  nctest = nc_open(GCMgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      tmpmat = apply(tmp,c(1,2),sum,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
    if(var=="tmax" | var=="tmin"){
      tmpmat = apply(tmp,c(1,2),mean,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
  for(s in 1:5){
    if(s<=4){
      regionidx = which(regionmasklist[[s]]==1)
      GCMhvarmat[i,s] = mean(tmpmat[regionidx],na.rm=TRUE)
    } else {
      GCMhvarmat[i,s] = mean(tmpmat,na.rm=TRUE)
    }
  }
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon");
  nc_close(nctest)
}


### GCM projected change
GCMpvarmat = matrix(NA,nrow=length(GCMpgroup),ncol=length(statelist)+1)
for(i in 1:length(GCMpgroup)){
  nctest = nc_open(GCMpgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      tmpmat = apply(tmp,c(1,2),sum,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
    if(var=="tmax" | var=="tmin"){
      tmpmat = apply(tmp,c(1,2),mean,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
for(s in 1:5){
  if(s<=4){
    regionidx = which(regionmasklist[[s]]==1)
    GCMpvarmat[i,s] = mean(tmpmat[regionidx],na.rm=TRUE)
  } else {
    GCMpvarmat[i,s] = mean(tmpmat,na.rm=TRUE)
  }
}
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

### LOCA historical + Livneh

LOCAhvarmat = matrix(NA,nrow=length(LOCAgroup),ncol=length(statelist)+1)
for(i in 1:length(LOCAgroup)){
  nctest = nc_open(LOCAgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)

    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      tmpmat = apply(tmp,c(1,2),sum,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
    if(var=="tmax" | var=="tmin"){
      tmpmat = apply(tmp,c(1,2),mean,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
  for(s in 1:5){
    if(s<=4){
      regionidx = which(regionmasklist[[s]]==1)
      LOCAhvarmat[i,s] = mean(tmpmat[regionidx],na.rm=TRUE)
    } else {
      LOCAhvarmat[i,s] = mean(tmpmat,na.rm=TRUE)
    }
  }
  
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}


### LOCA projected change

LOCApvarmat = matrix(NA,nrow=length(LOCApgroup),ncol=length(statelist)+1)
for(i in 1:length(LOCApgroup)){
  nctest = nc_open(LOCApgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      tmpmat = apply(tmp,c(1,2),sum,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
    if(var=="tmax" | var=="tmin"){
      tmpmat = apply(tmp,c(1,2),mean,na.rm=TRUE)
      tmpmat = ifelse(is.na(tmp[,,1])==FALSE,tmpmat,NA)
    } 
  for(s in 1:5){
    if(s<=4){
      regionidx = which(regionmasklist[[s]]==1)
      LOCApvarmat[i,s] = mean(tmpmat[regionidx],na.rm=TRUE)
    } else {
      LOCApvarmat[i,s] = mean(tmpmat,na.rm=TRUE)
    }
  }
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

#######
# projected changes

GCMchange = GCMpvarmat-GCMhvarmat[1:26]
LOCAchange = LOCApvarmat-LOCAhvarmat[1:26]
OBS = LOCAhvarmat[27,]

######
# prep weights

GCMweights$Wh = (GCMweights$Wuh*GCMweights$Wqh)/sum(GCMweights$Wuh*GCMweights$Wqh)
GCMweights$Wc = (GCMweights$Wuc*GCMweights$Wqc)/sum(GCMweights$Wuc*GCMweights$Wqc)
GCMweights$Ws = GCMweights$Wqh/sum(GCMweights$Wqh)

LOCAweights$Wh = (LOCAweights$Wuh*LOCAweights$Wqh)/sum(LOCAweights$Wuh*LOCAweights$Wqh)
LOCAweights$Wc = (LOCAweights$Wuc*LOCAweights$Wqc)/sum(LOCAweights$Wuc*LOCAweights$Wqc)
LOCAweights$Ws = LOCAweights$Wqh/sum(LOCAweights$Wqh)


######
# Calculate historical means (weighted and unweighted)

GCMhistmeanmat = LOCAhistmeanmat = matrix(NA,nrow=5,ncol=5)

for(i in 1:5){
  if(i==1){
    GCMhistmeanmat[i,] = apply(GCMhvarmat[1:26,],2,mean,na.rm=TRUE)
    LOCAhistmeanmat[i,] = apply(LOCAhvarmat[1:26,],2,mean,na.rm=TRUE)
  } else {
    if(i==2){
      GCMw = GCMweights$Ws
      LOCAw = LOCAweights$Ws
    }
    if(i==3){
      GCMw = GCMweights$Wh
      LOCAw = LOCAweights$Wh
    }
    if(i==4){
      GCMw = GCMweights$Wc
      LOCAw = LOCAweights$Wc
    }
    if(i==5){
      GCMw = GCMweights$BMA
      LOCAw = LOCAweights$BMA
    }
    for(s in 1:5){
      GCMhistmeanmat[i,s] = sum(GCMhvarmat[1:26,s]*GCMw)
      LOCAhistmeanmat[i,s] = sum(LOCAhvarmat[1:26,s]*LOCAw)
    }
  }
}

######
# Calculate change means (weighted and unweighted)

GCMchangemeanmat = LOCAchangemeanmat = matrix(NA,nrow=5,ncol=5)

for(i in 1:5){
  if(i==1){
    GCMchangemeanmat[i,] = apply(GCMchange,2,mean,na.rm=TRUE)
    LOCAchangemeanmat[i,] = apply(LOCAchange,2,mean,na.rm=TRUE)
  } else {
    if(i==2){
      GCMw = GCMweights$Ws
      LOCAw = LOCAweights$Ws
    }
    if(i==3){
      GCMw = GCMweights$Wh
      LOCAw = LOCAweights$Wh
    }
    if(i==4){
      GCMw = GCMweights$Wc
      LOCAw = LOCAweights$Wc
    }
    if(i==5){
      GCMw = GCMweights$BMA
      LOCAw = LOCAweights$BMA
    }
    for(s in 1:5){
      GCMchangemeanmat[i,s] = sum(GCMchange[1:26,s]*GCMw)
      LOCAchangemeanmat[i,s] = sum(LOCAchange[1:26,s]*LOCAw)
    }
  }
}

######
# Calculate change variance (weighted and unweighted)

GCMchangevarmat = LOCAchangevarmat = matrix(NA,nrow=5,ncol=5)

for(i in 1:5){
  if(i==1){
    for(s in 1:5){
      GCMchangevarmat[i,s] = var(GCMchange[,s],na.rm=TRUE)
      LOCAchangevarmat[i,s] = var(LOCAchange[,s],na.rm=TRUE)
    }
  } else {
    if(i==2){
      GCMw = GCMweights$Ws
      LOCAw = LOCAweights$Ws
    }
    if(i==3){
      GCMw = GCMweights$Wh
      LOCAw = LOCAweights$Wh
    }
    if(i==4){
      GCMw = GCMweights$Wc
      LOCAw = LOCAweights$Wc
    }
    if(i==5){
      GCMw = GCMweights$BMA
      LOCAw = LOCAweights$BMA
    }
    for(s in 1:5){
      GCMchangevarmat[i,s] = weighted.var(GCMchange[1:26,s],w=GCMw)
      LOCAchangevarmat[i,s] = weighted.var(LOCAchange[1:26,s],w=LOCAw)
    }
  }
}

######
######
# BMA posterior weights calcs

GCMBMAweights = read.table("posterior_BMA_combo.txt")
GCMBMAweightst = t(GCMBMAweights)

LOCABMAweights = read.table("posterior_BMA_combo_LOCA.txt")
LOCABMAweightst = t(LOCABMAweights)

GCMhistmeanmatBMA = LOCAhistmeanmatBMA = matrix(NA,nrow=100,ncol=5)

for(i in 1:100){
  GCMw = GCMBMAweightst[,i]
  LOCAw = LOCABMAweightst[,i]
  for(s in 1:5){
    GCMhistmeanmatBMA[i,s] = sum(GCMhvarmat[1:26,s]*GCMw)
    LOCAhistmeanmatBMA[i,s] = sum(LOCAhvarmat[1:26,s]*LOCAw)
  }
}


######
# Calculate change means (weighted and unweighted)

GCMchangemeanmatBMA = LOCAchangemeanmatBMA = matrix(NA,nrow=100,ncol=5)

for(i in 1:100){
  GCMw = GCMBMAweightst[,i]
  LOCAw = LOCABMAweightst[,i]
  for(s in 1:5){
    GCMchangemeanmatBMA[i,s] = sum(GCMchange[1:26,s]*GCMw)
    LOCAchangemeanmatBMA[i,s] = sum(LOCAchange[1:26,s]*LOCAw)
  }
}

######
# Calculate change variance (weighted and unweighted)

GCMchangevarmatBMA = LOCAchangevarmatBMA = matrix(NA,nrow=100,ncol=5)

for(i in 1:100){
  GCMw = GCMBMAweightst[,i]
  LOCAw = LOCABMAweightst[,i]
    for(s in 1:5){
      GCMchangevarmatBMA[i,s] = weighted.var(GCMchange[1:26,s],w=GCMw)
      LOCAchangevarmatBMA[i,s] = weighted.var(LOCAchange[1:26,s],w=LOCAw)
    }
  }

#####

GCMhistmeanmatBMA = data.frame(GCMhistmeanmatBMA)
names(GCMhistmeanmatBMA) = c(statelist,"full")
GCMhistmeanmatBMA$DS = "CMIP5"

LOCAhistmeanmatBMA = data.frame(LOCAhistmeanmatBMA)
names(LOCAhistmeanmatBMA) = c(statelist,"full")
LOCAhistmeanmatBMA$DS = "LOCA"

histmeanBMA = rbind(GCMhistmeanmatBMA,LOCAhistmeanmatBMA)


GCMchangemeanmatBMA = data.frame(GCMchangemeanmatBMA)
names(GCMchangemeanmatBMA) = c(statelist,"full")
GCMchangemeanmatBMA$DS = "CMIP5"

LOCAchangemeanmatBMA = data.frame(LOCAchangemeanmatBMA)
names(LOCAchangemeanmatBMA) = c(statelist,"full")
LOCAchangemeanmatBMA$DS = "LOCA"

changemeanBMA = rbind(GCMchangemeanmatBMA,LOCAchangemeanmatBMA)


GCMchangevarmatBMA = data.frame(GCMchangevarmatBMA)
names(GCMchangevarmatBMA) = c(statelist,"full")
GCMchangevarmatBMA$DS = "CMIP5"

LOCAchangevarmatBMA = data.frame(LOCAchangevarmatBMA)
names(LOCAchangevarmatBMA) = c(statelist,"full")
LOCAchangevarmatBMA$DS = "LOCA"

changevarBMA = rbind(GCMchangevarmatBMA,LOCAchangevarmatBMA)

#####

histmeanBMAbias = histmeanBMA
histmeanBMAbias[,1:5] = histmeanBMAbias[,1:5]-OBS

#####

BMAtable = NULL

for(i in 1:5){
  
  statename = names(histmeanBMA)[i]
  
  tmp1 = histmeanBMA[,c(i,6)]
  names(tmp1)[1] = "histmean"
  tmp2 = histmeanBMAbias[,c(i,6)]
  names(tmp2)[1] = "histmeanbias"
  tmp3 = changemeanBMA[,c(i,6)]
  names(tmp3)[1] = "changemean"
  tmp4 = changevarBMA[,c(i,6)]
  names(tmp4)[1] = "changevar"
  
  tmp = cbind(tmp1,tmp2[,1])
  tmp = cbind(tmp,tmp3[,1])
  tmp = cbind(tmp,tmp4[,1])
  
  names(tmp) = c("histmean","DS","histmeanbias","changemean","changevar")
  
  tmp$state = statename
  
  BMAtable = rbind(BMAtable,tmp)
  
}

######

BMAtable$state <- factor(BMAtable$state,levels = c('full','new mexico','texas','oklahoma','louisiana'),ordered = TRUE)

######

GCMhistmeanbias = GCMhistmeanmat
LOCAhistmeanbias = LOCAhistmeanmat
for(i in 1:5){
  GCMhistmeanbias[i,]=GCMhistmeanmat[i,]-OBS
  LOCAhistmeanbias[i,]=LOCAhistmeanmat[i,]-OBS
}


BMAsmalltable = NULL

group = c("unweighted","skill","SI-h","SI-c","BMA")

for(i in 1:5){
  GCMtmp = data.frame(GCMhistmeanmat[i,],GCMhistmeanbias[i,],GCMchangemeanmat[i,],GCMchangevarmat[i,])
  LOCAtmp = data.frame(LOCAhistmeanmat[i,],LOCAhistmeanbias[i,],LOCAchangemeanmat[i,],LOCAchangevarmat[i,])
  names(GCMtmp) = names(LOCAtmp) = c("histmean","histmeanbias","changemean","changevar")
  GCMtmp$DS = "CMIP5"
  LOCAtmp$DS = "LOCA"
  GCMtmp$group = rep(group[i],5)
  LOCAtmp$group = rep(group[i],5)
  GCMtmp$state = c(statelist,"full")
  LOCAtmp$state = c(statelist,"full")
  tmp = rbind(GCMtmp,LOCAtmp)
  BMAsmalltable = rbind(BMAsmalltable,tmp)
}


ggplot(BMAtable, aes(x=state, y=histmean, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=histmean,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("cyan", "yellow"))+
  ggtitle("Hist Ens. Mean")+xlab("")+ylab("Mean (mm)")+ylim(range(BMAtable$histmean,na.rm=TRUE)) 

ggplot(BMAtable, aes(x=state, y=histmeanbias, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=histmeanbias,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("cyan", "yellow"))+
  ggtitle("Hist Ens. Mean Bias")+xlab("")+ylab("Bias (mm)")+ylim(range(BMAtable$histmeanbias,na.rm=TRUE)) 

ggplot(BMAtable, aes(x=state, y=changemean, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=changemean,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("cyan", "yellow"))+
  ggtitle("Ens. Mean Projected Change")+xlab("")+ylab("Change (mm)")+ylim(range(BMAtable$changemean,na.rm=TRUE)) 

ggplot(BMAtable, aes(x=state, y=changevar, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=changevar,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("cyan", "yellow"))+
  ggtitle("Ens. Variance. Projected Change")+xlab("")+ylab("Variance (mm^2)")+ylim(range(BMAtable$changevar,na.rm=TRUE)) 

##########

ggplot(BMAtable, aes(x=state, y=histmean)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=histmean,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean")+xlab("")+ylab("Mean (mm)")+ylim(range(BMAtable$histmean,na.rm=TRUE)) +facet_grid(cols=vars(DS))

ggplot(BMAtable, aes(x=state, y=histmeanbias)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=histmeanbias,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean Bias")+xlab("")+ylab("Bias (mm)")+ylim(range(BMAtable$histmeanbias,na.rm=TRUE)) +facet_grid(cols=vars(DS))+geom_hline(yintercept=0,linetype="dashed")

ggplot(BMAtable, aes(x=state, y=changemean)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=changemean,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Ens. Mean Projected Change")+xlab("")+ylab("Change (mm)")+ylim(range(BMAtable$changemean,na.rm=TRUE)) +facet_grid(cols=vars(DS))+geom_hline(yintercept=0,linetype="dashed")

ggplot(BMAtable, aes(x=state, y=changevar)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=BMAsmalltable, mapping=aes(x=state,y=changevar,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Ens.Variance Projected Change")+xlab("")+ylab("Variance (mm^2)")+ylim(range(BMAtable$changevar,na.rm=TRUE)) +facet_grid(cols=vars(DS))+geom_hline(yintercept=0,linetype="dashed")

