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
# Gather data

ncvarname = paste(var,"climo",sep="")

### GCM hist + Livneh
GCMhvardatalist = list()
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
  if(type!="ann"){
    GCMhvardatalist[[i]] = tmp
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      GCMhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      GCMhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      GCMhvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      GCMhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist[[i]],NA)
    } 
  }
  
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon");
  nc_close(nctest)
}

sapply(GCMhvardatalist,mean,na.rm=TRUE)

### GCM projected change
GCMpvardatalist = list()
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
  if(type!="ann"){
    GCMpvardatalist[[i]] = tmp-GCMhvardatalist[[i]]
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      GCMpvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      GCMpvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      GCMpvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      GCMpvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist[[i]],NA)
    } 
  }
  
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(GCMpvardatalist,mean,na.rm=TRUE)


### LOCA historical + Livneh

LOCAhvardatalist = list()
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
  if(type!="ann"){
    GCMhvardatalist[[i]] = tmp
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      LOCAhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      LOCAhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      LOCAhvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      LOCAhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist[[i]],NA)
    } 
  }
  
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCAhvardatalist,mean,na.rm=TRUE)

### LOCA projected change

LOCApvardatalist = list()
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
  if(type!="ann"){
    LOCApvardatalist[[i]] = tmp-LOCAhvardatalist[[i]]
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      LOCApvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      LOCApvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      LOCApvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      LOCApvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist[[i]],NA)
    } 
  }
  
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCApvardatalist,mean,na.rm=TRUE)

#######
# projected changes

GCMchange = LOCAchange = GCMproj = LOCAproj = GCMhist = LOCAhist = array(NA,dim=c(length(lon),ncol=length(lat),26))
OBS = LOCAhvardatalist[[27]]
for(i in 1:26){
  GCMchange[,,i] = GCMpvardatalist[[i]]-GCMhvardatalist[[i]]
  LOCAchange[,,i] = LOCApvardatalist[[i]]-LOCAhvardatalist[[i]]
  GCMproj[,,i] = GCMpvardatalist[[i]]
  LOCAproj[,,i] = LOCApvardatalist[[i]]
  GCMhist[,,i] = GCMhvardatalist[[i]]
  LOCAhist[,,i] = LOCAhvardatalist[[i]]
}

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

GCMunweightedmean_hist = apply(GCMhist,c(1,2),mean,na.rm=TRUE)
LOCAunweightedmean_hist = apply(LOCAhist,c(1,2),mean,na.rm=TRUE)

GCMskillmean_hist = GCMSIhmean_hist = GCMSIcmean_hist = GCMBMAmean_hist = GCMunweightedmean_hist
LOCAskillmean_hist = LOCASIhmean_hist = LOCASIcmean_hist = LOCABMAmean_hist = LOCAunweightedmean_hist

for(i in 1:26){
  
  ## skill mean
  tmpG = GCMhist[,,i]*GCMweights$Ws[i]
  tmpL = LOCAhist[,,i]*LOCAweights$Ws[i]
  if(i==1){
    GCMskillmean_hist = tmpG
    LOCAskillmean_hist = tmpL
  } else {
    GCMskillmean_hist = GCMskillmean_hist+tmpG
    LOCAskillmean_hist = LOCAskillmean_hist+tmpL
  }
  
  ## skill+ind hist only
  tmpG = GCMhist[,,i]*GCMweights$Wh[i]
  tmpL = LOCAhist[,,i]*LOCAweights$Wh[i]
  if(i==1){
    GCMSIhmean_hist = tmpG
    LOCASIhmean_hist = tmpL
  } else {
    GCMSIhmean_hist = GCMSIhmean_hist+tmpG
    LOCASIhmean_hist = LOCASIhmean_hist+tmpL
  }
  
  ## skill+ind hist and change
  tmpG = GCMhist[,,i]*GCMweights$Wc[i]
  tmpL = LOCAhist[,,i]*LOCAweights$Wc[i]
  if(i==1){
    GCMSIcmean_hist = tmpG
    LOCASIcmean_hist = tmpL
  } else {
    GCMSIcmean_hist = GCMSIcmean_hist+tmpG
    LOCASIcmean_hist = LOCASIcmean_hist+tmpL
  }
  
  ## BMA
  tmpG = GCMhist[,,i]*GCMweights$BMA[i]
  tmpL = LOCAhist[,,i]*LOCAweights$BMA[i]
  if(i==1){
    GCMBMAmean_hist = tmpG
    LOCABMAmean_hist = tmpL
  } else {
    GCMBMAmean_hist = GCMBMAmean_hist+tmpG
    LOCABMAmean_hist = LOCABMAmean_hist+tmpL
  }
  
}

LOCAunweightedmean_bias = LOCAunweightedmean_hist-OBS
LOCAskillmean_bias = LOCAskillmean_hist-OBS
LOCASIhmean_bias = LOCASIhmean_hist-OBS
LOCASIcmean_bias = LOCASIcmean_hist-OBS
LOCABMAmean_bias = LOCABMAmean_hist-OBS

GCMunweightedmean_bias = GCMunweightedmean_hist-OBS
GCMskillmean_bias = GCMskillmean_hist-OBS
GCMSIhmean_bias = GCMSIhmean_hist-OBS
GCMSIcmean_bias = GCMSIcmean_hist-OBS
GCMBMAmean_bias = GCMBMAmean_hist-OBS

######
# Calculate change means (weighted and unweighted)

GCMunweightedmean_change = apply(GCMchange,c(1,2),mean,na.rm=TRUE)
LOCAunweightedmean_change = apply(LOCAchange,c(1,2),mean,na.rm=TRUE)

GCMskillmean_change = GCMSIhmean_change = GCMSIcmean_change = GCMBMAmean_change = GCMunweightedmean_change
LOCAskillmean_change = LOCASIhmean_change = LOCASIcmean_change = LOCABMAmean_change = LOCAunweightedmean_change

for(i in 1:26){
  
  ## skill mean
  tmpG = GCMchange[,,i]*GCMweights$Ws[i]
  tmpL = LOCAchange[,,i]*LOCAweights$Ws[i]
  if(i==1){
    GCMskillmean_change = tmpG
    LOCAskillmean_change = tmpL
  } else {
    GCMskillmean_change = GCMskillmean_change+tmpG
    LOCAskillmean_change = LOCAskillmean_change+tmpL
  }
  
  ## skill+ind hist only
  tmpG = GCMchange[,,i]*GCMweights$Wh[i]
  tmpL = LOCAchange[,,i]*LOCAweights$Wh[i]
  if(i==1){
    GCMSIhmean_change = tmpG
    LOCASIhmean_change = tmpL
  } else {
    GCMSIhmean_change = GCMSIhmean_change+tmpG
    LOCASIhmean_change = LOCASIhmean_change+tmpL
  }
  
  ## skill+ind hist and change
  tmpG = GCMchange[,,i]*GCMweights$Wc[i]
  tmpL = LOCAchange[,,i]*LOCAweights$Wc[i]
  if(i==1){
    GCMSIcmean_change = tmpG
    LOCASIcmean_change = tmpL
  } else {
    GCMSIcmean_change = GCMSIcmean_change+tmpG
    LOCASIcmean_change = LOCASIcmean_change+tmpL
  }
  
  ## BMA
  tmpG = GCMchange[,,i]*GCMweights$BMA[i]
  tmpL = LOCAchange[,,i]*LOCAweights$BMA[i]
  if(i==1){
    GCMBMAmean_change = tmpG
    LOCABMAmean_change = tmpL
  } else {
    GCMBMAmean_change = GCMBMAmean_change+tmpG
    LOCABMAmean_change = LOCABMAmean_change+tmpL
  }
  
}

######
# Calculate change variance (weighted and unweighted)

GCMunweightedvar_change = GCMskillvar_change = GCMSIhvar_change = GCMSIcvar_change = GCMBMAvar_change = GCMunweightedmean_change
LOCAunweightedvar_change = LOCAskillvar_change = LOCASIhvar_change = LOCASIcvar_change = LOCABMAvar_change = LOCAunweightedmean_change

for(R in 1:length(lon)){
  for(C in 1:length(lat)){
    if(all(is.na(GCMchange[R,C,])==TRUE)==FALSE){
      GCMunweightedvar_change[R,C] = var(GCMchange[R,C,],na.rm=TRUE)
      LOCAunweightedvar_change[R,C] = var(LOCAchange[R,C,],na.rm=TRUE)
      
      GCMskillvar_change[R,C] = weighted.var(x=GCMchange[R,C,],w=GCMweights$Ws,na.rm=TRUE)
      LOCAskillvar_change[R,C] = weighted.var(x=LOCAchange[R,C,],w=LOCAweights$Ws,na.rm=TRUE)
      
      GCMSIhvar_change[R,C] = weighted.var(x=GCMchange[R,C,],w=GCMweights$Wh,na.rm=TRUE)
      LOCASIhvar_change[R,C] = weighted.var(x=LOCAchange[R,C,],w=LOCAweights$Wh,na.rm=TRUE)
      
      GCMSIcvar_change[R,C] = weighted.var(x=GCMchange[R,C,],w=GCMweights$Wc,na.rm=TRUE)
      LOCASIcvar_change[R,C] = weighted.var(x=LOCAchange[R,C,],w=LOCAweights$Wc,na.rm=TRUE)
      
      GCMBMAvar_change[R,C] = weighted.var(x=GCMchange[R,C,],w=GCMweights$BMA,na.rm=TRUE)
      LOCABMAvar_change[R,C] = weighted.var(x=LOCAchange[R,C,],w=LOCAweights$BMA,na.rm=TRUE)
      message("Finished calcs for R: ",R," and C: ",C)
    }
  }
}

#######

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
# gather means
meansdat = NULL
for(s in 1:5){
  if(s<=4){
    regionidx = which(regionmasklist[[s]]==1)
    histmeans_GCM = c(mean(GCMunweightedmean_hist[regionidx],na.rm=TRUE),mean(GCMskillmean_hist[regionidx],na.rm=TRUE),mean(GCMSIhmean_hist[regionidx],na.rm=TRUE),mean(GCMSIcmean_hist[regionidx],na.rm=TRUE),mean(GCMBMAmean_hist[regionidx],na.rm=TRUE)) 
    histmeans_LOCA = c(mean(LOCAunweightedmean_hist[regionidx],na.rm=TRUE),mean(LOCAskillmean_hist[regionidx],na.rm=TRUE),mean(LOCASIhmean_hist[regionidx],na.rm=TRUE),mean(LOCASIcmean_hist[regionidx],na.rm=TRUE),mean(LOCABMAmean_hist[regionidx],na.rm=TRUE)) 
    
    changemeans_GCM = c(mean(GCMunweightedmean_change[regionidx],na.rm=TRUE),mean(GCMskillmean_change[regionidx],na.rm=TRUE),mean(GCMSIhmean_change[regionidx],na.rm=TRUE),mean(GCMSIcmean_change[regionidx],na.rm=TRUE),mean(GCMBMAmean_change[regionidx],na.rm=TRUE)) 
    changemeans_LOCA = c(mean(LOCAunweightedmean_change[regionidx],na.rm=TRUE),mean(LOCAskillmean_change[regionidx],na.rm=TRUE),mean(LOCASIhmean_change[regionidx],na.rm=TRUE),mean(LOCASIcmean_change[regionidx],na.rm=TRUE),mean(LOCABMAmean_change[regionidx],na.rm=TRUE)) 
    
    changevars_GCM = c(mean(GCMunweightedvar_change[regionidx],na.rm=TRUE),mean(GCMskillvar_change[regionidx],na.rm=TRUE),mean(GCMSIhvar_change[regionidx],na.rm=TRUE),mean(GCMSIcvar_change[regionidx],na.rm=TRUE),mean(GCMBMAvar_change[regionidx],na.rm=TRUE)) 
    changevars_LOCA = c(mean(LOCAunweightedvar_change[regionidx],na.rm=TRUE),mean(LOCAskillvar_change[regionidx],na.rm=TRUE),mean(LOCASIhvar_change[regionidx],na.rm=TRUE),mean(LOCASIcvar_change[regionidx],na.rm=TRUE),mean(LOCABMAvar_change[regionidx],na.rm=TRUE)) 
    
    obs = rep(mean(OBS[regionidx],na.rm=TRUE),10)
    region = rep(statelist[s],10)
    
    RMSE_GCM = c(sqrt(mean((GCMunweightedmean_bias[regionidx])^2,na.rm=TRUE)),
                 sqrt(mean((GCMskillmean_bias[regionidx])^2,na.rm=TRUE)),
                 sqrt(mean((GCMSIhmean_bias[regionidx])^2,na.rm=TRUE)),
                 sqrt(mean((GCMSIcmean_bias[regionidx])^2,na.rm=TRUE)),
                 sqrt(mean((GCMBMAmean_bias[regionidx])^2,na.rm=TRUE)))
    
    RMSE_LOCA = c(sqrt(mean((LOCAunweightedmean_bias[regionidx])^2,na.rm=TRUE)),
                  sqrt(mean((LOCAskillmean_bias[regionidx])^2,na.rm=TRUE)),
                  sqrt(mean((LOCASIhmean_bias[regionidx])^2,na.rm=TRUE)),
                  sqrt(mean((LOCASIcmean_bias[regionidx])^2,na.rm=TRUE)),
                  sqrt(mean((LOCABMAmean_bias[regionidx])^2,na.rm=TRUE)))
    
  } else {
    histmeans_GCM = c(mean(GCMunweightedmean_hist,na.rm=TRUE),mean(GCMskillmean_hist,na.rm=TRUE),mean(GCMSIhmean_hist,na.rm=TRUE),mean(GCMSIcmean_hist,na.rm=TRUE),mean(GCMBMAmean_hist,na.rm=TRUE)) 
    histmeans_LOCA = c(mean(LOCAunweightedmean_hist,na.rm=TRUE),mean(LOCAskillmean_hist,na.rm=TRUE),mean(LOCASIhmean_hist,na.rm=TRUE),mean(LOCASIcmean_hist,na.rm=TRUE),mean(LOCABMAmean_hist,na.rm=TRUE)) 
    changemeans_GCM = c(mean(GCMunweightedmean_change,na.rm=TRUE),mean(GCMskillmean_change,na.rm=TRUE),mean(GCMSIhmean_change,na.rm=TRUE),mean(GCMSIcmean_change,na.rm=TRUE),mean(GCMBMAmean_change,na.rm=TRUE)) 
    changemeans_LOCA = c(mean(LOCAunweightedmean_change,na.rm=TRUE),mean(LOCAskillmean_change,na.rm=TRUE),mean(LOCASIhmean_change,na.rm=TRUE),mean(LOCASIcmean_change,na.rm=TRUE),mean(LOCABMAmean_change,na.rm=TRUE)) 
    changevars_GCM = c(mean(GCMunweightedvar_change,na.rm=TRUE),mean(GCMskillvar_change,na.rm=TRUE),mean(GCMSIhvar_change,na.rm=TRUE),mean(GCMSIcvar_change,na.rm=TRUE),mean(GCMBMAvar_change,na.rm=TRUE)) 
    changevars_LOCA = c(mean(LOCAunweightedvar_change,na.rm=TRUE),mean(LOCAskillvar_change,na.rm=TRUE),mean(LOCASIhvar_change,na.rm=TRUE),mean(LOCASIcvar_change,na.rm=TRUE),mean(LOCABMAvar_change,na.rm=TRUE)) 
    
    obs = rep(mean(OBS,na.rm=TRUE),10)
    region = rep("full",10)
    
    RMSE_GCM = c(sqrt(mean((GCMunweightedmean_bias)^2,na.rm=TRUE)),
                 sqrt(mean((GCMskillmean_bias)^2,na.rm=TRUE)),
                 sqrt(mean((GCMSIhmean_bias)^2,na.rm=TRUE)),
                 sqrt(mean((GCMSIcmean_bias)^2,na.rm=TRUE)),
                 sqrt(mean((GCMBMAmean_bias)^2,na.rm=TRUE)))
    
    RMSE_LOCA = c(sqrt(mean((LOCAunweightedmean_bias)^2,na.rm=TRUE)),
                  sqrt(mean((LOCAskillmean_bias)^2,na.rm=TRUE)),
                  sqrt(mean((LOCASIhmean_bias)^2,na.rm=TRUE)),
                  sqrt(mean((LOCASIcmean_bias)^2,na.rm=TRUE)),
                  sqrt(mean((LOCABMAmean_bias)^2,na.rm=TRUE)))
    
  }
  
  histmeans = c(histmeans_GCM,histmeans_LOCA)
  changemeans = c(changemeans_GCM,changemeans_LOCA)
  changevars = c(changevars_GCM,changevars_LOCA)
  rmse = c(RMSE_GCM,RMSE_LOCA)
  group = rep(c("unweighted","skill","SI-h","SI-c","BMA"),2)
  DS = rep(c("CMIP5","LOCA"),each=5)
  
  meansframe = data.frame(group,region,DS,histmeans,obs,changemeans,changevars,rmse)
  meansdat = rbind(meansdat,meansframe)
}

meansdat$bias = meansdat$histmeans-meansdat$obs
#save(list="meansdat",file="WeightedMeansVars_prann_pronly_modi.Rdata")
meansdat$region <- factor(meansdat$region,levels = c('full','new mexico','texas','oklahoma','louisiana'),ordered = TRUE)

ggplot(meansdat, aes(x=region, y=histmeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) + geom_point(data=meansdat, mapping=aes(x=region,y=obs),size=5,shape=23,fill="blue") +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Historical Ens. Means against Livneh")+xlab("Region")+ylab("Historical Mean (mm)")
ggplot(meansdat, aes(x=region, y=bias)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Bias of Ens. Means against Livneh")+xlab("Region")+ylab("Bias (mm)")
ggplot(meansdat, aes(x=region, y=rmse)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("RMSE of Ens. Means against Livneh")+xlab("Region")+ylab("RMSE (mm)")
ggplot(meansdat, aes(x=region, y=changemeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Mean Projected Changes - RCP8.5 2070-2099")+xlab("Region")+ylab("Change (mm)")
ggplot(meansdat, aes(x=region, y=changevars)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Variance Projected Changes - RCP8.5 2070-2099")+xlab("Region")+ylab("Var (mm^2)")

####
# BMA all grab

load("/home/woot0002/DS_ind/BMAposterior_meansandvars.Rdata")

meansdatBMA = NULL
for(i in 1:100){
  GCMBMAhm = GCMBMAhistmean[,,i]
  GCMBMAcm = GCMBMAchangemean[,,i]
  GCMBMAcv = GCMBMAchangevar[,,i]
  
  LOCABMAhm = LOCABMAhistmean[,,i]
  LOCABMAcm = LOCABMAchangemean[,,i]
  LOCABMAcv = LOCABMAchangevar[,,i]
  
for(s in 1:5){
  if(s<=4){
    regionidx = which(regionmasklist[[s]]==1)
    histmeans_GCM = mean(GCMBMAhm[regionidx],na.rm=TRUE)
    histmeans_LOCA = mean(LOCABMAhm[regionidx],na.rm=TRUE)
    
    changemeans_GCM = mean(GCMBMAcm[regionidx],na.rm=TRUE)
    changemeans_LOCA = mean(LOCABMAcm[regionidx],na.rm=TRUE)
    
    changevars_GCM = mean(GCMBMAcv[regionidx],na.rm=TRUE)
    changevars_LOCA = mean(LOCABMAcv[regionidx],na.rm=TRUE)
    
    obs = mean(OBS[regionidx],na.rm=TRUE)
    region = statelist[s]
    
    RMSE_GCM = sqrt(mean((GCMBMAhm[regionidx]-OBS[regionidx])^2,na.rm=TRUE))
    RMSE_LOCA = sqrt(mean((LOCABMAhm[regionidx]-OBS[regionidx])^2,na.rm=TRUE))
    
  } else {
    histmeans_GCM = mean(GCMBMAhm,na.rm=TRUE)
    histmeans_LOCA = mean(LOCABMAhm,na.rm=TRUE)
    
    changemeans_GCM = mean(GCMBMAcm,na.rm=TRUE)
    changemeans_LOCA = mean(LOCABMAcm,na.rm=TRUE)
    
    changevars_GCM = mean(GCMBMAcv,na.rm=TRUE)
    changevars_LOCA = mean(LOCABMAcv,na.rm=TRUE)
    
    obs = mean(OBS,na.rm=TRUE)
    region="full"
    RMSE_GCM = sqrt(mean((GCMBMAhm-OBS)^2,na.rm=TRUE))
    RMSE_LOCA = sqrt(mean((LOCABMAhm-OBS)^2,na.rm=TRUE))
  }
  
  histmeans = c(histmeans_GCM,histmeans_LOCA)
  changemeans = c(changemeans_GCM,changemeans_LOCA)
  changevars = c(changevars_GCM,changevars_LOCA)
  rmse = c(RMSE_GCM,RMSE_LOCA)
  DS = c("CMIP5","LOCA")
  BMAmem = i
  
  meansframe = data.frame(BMAmem,region,DS,histmeans,obs,changemeans,changevars,rmse)
  meansdatBMA = rbind(meansdatBMA,meansframe)
}
}

meansdatBMA$changesd = sqrt(meansdatBMA$changevars)
meansdatBMA$bias = meansdatBMA$histmeans-meansdatBMA$obs

meansdatBMA$region <- factor(meansdatBMA$region,levels = c('full','new mexico','texas','oklahoma','louisiana'),ordered = TRUE)

meansdat$changesd = sqrt(meansdat$changevars)


ggplot(meansdatBMA, aes(x=region, y=histmeans)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat, mapping=aes(x=region,y=histmeans,colour = factor(group),shape=factor(DS)),size=2.5) + geom_point(data=meansdat, mapping=aes(x=region,y=obs),size=2.5,shape=23,fill="blue") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean")+xlab("")+ylab("Mean (mm)")+ylim(0,1550) +facet_grid(cols=vars(DS))

ggplot(meansdatBMA, aes(x=region, y=bias)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat, mapping=aes(x=region,y=bias,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean Bias")+xlab("")+ylab("Bias (mm)")+ylim(range(meansdatBMA$bias,na.rm=TRUE)) +facet_grid(cols=vars(DS))+geom_hline(yintercept=0,linetype="dashed")

ggplot(meansdatBMA, aes(x=region, y=rmse)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat, mapping=aes(x=region,y=rmse,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean RMSE")+xlab("")+ylab("RMSE (mm)")+ylim(range(meansdatBMA$rmse,na.rm=TRUE)) +facet_grid(cols=vars(DS))+geom_hline(yintercept=0,linetype="dashed")

ggplot(meansdatBMA, aes(x=region, y=changemeans)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat, mapping=aes(x=region,y=changemeans,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Ens. Mean Projected Change")+xlab("")+ylab("Change (mm)")+ylim(range(meansdatBMA$changemeans,na.rm=TRUE)) +facet_grid(cols=vars(DS))+geom_hline(yintercept=0,linetype="dashed")

ggplot(meansdatBMA, aes(x=region, y=changesd)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat, mapping=aes(x=region,y=changesd,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Ens.Std. Deviation Projected Change")+xlab("")+ylab("Std. Deviation (mm)")+ylim(range(meansdatBMA$changesd,na.rm=TRUE)) +facet_grid(cols=vars(DS))+geom_hline(yintercept=0,linetype="dashed")

######

######
# gather means
memsdat = NULL

for(m in 1:26){
  
  GCMhtmp = GCMhvardatalist[[m]]
  LOCAhtmp = LOCAhvardatalist[[m]]
  GCMptmp = GCMpvardatalist[[m]]
  LOCAptmp = LOCApvardatalist[[m]]
  
  for(s in 1:5){
    if(s<=4){
      regionidx = which(regionmasklist[[s]]==1)
      histmeans_GCM = mean(GCMhtmp[regionidx],na.rm=TRUE)
      histmeans_LOCA = mean(LOCAhtmp[regionidx],na.rm=TRUE)
      
      changemeans_GCM = mean(GCMptmp[regionidx]-GCMhtmp[regionidx],na.rm=TRUE)
      changemeans_LOCA = mean(LOCAptmp[regionidx]-LOCAhtmp[regionidx],na.rm=TRUE)
      
      obs = rep(mean(OBS[regionidx],na.rm=TRUE),2)
      region = rep(statelist[s],2)
      
      GCMbias = GCMhtmp[regionidx]-OBS[regionidx]
      LOCAbias = LOCAhtmp[regionidx]-OBS[regionidx]
      
      RMSE_GCM = sqrt(mean((GCMbias)^2,na.rm=TRUE))
      RMSE_LOCA = sqrt(mean((LOCAbias)^2,na.rm=TRUE))
      
    } else {
      histmeans_GCM = mean(GCMhtmp,na.rm=TRUE)
      histmeans_LOCA = mean(LOCAhtmp,na.rm=TRUE)
      
      changemeans_GCM = mean(GCMptmp-GCMhtmp,na.rm=TRUE)
      changemeans_LOCA = mean(LOCAptmp-LOCAhtmp,na.rm=TRUE)
      
      obs = rep(mean(OBS,na.rm=TRUE),2)
      region = rep("full",2)
      
      GCMbias = GCMhtmp-OBS
      LOCAbias = LOCAhtmp-OBS
      
      RMSE_GCM = sqrt(mean((GCMbias)^2,na.rm=TRUE))
      RMSE_LOCA = sqrt(mean((LOCAbias)^2,na.rm=TRUE))
    }
    
    histmeans = c(histmeans_GCM,histmeans_LOCA)
    changemeans = c(changemeans_GCM,changemeans_LOCA)
    bias = histmeans-obs
    rmse = c(RMSE_GCM,RMSE_LOCA)
    GCM = rep(GCMweights$GCM[m],2)
    DS = c("CMIP5","LOCA")
    memsframe = data.frame(GCM,region,DS,histmeans,obs,changemeans,bias,rmse)
    memsdat = rbind(memsdat,memsframe)
  }
message("GCM calcs finished for GCM ",m," / 26")
}

memsdat$region <- factor(memsdat$region,levels = c('full','new mexico','texas','oklahoma','louisiana'),ordered = TRUE)


ggplot(memsdat, aes(x=region, y=histmeans, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=memsdat, mapping=aes(x=region,y=obs),size=2.5,shape=23,fill="blue") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("Hist Ens. Mean")+xlab("")+ylab("Mean (mm)")+ylim(range(memsdat$histmeans))

ggplot(memsdat, aes(x=region, y=bias, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("Bias of Hist Ens. Mean")+xlab("")+ylab("Bias (mm)")+ylim(range(memsdat$bias))

ggplot(memsdat, aes(x=region, y=rmse, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("RMSE of Hist Ens. Mean")+xlab("")+ylab("RMSE (mm)")+ylim(range(memsdat$rmse))

ggplot(memsdat, aes(x=region, y=changemeans, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("Ens. Mean Projected Change")+xlab("")+ylab("Change (mm)")+ylim(range(memsdat$changemeans))

#######
# Get SD of change means for BMA sims and raw ensembles

aggregate(memsdat$changemeans,by=list(regionDS = paste(memsdat$region,memsdat$DS,sep="_")),sd,na.rm=TRUE)
aggregate(meansdatBMA$changemeans,by=list(regionDS = paste(meansdatBMA$region,meansdatBMA$DS,sep="_")),sd,na.rm=TRUE)

aggregate(memsdat$histmeans,by=list(regionDS = paste(memsdat$region,memsdat$DS,sep="_")),sd,na.rm=TRUE)
aggregate(meansdatBMA$histmeans,by=list(regionDS = paste(meansdatBMA$region,meansdatBMA$DS,sep="_")),sd,na.rm=TRUE)
