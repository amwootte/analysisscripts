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

#weighted.var <- function(x, w, na.rm = FALSE) {
#  if (na.rm) {
#    w <- w[i <- !is.na(x)]
#    x <- x[i]
#  }
#  sum.w <- sum(w)
#  sum.w2 <- sum(w^2)
#  mean.w <- sum(x * w) / sum(w)
#  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
#                                       na.rm)
#}

setwd("/home/woot0002/DS_ind/")

load(file="Sanderson_EnsembleWeights_v4pronly_ann.Rdata")
BMAweights = read.table("posterior_BMA_combo.txt")
BMAweightst = t(BMAweights)

GCMhdat = cbind(GCMhdat,BMAweightst) 

BMAweights = read.table("posterior_BMA_combo_LOCA.txt")
BMAweightst = t(BMAweights)
LOCAhdat = cbind(LOCAhdat,BMAweightst)

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
# Calculate historical means (weighted and unweighted)

GCMBMAhistmean = LOCABMAhistmean = array(NA,dim=c(length(lon),length(lat),100))

for(b in 10:109){
for(i in 1:26){
  ## BMA
  tmpG = GCMhist[,,i]*GCMweights[i,b]
  tmpL = LOCAhist[,,i]*LOCAweights[i,b]
  if(i==1){
    GCMBMAmean_hist = tmpG
    LOCABMAmean_hist = tmpL
  } else {
    GCMBMAmean_hist = GCMBMAmean_hist+tmpG
    LOCABMAmean_hist = LOCABMAmean_hist+tmpL
  }
}
  GCMBMAhistmean[,,(b-9)]=GCMBMAmean_hist
  LOCABMAhistmean[,,(b-9)]=LOCABMAmean_hist
}

######
# Calculate change means (weighted and unweighted)

GCMBMAchangemean = LOCABMAchangemean = array(NA,dim=c(length(lon),length(lat),100))

for(b in 10:109){
  for(i in 1:26){
    ## BMA
    tmpG = GCMchange[,,i]*GCMweights[i,b]
    tmpL = LOCAchange[,,i]*LOCAweights[i,b]
    if(i==1){
      GCMBMAmean_change = tmpG
      LOCABMAmean_change = tmpL
    } else {
      GCMBMAmean_change = GCMBMAmean_change+tmpG
      LOCABMAmean_change = LOCABMAmean_change+tmpL
    }
  }
  GCMBMAchangemean[,,(b-9)]=GCMBMAmean_change
  LOCABMAchangemean[,,(b-9)]=LOCABMAmean_change
}

######
# ensemble variance 

GCMBMAchangevar = LOCABMAchangevar = array(NA,dim=c(length(lon),length(lat),100))

for(b in 10:109){
  
  GCMBMAvar_change = LOCABMAvar_change = matrix(NA,nrow=length(lon),ncol=length(lat))
  
  for(R in 1:length(lon)){
    for(C in 1:length(lat)){
      if(all(is.na(GCMchange[R,C,])==TRUE)==FALSE){
        GCMBMAvar_change[R,C] = weighted.var(GCMchange[R,C,],w=GCMweights[,b],na.rm=TRUE)
        LOCABMAvar_change[R,C] = weighted.var(LOCAchange[R,C,],w=LOCAweights[,b],na.rm=TRUE)
        message("Finished calcs for R: ",R," and C: ",C)
      }
    }
  }
  
  GCMBMAchangevar[,,(b-9)]=GCMBMAvar_change
  LOCABMAchangevar[,,(b-9)]=LOCABMAvar_change
  message("Finished calcs for BMA member: ",b-9," / ",100)
}

save(list=c("GCMBMAhistmean","GCMBMAchangemean","GCMBMAchangevar","LOCABMAhistmean","LOCABMAchangemean","LOCABMAchangevar"),file="/home/woot0002/DS_ind/BMAposterior_meansandvars.Rdata")

