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

var = varin = "pr"
type="ann"
weightingused = "full"
stateapplied = "full"

if(weightingused=="full"){
  load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,".Rdata",sep=""))
  BMAweightsGCM = read.table("posterior_BMA_combo.txt")
  BMAweightsLOCA = read.table("posterior_BMA_combo_LOCA.txt")
} else {
  load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,"_",weightingused,".Rdata",sep=""))
  BMAweightsGCM = read.table(paste("posterior_BMA_combo_",var,"_",weightingused,".txt",sep=""))
  BMAweightsLOCA = read.table(paste("posterior_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
}

BMAweightst = t(BMAweightsGCM)
GCMhdat = cbind(GCMhdat,BMAweightst) 

BMAweightst = t(BMAweightsLOCA)
LOCAhdat = cbind(LOCAhdat,BMAweightst)

GCMweights= GCMhdat
LOCAweights = LOCAhdat


# precip files
GCMfiles_pr = system("ls /home/woot0002/GCMs/regrid/pr_*histclimo*.nc",intern=TRUE)
LOCAfiles_pr = system("ls /home/woot0002/LOCA/regrid/pr_*histclimo*.nc",intern=TRUE)
GCMprojfiles_pr = system("ls /home/woot0002/GCMs/regrid/pr_*projclimo*.nc",intern=TRUE)
LOCAprojfiles_pr = system("ls /home/woot0002/LOCA/regrid/pr_*projclimo*.nc",intern=TRUE)
LIVNEHfile_pr = system("ls /home/woot0002/monthlyclimo/pr_day*livneh*.nc",intern=TRUE)

# tasmax files
GCMfiles_tmax = system("ls /home/woot0002/GCMs/regrid/tasmax_*histclimo*.nc",intern=TRUE)
LOCAfiles_tmax = system("ls /home/woot0002/LOCA/regrid/tasmax_*histclimo*.nc",intern=TRUE)
GCMprojfiles_tmax = system("ls /home/woot0002/GCMs/regrid/tasmax_*projclimo*.nc",intern=TRUE)
LOCAprojfiles_tmax = system("ls /home/woot0002/LOCA/regrid/tasmax_*projclimo*.nc",intern=TRUE)
LIVNEHfile_tmax = system("ls /home/woot0002/monthlyclimo/tasmax_day*livneh*.nc",intern=TRUE)

# subset files down
load("/home/woot0002/DS_ind/manuscript1/GCMlist.Rdata")

GCM_hfiles_pr = GCM_pfiles_pr = LOCA_hfiles_pr = LOCA_pfiles_pr = c()
GCM_hfiles_tmax = GCM_pfiles_tmax = LOCA_hfiles_tmax = LOCA_pfiles_tmax = c()

for(i in 1:length(GCMlist)){
  #pr 
  GCM_hfiles_pr[i] = GCMfiles_pr[grep(paste(GCMlist[i],"_",sep=""),GCMfiles_pr)]
  GCM_pfiles_pr[i] = GCMprojfiles_pr[grep(paste(GCMlist[i],"_",sep=""),GCMprojfiles_pr)]
  LOCA_hfiles_pr[i] = LOCAfiles_pr[grep(paste(GCMlist[i],"_",sep=""),LOCAfiles_pr)]
  LOCA_pfiles_pr[i] = LOCAprojfiles_pr[grep(paste(GCMlist[i],"_",sep=""),LOCAprojfiles_pr)]
  #tmax
  GCM_hfiles_tmax[i] = GCMfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),GCMfiles_tmax)]
  GCM_pfiles_tmax[i] = GCMprojfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),GCMprojfiles_tmax)]
  LOCA_hfiles_tmax[i] = LOCAfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),LOCAfiles_tmax)]
  LOCA_pfiles_tmax[i] = LOCAprojfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),LOCAprojfiles_tmax)]
}

###
# create full filelist + metadata table - historical

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_hfiles_pr,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMhdat = filelist2[,c(2,3,4,6)]
names(GCMhdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_hfiles_pr,"/",fixed=TRUE))
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
GCMgroup_pr = c(GCM_hfiles_pr,LIVNEHfile_pr)
LOCAgroup_pr = c(LOCA_hfiles_pr,LIVNEHfile_pr)

GCMgroup_tmax = c(GCM_hfiles_tmax,LIVNEHfile_tmax)
LOCAgroup_tmax = c(LOCA_hfiles_tmax,LIVNEHfile_tmax)


###
# create full filelist + metadata table - projected

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_pfiles_pr,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMpdat = filelist2[,c(2,3,4,6)]
names(GCMpdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_pfiles_pr,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCApdat = filelist2[,c(2,3,4,6)]
names(LOCApdat) = names(GCMpdat)

# all files
GCMpgroup_pr = GCM_pfiles_pr
LOCApgroup_pr = LOCA_pfiles_pr
GCMpgroup_tmax = GCM_pfiles_tmax
LOCApgroup_tmax = LOCA_pfiles_tmax

######
# Gather data

ncvarname = "prclimo"
### GCM hist + Livneh - pr
GCMhvardatalist_pr = list()
for(i in 1:length(GCMgroup_pr)){
  nctest = nc_open(GCMgroup_pr[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  GCMhvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  if(stateapplied=="full"){
    GCMhvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist_pr[[i]],NA)
  } else {
    GCMhvardatalist_pr[[i]] = ifelse(regionmask==1,GCMhvardatalist_pr[[i]],NA)
  }
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon");
  nc_close(nctest)
}

sapply(GCMhvardatalist_pr,mean,na.rm=TRUE)

### GCM projected change - pr
GCMpvardatalist_pr = list()
for(i in 1:length(GCMpgroup_pr)){
  nctest = nc_open(GCMpgroup_pr[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  GCMpvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  if(stateapplied=="full"){
    GCMpvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist_pr[[i]],NA)
  } else {
    GCMpvardatalist_pr[[i]] = ifelse(regionmask==1,GCMpvardatalist_pr[[i]],NA)
  }
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(GCMpvardatalist_pr,mean,na.rm=TRUE)


### LOCA historical + Livneh - pr

LOCAhvardatalist_pr = list()
for(i in 1:length(LOCAgroup_pr)){
  nctest = nc_open(LOCAgroup_pr[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  LOCAhvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  if(stateapplied=="full"){
    LOCAhvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist_pr[[i]],NA)
  } else{
    LOCAhvardatalist_pr[[i]] = ifelse(regionmask==1,LOCAhvardatalist_pr[[i]],NA)
  }
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCAhvardatalist_pr,mean,na.rm=TRUE)

### LOCA projected change - pr

LOCApvardatalist_pr = list()
for(i in 1:length(LOCApgroup_pr)){
  nctest = nc_open(LOCApgroup_pr[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  LOCApvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  if(stateapplied=="full"){
    LOCApvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist_pr[[i]],NA)
  } else{
    LOCApvardatalist_pr[[i]] = ifelse(regionmask==1,LOCApvardatalist_pr[[i]],NA)
  }
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCApvardatalist_pr,mean,na.rm=TRUE)


######
# Gather Data 2

ncvarname = "tmaxclimo"
### GCM hist + Livneh - tmax
GCMhvardatalist_tmax = list()
for(i in 1:length(GCMgroup_tmax)){
  nctest = nc_open(GCMgroup_tmax[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  GCMhvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
  if(stateapplied=="full"){
    GCMhvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist_tmax[[i]],NA)
  } else{
    GCMhvardatalist_tmax[[i]] = ifelse(regionmask==1,GCMhvardatalist_tmax[[i]],NA)
  }
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon");
  nc_close(nctest)
}

sapply(GCMhvardatalist_tmax,mean,na.rm=TRUE)

### GCM projected change - tmax
GCMpvardatalist_tmax = list()
for(i in 1:length(GCMpgroup_tmax)){
  nctest = nc_open(GCMpgroup_tmax[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  GCMpvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
  if(stateapplied=="full"){
    GCMpvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist_tmax[[i]],NA)
  } else{
    GCMpvardatalist_tmax[[i]] = ifelse(regionmask==1,GCMpvardatalist_tmax[[i]],NA)
  }
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(GCMpvardatalist_tmax,mean,na.rm=TRUE)

### LOCA historical + Livneh - tmax

LOCAhvardatalist_tmax = list()
for(i in 1:length(LOCAgroup_tmax)){
  nctest = nc_open(LOCAgroup_tmax[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  LOCAhvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
  if(stateapplied=="full"){
    LOCAhvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist_tmax[[i]],NA)
  } else{
    LOCAhvardatalist_tmax[[i]] = ifelse(regionmask==1,LOCAhvardatalist_tmax[[i]],NA)
  }
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCAhvardatalist_tmax,mean,na.rm=TRUE)

### LOCA projected change - tmax

LOCApvardatalist_tmax = list()
for(i in 1:length(LOCApgroup_tmax)){
  nctest = nc_open(LOCApgroup_tmax[i])
  idx = which(names(nctest$var)==ncvarname)
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  LOCApvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
  if(stateapplied=="full"){
    LOCApvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist_tmax[[i]],NA)
  } else{
    LOCApvardatalist_tmax[[i]] = ifelse(regionmask==1,LOCApvardatalist_tmax[[i]],NA)
  }
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCApvardatalist_tmax,mean,na.rm=TRUE)


#######
# projected changes - _pr

GCMchange_pr = LOCAchange_pr = GCMproj_pr = LOCAproj_pr = GCMhist_pr = LOCAhist_pr = array(NA,dim=c(length(lon),ncol=length(lat),26))
OBS_pr = LOCAhvardatalist_pr[[27]]
for(i in 1:26){
  GCMchange_pr[,,i] = GCMpvardatalist_pr[[i]]-GCMhvardatalist_pr[[i]]
  LOCAchange_pr[,,i] = LOCApvardatalist_pr[[i]]-LOCAhvardatalist_pr[[i]]
  GCMproj_pr[,,i] = GCMpvardatalist_pr[[i]]
  LOCAproj_pr[,,i] = LOCApvardatalist_pr[[i]]
  GCMhist_pr[,,i] = GCMhvardatalist_pr[[i]]
  LOCAhist_pr[,,i] = LOCAhvardatalist_pr[[i]]
}

#######
# projected changes - _tmax

GCMchange_tmax = LOCAchange_tmax = GCMproj_tmax = LOCAproj_tmax = GCMhist_tmax = LOCAhist_tmax = array(NA,dim=c(length(lon),ncol=length(lat),26))
OBS_tmax = LOCAhvardatalist_tmax[[27]]
for(i in 1:26){
  GCMchange_tmax[,,i] = GCMpvardatalist_tmax[[i]]-GCMhvardatalist_tmax[[i]]
  LOCAchange_tmax[,,i] = LOCApvardatalist_tmax[[i]]-LOCAhvardatalist_tmax[[i]]
  GCMproj_tmax[,,i] = GCMpvardatalist_tmax[[i]]
  LOCAproj_tmax[,,i] = LOCApvardatalist_tmax[[i]]
  GCMhist_tmax[,,i] = GCMhvardatalist_tmax[[i]]
  LOCAhist_tmax[,,i] = LOCAhvardatalist_tmax[[i]]
}

######
# Calculate historical means (weighted and unweighted)

GCMBMAhistmean_pr = LOCABMAhistmean_pr = array(NA,dim=c(length(lon),length(lat),100))
GCMBMAhistmean_tmax = LOCABMAhistmean_tmax = array(NA,dim=c(length(lon),length(lat),100))

for(b in 10:109){
for(i in 1:26){
  ## BMA
  tmpG = GCMhist_pr[,,i]*GCMweights[i,b]
  tmpL = LOCAhist_pr[,,i]*LOCAweights[i,b]
  if(i==1){
    GCMBMAmean_hist = tmpG
    LOCABMAmean_hist = tmpL
  } else {
    GCMBMAmean_hist = GCMBMAmean_hist+tmpG
    LOCABMAmean_hist = LOCABMAmean_hist+tmpL
  }
}
  GCMBMAhistmean_pr[,,(b-9)]=GCMBMAmean_hist
  LOCABMAhistmean_pr[,,(b-9)]=LOCABMAmean_hist
}

for(b in 10:109){
  for(i in 1:26){
    ## BMA
    tmpG = GCMhist_tmax[,,i]*GCMweights[i,b]
    tmpL = LOCAhist_tmax[,,i]*LOCAweights[i,b]
    if(i==1){
      GCMBMAmean_hist = tmpG
      LOCABMAmean_hist = tmpL
    } else {
      GCMBMAmean_hist = GCMBMAmean_hist+tmpG
      LOCABMAmean_hist = LOCABMAmean_hist+tmpL
    }
  }
  GCMBMAhistmean_tmax[,,(b-9)]=GCMBMAmean_hist
  LOCABMAhistmean_tmax[,,(b-9)]=LOCABMAmean_hist
}



######
# Calculate change means (weighted and unweighted)

GCMBMAchangemean_pr = LOCABMAchangemean_pr = array(NA,dim=c(length(lon),length(lat),100))
GCMBMAchangemean_tmax = LOCABMAchangemean_tmax = array(NA,dim=c(length(lon),length(lat),100))

for(b in 10:109){
  for(i in 1:26){
    ## BMA
    tmpG = GCMchange_pr[,,i]*GCMweights[i,b]
    tmpL = LOCAchange_pr[,,i]*LOCAweights[i,b]
    if(i==1){
      GCMBMAmean_change = tmpG
      LOCABMAmean_change = tmpL
    } else {
      GCMBMAmean_change = GCMBMAmean_change+tmpG
      LOCABMAmean_change = LOCABMAmean_change+tmpL
    }
  }
  GCMBMAchangemean_pr[,,(b-9)]=GCMBMAmean_change
  LOCABMAchangemean_pr[,,(b-9)]=LOCABMAmean_change
}


for(b in 10:109){
  for(i in 1:26){
    ## BMA
    tmpG = GCMchange_tmax[,,i]*GCMweights[i,b]
    tmpL = LOCAchange_tmax[,,i]*LOCAweights[i,b]
    if(i==1){
      GCMBMAmean_change = tmpG
      LOCABMAmean_change = tmpL
    } else {
      GCMBMAmean_change = GCMBMAmean_change+tmpG
      LOCABMAmean_change = LOCABMAmean_change+tmpL
    }
  }
  GCMBMAchangemean_tmax[,,(b-9)]=GCMBMAmean_change
  LOCABMAchangemean_tmax[,,(b-9)]=LOCABMAmean_change
}



######
# ensemble variance 

GCMBMAchangevar_pr = LOCABMAchangevar_pr = array(NA,dim=c(length(lon),length(lat),100))
GCMBMAchangevar_tmax = LOCABMAchangevar_tmax = array(NA,dim=c(length(lon),length(lat),100))

for(b in 10:109){
  
  GCMBMAvar_change = LOCABMAvar_change = matrix(NA,nrow=length(lon),ncol=length(lat))
  
  for(R in 1:length(lon)){
    for(C in 1:length(lat)){
      if(all(is.na(GCMchange_pr[R,C,])==TRUE)==FALSE){
        GCMBMAvar_change[R,C] = weighted.var(GCMchange_pr[R,C,],w=GCMweights[,b],na.rm=TRUE)
        LOCABMAvar_change[R,C] = weighted.var(LOCAchange_pr[R,C,],w=LOCAweights[,b],na.rm=TRUE)
        message("Finished calcs for R: ",R," and C: ",C)
      }
    }
  }
  
  GCMBMAchangevar_pr[,,(b-9)]=GCMBMAvar_change
  LOCABMAchangevar_pr[,,(b-9)]=LOCABMAvar_change
  message("Finished calcs for BMA member: ",b-9," / ",100)
}


for(b in 10:109){
  
  GCMBMAvar_change = LOCABMAvar_change = matrix(NA,nrow=length(lon),ncol=length(lat))
  
  for(R in 1:length(lon)){
    for(C in 1:length(lat)){
      if(all(is.na(GCMchange_tmax[R,C,])==TRUE)==FALSE){
        GCMBMAvar_change[R,C] = weighted.var(GCMchange_tmax[R,C,],w=GCMweights[,b],na.rm=TRUE)
        LOCABMAvar_change[R,C] = weighted.var(LOCAchange_tmax[R,C,],w=LOCAweights[,b],na.rm=TRUE)
        message("Finished calcs for R: ",R," and C: ",C)
      }
    }
  }
  
  GCMBMAchangevar_tmax[,,(b-9)]=GCMBMAvar_change
  LOCABMAchangevar_tmax[,,(b-9)]=LOCABMAvar_change
  message("Finished calcs for BMA member: ",b-9," / ",100)
}


save(list=c("GCMBMAhistmean_pr","GCMBMAchangemean_pr","GCMBMAchangevar_pr","LOCABMAhistmean_pr","LOCABMAchangemean_pr","LOCABMAchangevar_pr",
            "GCMBMAhistmean_tmax","GCMBMAchangemean_tmax","GCMBMAchangevar_tmax","LOCABMAhistmean_tmax","LOCABMAchangemean_tmax","LOCABMAchangevar_tmax"),file=paste("/home/woot0002/DS_ind/BMAposterior_meansandvars_",var,"_WU",weightingused,".Rdata",sep=""))

