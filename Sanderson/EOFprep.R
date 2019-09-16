###########################
#
# Sanderson - EOF preparation
#
###########################
#
# 1) Area Weighting
# 2) Normalization
# 3) Concatenation
# 4) Save to Rdata files
#
###########################

##########################
# Get files and information

library(ncdf4)
library(sp)
library(fields)
#library(smacof)

obsfile = "/home/woot0002/obs/regrid/METDATA_histclimo.nc"
LOCAfiles = system("ls /home/woot0002/LOCA/regrid/*.nc",intern=TRUE)
GCMfiles = system("ls /home/woot0002/GCMs/regrid/*.nc",intern=TRUE)

GCMdata = read.table("prCMIP5grab.csv",sep=",",header=TRUE)
GCMs = unique(GCMdata$model)

###
# obsdata
nctest = nc_open(obsfile)
lonobs = ncvar_get(nctest,"lon")
latobs = ncvar_get(nctest,"lat")
obstmax = ncvar_get(nctest,"tmaxclimo")
obstmax95 = ncvar_get(nctest,"tmax95climo")
obstmin = ncvar_get(nctest,"tminclimo")
obstmin32 = ncvar_get(nctest,"tmin32climo")
obspr = ncvar_get(nctest,"prclimo")
obspr50 = ncvar_get(nctest,"pr50climo")
nc_close(nctest)

###
# LOCA

LOCAtmax = LOCAtmax95 = LOCAtmin = LOCAtmin32 = LOCApr = LOCApr50 = list()

for(g in 1:length(GCMs)){
  
  filesin = LOCAfiles[grep(paste(GCMs[g],"_",sep=""),LOCAfiles)]
  
  if(length(filesin)>0){
    prfilein = filesin[grep("pr_",filesin)]
    tmaxfilein = filesin[grep("tasmax_",filesin)]
    tminfilein = filesin[grep("tasmin_",filesin)]
    
    if(length(prfilein)>0){
    nctest = nc_open(prfilein) # pr files
    if(g==1){
      lonLOCA = ncvar_get(nctest,"lon")
      latLOCA = ncvar_get(nctest,"lat")
    }
    LOCApr[[g]] = ncvar_get(nctest,"prclimo")
    LOCApr50[[g]] = ncvar_get(nctest,"pr50climo")
    nc_close(nctest)
    } else {
      LOCApr[[g]] = LOCApr50[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
    }
    
    if(length(tmaxfilein)>0){
    nctest = nc_open(tmaxfilein) # tmax files
    LOCAtmax[[g]] = ncvar_get(nctest,"tmaxclimo")
    LOCAtmax95[[g]] = ncvar_get(nctest,"tmax95climo")
    nc_close(nctest)
    } else {
      LOCAtmax[[g]] = LOCAtmax95[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
    }
    
    if(length(tminfilein)>0){
    nctest = nc_open(tminfilein) # tmin files
    LOCAtmin[[g]] = ncvar_get(nctest,"tminclimo")
    LOCAtmin32[[g]] = ncvar_get(nctest,"tmin32climo")
    nc_close(nctest)
    } else {
      LOCAtmin[[g]] = LOCAtmin32[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
    }
  } else {
    LOCApr[[g]] = LOCApr50[[g]] = LOCAtmax[[g]] = LOCAtmax95[[g]] = LOCAtmin[[g]] = LOCAtmin32[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
  }
}


###
# GCMs

GCMtmax = GCMtmax95 = GCMtmin = GCMtmin32 = GCMpr = GCMpr50 = list()

for(g in 1:length(GCMs)){
  
  filesin = GCMfiles[grep(paste(GCMs[g],"_",sep=""),GCMfiles)]
  
  if(length(filesin)>0){
    prfilein = filesin[grep("pr_",filesin)]
    tmaxfilein = filesin[grep("tasmax_",filesin)]
    tminfilein = filesin[grep("tasmin_",filesin)]
    
    if(length(prfilein)>0){
      nctest = nc_open(prfilein) # pr files
      if(g==1){
        lonGCM = ncvar_get(nctest,"lon")
        latGCM = ncvar_get(nctest,"lat")
      }
      GCMpr[[g]] = ncvar_get(nctest,"prclimo")
      GCMpr50[[g]] = ncvar_get(nctest,"pr50climo")
      nc_close(nctest)
    } else {
      GCMpr[[g]] = GCMpr50[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
    }
    
    if(length(tmaxfilein)>0){
      nctest = nc_open(tmaxfilein) # tmax files
      GCMtmax[[g]] = ncvar_get(nctest,"tmaxclimo")
      GCMtmax95[[g]] = ncvar_get(nctest,"tmax95climo")
      nc_close(nctest)
    } else {
      GCMtmax[[g]] = GCMtmax95[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
    }
    
    
    if(length(tminfilein)>0){
      nctest = nc_open(tminfilein) # tmin files
      GCMtmin[[g]] = ncvar_get(nctest,"tminclimo")
      GCMtmin32[[g]] = ncvar_get(nctest,"tmin32climo")
      nc_close(nctest)
    } else {
      GCMtmin[[g]] = GCMtmin32[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
    }
  }else {
    GCMpr[[g]] = GCMpr50[[g]] = GCMtmax[[g]] = GCMtmax95[[g]] = GCMtmin[[g]] = GCMtmin32[[g]] = array(NA,dim=c(length(lonobs),length(latobs),12))
  }
}

#####################
# 1) Area weighting

# calculate weights by latitude
latweights = cos((latobs*pi)/180)
latweightmat = matrix(rep(latweights,each=length(lonobs)),nrow=length(lonobs),ncol=length(latweights))
latweightvector  = as.vector(latweightmat)

###
# Weight function for later
weightfunc = function(datafile,latweightmat){
  results=array(data=NA,dim=c(nrow(datafile),ncol(datafile),12))
  for(i in 1:12){
    results[,,i]=datafile[,,i]*latweightmat
  }
  return(results)
}

###
# apply area weighting
LOCAprW = lapply(LOCApr,weightfunc,latweightmat)
LOCApr50W = lapply(LOCApr50,weightfunc,latweightmat)
LOCAtmaxW = lapply(LOCAtmax,weightfunc,latweightmat)
LOCAtmax95W = lapply(LOCAtmax95,weightfunc,latweightmat)
LOCAtminW = lapply(LOCAtmin,weightfunc,latweightmat)
LOCAtmin32W = lapply(LOCAtmin32,weightfunc,latweightmat)

GCMprW = lapply(GCMpr,weightfunc,latweightmat)
GCMpr50W = lapply(GCMpr50,weightfunc,latweightmat)
GCMtmaxW = lapply(GCMtmax,weightfunc,latweightmat)
GCMtmax95W = lapply(GCMtmax95,weightfunc,latweightmat)
GCMtminW = lapply(GCMtmin,weightfunc,latweightmat)
GCMtmin32W = lapply(GCMtmin32,weightfunc,latweightmat)

obsprW = weightfunc(obspr,latweightmat)
obspr50W = weightfunc(obspr50,latweightmat)
obstmaxW = weightfunc(obstmax,latweightmat)
obstmax95W = weightfunc(obstmax95,latweightmat)
obstminW = weightfunc(obstmin,latweightmat)
obstmin32W = weightfunc(obstmin32,latweightmat)

#############
# 2) Normalization factor

###
# calculate factors from obs

prfac = mean(apply(obsprW,c(1,2),var),na.rm=TRUE)
pr50fac = mean(apply(obspr50W,c(1,2),var),na.rm=TRUE)
tmaxfac = mean(apply(obstmaxW,c(1,2),var),na.rm=TRUE)
tmax95fac = mean(apply(obstmax95W,c(1,2),var),na.rm=TRUE)
tminfac = mean(apply(obstminW,c(1,2),var),na.rm=TRUE)
tmin32fac = mean(apply(obstmin32W,c(1,2),var),na.rm=TRUE)

###
# apply normalization factor

LOCAprnorm = lapply(LOCAprW,"/",prfac)
LOCApr50norm = lapply(LOCApr50W,"/",pr50fac)
LOCAtmaxnorm = lapply(LOCAtmaxW,"/",tmaxfac)
LOCAtmax95norm = lapply(LOCAtmax95W,"/",tmax95fac)
LOCAtminnorm = lapply(LOCAtminW,"/",tminfac)
LOCAtmin32norm = lapply(LOCAtmin32W,"/",tmin32fac)

GCMprnorm = lapply(GCMprW,"/",prfac)
GCMpr50norm = lapply(GCMpr50W,"/",pr50fac)
GCMtmaxnorm = lapply(GCMtmaxW,"/",tmaxfac)
GCMtmax95norm = lapply(GCMtmax95W,"/",tmax95fac)
GCMtminnorm = lapply(GCMtminW,"/",tminfac)
GCMtmin32norm = lapply(GCMtmin32W,"/",tmin32fac)

obsprnorm = obsprW/prfac
obspr50norm = obspr50W/pr50fac
obstmaxnorm = obstmaxW/tmaxfac
obstmax95norm = obstmax95W/tmax95fac
obstminnorm = obstminW/tminfac
obstmin32norm = obstmin32W/tmin32fac

#############
# 3) Concatenation

LOCAmat = GCMmat = matrix(NA,nrow=length(GCMs)+1,ncol=12*6*length(lonobs)*length(latobs))

LOCAtmp1 = LOCAtmp2 = LOCAtmp3 = LOCAtmp4 = LOCAtmp5 = LOCAtmp6 = matrix(NA,nrow=length(GCMs)+1,ncol=12*length(lonobs)*length(latobs))
GCMtmp1 = GCMtmp2 = GCMtmp3 = GCMtmp4 = GCMtmp5 = GCMtmp6 = matrix(NA,nrow=length(GCMs)+1,ncol=12*length(lonobs)*length(latobs))

endpoints = 1:12*length(lonobs)*length(latobs)
startpoints = c(1,endpoints+1)[1:12]

for(g in 1:(length(GCMs)+1)){
  for(m in 1:12){
    if(g<=30){
      LOCAtmp1[g,startpoints[m]:endpoints[m]] = LOCAprnorm[[g]][,,m] 
      LOCAtmp2[g,startpoints[m]:endpoints[m]] = LOCApr50norm[[g]][,,m]
      LOCAtmp3[g,startpoints[m]:endpoints[m]] = LOCAtmaxnorm[[g]][,,m]
      LOCAtmp4[g,startpoints[m]:endpoints[m]] = LOCAtmax95norm[[g]][,,m]
      LOCAtmp5[g,startpoints[m]:endpoints[m]] = LOCAtminnorm[[g]][,,m]
      LOCAtmp6[g,startpoints[m]:endpoints[m]] = LOCAtmin32norm[[g]][,,m]
    
      GCMtmp1[g,startpoints[m]:endpoints[m]] = GCMprnorm[[g]][,,m] 
      GCMtmp2[g,startpoints[m]:endpoints[m]] = GCMpr50norm[[g]][,,m]
      GCMtmp3[g,startpoints[m]:endpoints[m]] = GCMtmaxnorm[[g]][,,m]
      GCMtmp4[g,startpoints[m]:endpoints[m]] = GCMtmax95norm[[g]][,,m]
      GCMtmp5[g,startpoints[m]:endpoints[m]] = GCMtminnorm[[g]][,,m]
      GCMtmp6[g,startpoints[m]:endpoints[m]] = GCMtmin32norm[[g]][,,m]
    } else {
      LOCAtmp1[g,startpoints[m]:endpoints[m]] = obsprnorm[,,m] 
      LOCAtmp2[g,startpoints[m]:endpoints[m]] = obspr50norm[,,m]
      LOCAtmp3[g,startpoints[m]:endpoints[m]] = obstmaxnorm[,,m]
      LOCAtmp4[g,startpoints[m]:endpoints[m]] = obstmax95norm[,,m]
      LOCAtmp5[g,startpoints[m]:endpoints[m]] = obstminnorm[,,m]
      LOCAtmp6[g,startpoints[m]:endpoints[m]] = obstmin32norm[,,m]
      
      GCMtmp1[g,startpoints[m]:endpoints[m]] = obsprnorm[,,m] 
      GCMtmp2[g,startpoints[m]:endpoints[m]] = obspr50norm[,,m]
      GCMtmp3[g,startpoints[m]:endpoints[m]] = obstmaxnorm[,,m]
      GCMtmp4[g,startpoints[m]:endpoints[m]] = obstmax95norm[,,m]
      GCMtmp5[g,startpoints[m]:endpoints[m]] = obstminnorm[,,m]
      GCMtmp6[g,startpoints[m]:endpoints[m]] = obstmin32norm[,,m]
    }
  message("Finished month ",m," concatenation")
  }
  message("Finished GCM ",GCMs[g]," concatenation")
}

endpoints = 1:6*ncol(LOCAtmp1)
startpoints = c(1,endpoints+1)[1:6]

LOCAmat[,startpoints[1]:endpoints[1]] = LOCAtmp1
LOCAmat[,startpoints[2]:endpoints[2]] = LOCAtmp2
LOCAmat[,startpoints[3]:endpoints[3]] = LOCAtmp3
LOCAmat[,startpoints[4]:endpoints[4]] = LOCAtmp4
LOCAmat[,startpoints[5]:endpoints[5]] = LOCAtmp5
LOCAmat[,startpoints[6]:endpoints[6]] = LOCAtmp6

GCMmat[,startpoints[1]:endpoints[1]] = GCMtmp1
GCMmat[,startpoints[2]:endpoints[2]] = GCMtmp2
GCMmat[,startpoints[3]:endpoints[3]] = GCMtmp3
GCMmat[,startpoints[4]:endpoints[4]] = GCMtmp4
GCMmat[,startpoints[5]:endpoints[5]] = GCMtmp5
GCMmat[,startpoints[6]:endpoints[6]] = GCMtmp6

rm(list = c("LOCAtmp1","LOCAtmp2","LOCAtmp3","LOCAtmp4","LOCAtmp5","LOCAtmp6","GCMtmp1","GCMtmp2","GCMtmp3","GCMtmp4","GCMtmp5","GCMtmp6"))

obsNAs = which(is.na(LOCAmat[31,])==TRUE)
LOCAmat = LOCAmat[,-obsNAs]
GCMmat = GCMmat[,-obsNAs]

LOCAcheck = sapply(LOCAtmin,max,na.rm=TRUE)
GCMout = which(LOCAcheck== -Inf)

GCMcheck = sapply(GCMtmin,max,na.rm=TRUE)
GCMout = c(GCMout,which(GCMcheck== -Inf))
GCMout = unique(GCMout)

LOCAmat = LOCAmat[-GCMout,]
GCMmat = GCMmat[-GCMout,]

colcheck1 = apply(LOCAmat[1:26,],2,max)
colcheck1 = which(is.na(colcheck1)==TRUE)

LOCAmat = LOCAmat[,-colcheck1]

NAcheck = which(is.na(GCMmat)==TRUE,arr.ind=TRUE)
GCMmat = GCMmat[,-unique(NAcheck[,2])]
LOCAmat2 = LOCAmat[,-unique(NAcheck[,2])]

####
# anomaly calculation
GCMmat = GCMmat-apply(GCMmat,2,mean,na.rm=TRUE)
LOCAmat = LOCAmat-apply(LOCAmat,2,mean,na.rm=TRUE)
LOCAmat2 = LOCAmat2-apply(LOCAmat2,2,mean,na.rm=TRUE)

GCMs = GCMs[-GCMout]
models = c(as.character(GCMs),"OBS")
############
# Save out concatenated matrices

save(list=c("GCMmat","LOCAmat","LOCAmat2","models"),file="CPO_SKCanomalymatrices.Rdata")

