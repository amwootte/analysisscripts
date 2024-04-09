rm(list=ls(all=TRUE))
gc()

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
library(parallel)
library(foreach)
library(doParallel)


#####
# Version of this script to run parallel but also one domain as a time.

var = "pr"
enssize = 5
ensgroup = "CMIP6"
mask = "SGP-NCA"
climotype = "annual"
inputfile = paste("/home/woot0002/RCMES/",var,"_",ensgroup,"_climo_",mask,"_",climotype,".nc",sep="")



#####

nctest = nc_open(inputfile)
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")

if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
varnames = c()
outputdata = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnames[i]=nctest$var[[i]]$name
  outputdata[,,i]=ncvar_get(nctest,varnames[i])
}
} 

if(climotype=="monthly"){
  varnames = c()
  outputdata = array(NA,dim=c(length(lon),length(lat),12,length(nctest$var)))
  for(i in 1:length(nctest$var)){
    varnames[i]=nctest$var[[i]]$name
    outputdata[,,,i]=ncvar_get(nctest,varnames[i])
  }
} 

varunits = nctest$var[[1]]$units

dimunits=c()
for(i in 1:length(nctest$dim)){
  dimunits[i]=nctest$dim[[i]]$units
}

nc_close(nctest)

#####
# RMSE calcs 

idxuse = 2:(length(varnames)-1)

GCMnamelist = list()
RMSEout = NULL

combosused = combn(idxuse,enssize,simplify=TRUE)

if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
  OBS = outputdata[,,1]
}
if(climotype=="monthly"){
  OBS = outputdata[,,,1]
}

#RMSEcalculator = function(arraydat,combosused,i,OBS){
#  meanval = apply(arraydat[,,combosused[,i]],c(1,2),mean,na.rm=TRUE)
#  sqrt(mean((meanval-OBS)^2,na.rm=TRUE))
#}
  
  RMSEcalculator = function(arraydat,combosused,i,OBS,climotype="annual"){
    
    if(climotype=="annual"  | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
      testarray = apply(arraydat[,,combosused[,i]],3,c)
    }
    if(climotype=="monthly"){
      testarray = apply(arraydat[,,,combosused[,i]],4,c)
    }
    meanval = rowMeans(testarray,na.rm=TRUE)
    OBSarray = as.vector(OBS)
    sqrt(mean((meanval-OBS)^2,na.rm=TRUE))
  }
  

  idxs = 1:ncol(combosused)

library(parallel)
#system.time({
  n_cores <- detectCores(logical=FALSE)
  cl <- makeCluster(n_cores-1, type = "FORK",methods=FALSE,useXDR=FALSE)
  #RMSEs = c()
  #for(i in 1:10){
  #  meanval=parApply(cl, outputdata[,,combosused[,i]], c(1,2), mean,na.rm=TRUE)
  #  RMSEs[i]=sqrt(mean((meanval-OBS)^2,na.rm=TRUE))
  #}
  
  RMSEs <- parLapply(cl , idxs,RMSEcalculator,arraydat=outputdata,combosused=combosused,OBS,climotype=climotype) # dimnum=3 for climotype annual, dimnum=4 for climotype monthly
#})
stopCluster(cl)

RMSEs = do.call("c",RMSEs)

NRMSEs=RMSEs/sd(OBS,na.rm=TRUE)

group = 1:ncol(combosused)
#RMSEout  = data.frame(enssize,group,mask,RMSEs,NRMSEs)
#message("Finished calcs for ",mask)

fileout = paste("/home/woot0002/RCMES/",var,"_",ensgroup,"_RMSEdat_enssize",enssize,"_dom",mask,"_",climotype,"_reduced.Rdata",sep="")
save(list=c("varnames","idxs","combosused","RMSEs","NRMSEs","group","enssize"),file=fileout)

#save(list=c("idxs"),file=paste("RCMES/",var,"_LOCA_RMSEdat_enssize",enssize,"_dom",mask,"_idxs.Rdata",sep=""))

