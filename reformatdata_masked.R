#######################
#
# Data reformatting

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = varin = "pr"
models = "LOCA"
period = "hist"
statemask = "louisiana"

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

test = nc_open(paste("/home/woot0002/DS_ind/",statemask,"_mask.nc",sep=""))
regionmask = ncvar_get(test,"mask")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

if(models=="GCMs" | models=="LOCA"){

  if(period=="hist"){
    files = system(paste("ls /home/woot0002/",models,"/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)
  }
  if(period=="proj"){
    files = system(paste("ls /home/woot0002/",models,"/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)
  }
load("/home/woot0002/DS_ind/manuscript1/GCMlist.Rdata")

hfiles = c()

for(i in 1:length(GCMlist)){
  hfiles[i] = files[grep(paste(GCMlist[i],"_",sep=""),files)]
}

ncvarname = paste(var,"climo",sep="")

### File data
for(i in 1:length(GCMlist)){
  nctest = nc_open(hfiles[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  vardata= ncvar_get(nctest,nctest$var[[idx]]$name)
  nc_close(nctest)
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) hvarmat = matrix(NA,ncol=length(GCMlist),nrow=dim(vardata)[1]*dim(vardata)[2])
  if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
    vardat = apply(vardata,c(1,2),sum,na.rm=TRUE)
    vardat = ifelse(regionmask==1,vardat,NA)
  } 
  if(var=="tmax" | var=="tmin"){
    vardat = apply(vardata,c(1,2),mean,na.rm=TRUE)
    vardat = ifelse(regionmask==1,vardat,NA)
  } 
  
  hvarmat[,i]=vardat
  message("Finished data pull for GCM ",GCMlist[i])
}

hvarmat = data.frame(hvarmat)
names(hvarmat) = GCMlist

if(period=="hist"){
  write.table(hvarmat,paste("/home/woot0002/DS_ind/",models,"_hist_",var,"_",statemask,".csv",sep=""),sep=",",row.names=FALSE)
}
if(period=="proj"){
  write.table(hvarmat,paste("/home/woot0002/DS_ind/",models,"_proj_",var,"_",statemask,".csv",sep=""),sep=",",row.names=FALSE)
}

} 

if(models=="OBS"){
  
  mask = read.table("/home/woot0002/DS_ind/LOCA_hist_tmax.csv",sep=",",header=TRUE)
  
  
  hfiles = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)
  
  ncvarname = paste(var,"climo",sep="")
  
  ### File data
    nctest = nc_open(hfiles)
    idx = which(names(nctest$var)==ncvarname)
    if(length(idx)==0){
      idx = which(names(nctest$var)==var)
      if(length(idx)==0){
        idx = which(names(nctest$var)==varin)
      }
    }
    vardata= ncvar_get(nctest,nctest$var[[idx]]$name)
    nc_close(nctest)
    #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
    hvarmat = matrix(NA,ncol=1,nrow=dim(vardata)[1]*dim(vardata)[2])
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      vardat = apply(vardata,c(1,2),sum,na.rm=TRUE)
      vardat = ifelse(regionmask==1,vardat,NA)
    } 
    if(var=="tmax" | var=="tmin"){
      vardat = apply(vardata,c(1,2),mean,na.rm=TRUE)
      vardat = ifelse(regionmask==1,vardat,NA)
    } 
    
    hvarmat[,1]=vardat
    hvarmat[,1] = ifelse(regionmask==1,hvarmat[,1],NA)
  hvarmat = data.frame(hvarmat)
  names(hvarmat) = "LIVNEH"
  
  write.table(hvarmat,paste("/home/woot0002/DS_ind/LIVNEH_hist_",var,"_",statemask,".csv",sep=""),sep=",",row.names=FALSE)
  
}