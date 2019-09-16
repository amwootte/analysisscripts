#####################
#
# RMSE calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = "pr"

files = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*historical*.nc",sep=""),intern=TRUE)
###
# reorder to make obs data last in the files and filelist
files = files[c(2:19,1,20:21)]
filelist = do.call("rbind",strsplit(files,"/",fixed=TRUE))

####
# Extra Functions
weightfunc = function(datafile,latweightmat){
  results=array(data=NA,dim=c(nrow(datafile),ncol(datafile),12))
  for(i in 1:12){
    results[,,i]=datafile[,,i]*latweightmat
  }
  return(results)
}

RMSE = function(exp,obs){
  err = exp-obs
  errsq = err^2
  merr = mean(errsq,na.rm=TRUE)
  if(merr!=0){
    rmse = sqrt(merr)
  } else {
    rmse =0
  }
  rmse
}

####
# Apply the area weighting
prmon_latweight = list()

for(i in 1:length(files)){
  test = nc_open(files[i])
  lat = ncvar_get(test,"lat")
  lon = ncvar_get(test,"lon")
  prmon = ncvar_get(test,"prmon")
  nc_close(test)
  latweights = cos((lat*pi)/180)
  latweightmat = matrix(rep(latweights,each=length(lon)),nrow=length(lon),ncol=length(latweights))
  latweightvector = as.vector(latweightmat)
  prmon_latweight[[i]]=weightfunc(prmon,latweightmat = latweightmat)
}

####
# RMSE matrix calculation

RMSEmat = matrix(NA,nrow=length(files),ncol=(length(files)))

for(i in 1:length(files)){
  for(j in 1:length(files)){
    if(i!=j){
      RMSEmat[i,j] = RMSE(exp=as.vector(prmon_latweight[[i]]),obs=as.vector(prmon_latweight[[j]]))
    }
  }
  message("RMSE calculated for ",i," / ",length(files))
}

####
# normalize RMSE matrix

normRMSEmat = RMSEmat/mean(RMSEmat,na.rm=TRUE)

####
# save output in Rdata files

save(list=c("files","RMSEmat","normRMSEmat"),file=paste("/home/woot0002/RMSEfiles/",var,"_RMSEmats.Rdata",sep=""))



