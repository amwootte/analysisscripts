source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

load("/home/woot0002/DS_ind/GCMlist.Rdata")

varin="pr"

diff1=diff2 = var1old = var2old = var1new = var2new= list()

for(i in 1:length(GCMlist)){
  
  GCMfiles_old = system(paste("ls /home/woot0002/GCMs/regrid/notuse2/",varin,"*",GCMlist[i],"_*histclimo*.nc",sep=""),intern=TRUE)
  GCMfiles_new = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"*",GCMlist[i],"_*histclimo*.nc",sep=""),intern=TRUE)
  if(varin=="pr"){
  ncvarname1 = "prclimo"
  ncvarname2 = "pr50climo"
  }
  if(varin=="tmax"){
    ncvarname1 = "tasmaxclimo"
    ncvarname2 = "tmax95climo"
  }
  if(varin=="tmin"){
    ncvarname1 = "tasminclimo"
    ncvarname2 = "tmin32climo"
  }
  
    nctest = nc_open(GCMfiles_old)
    var1old[[i]] = ncvar_get(nctest,ncvarname1)
    var2old[[i]] = ncvar_get(nctest,ncvarname2)
    if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
    nc_close(nctest)
    
    nctest = nc_open(GCMfiles_new)
    var1new[[i]] = ncvar_get(nctest,ncvarname1)
    var2new[[i]] = ncvar_get(nctest,ncvarname2)
    nc_close(nctest)
    
    diff1[[i]]=var1new[[i]]-var1old[[i]]
    diff2[[i]]=var2new[[i]]-var2old[[i]]
    
}

sapply(diff1,range,na.rm=TRUE)
sapply(diff2,range,na.rm=TRUE)

testsfc = list(x=lon,y=lat,z=diff1[[18]][,,7])
surface(testsfc,type="I")

testsfc = list(x=lon,y=lat,z=var1old[[18]][,,7])
surface(testsfc,type="I")

testsfc = list(x=lon,y=lat,z=var1new[[18]][,,7])
surface(testsfc,type="I")

