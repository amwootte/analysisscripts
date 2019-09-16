
library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("/home/woot0002/analysisfunctions.R")

###########
# 1. Data Gather and conversion

  tmaxhistfilelist = system("ls /data4/data/OBS/METDATA/tasmax/tmmx_*.nc",intern=T)
  tminhistfilelist = system("ls /data4/data/OBS/METDATA/tasmin/tmmn_*.nc",intern=T)
 
dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
years = 1981:2005

for(i in 1:length(tmaxhistfilelist)){
  ptm = proc.time()
  message("Starting work on file ",tmaxhistfilelist[i])
  
  yearidx = which(as.numeric(substr(dates,1,4))==years[i])
  
  test1 = nc_open(tmaxhistfilelist[i])
  if(i==1){
    lon = ncvar_get(test1,"lon")
    lat = ncvar_get(test1,"lat")
    tmaxarray = tminarray = array(NA,dim=c(length(lon),length(lat),length(dates)))
  }
  test2 = nc_open(tminhistfilelist[i])
  tmaxarray[,,yearidx] = ncvar_get(test1,"air_temperature")
  tminarray[,,yearidx] = ncvar_get(test2,"air_temperature")
  nc_close(test1)
  nc_close(test2)
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(tmaxhistfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}  
  
  
  
  
  
  message("starting calcs for q95")
  for(R in 1:length(lon)){
    for(C in 1:length(lat)){
      tmp1 = ncvar_get(test1,"tasmax",start = c(R,C,1),count=c(1,1,-1))
      tmp2 = ncvar_get(test2,"tasmin",start = c(R,C,1),count=c(1,1,-1))
      if(all(is.na(tmp1)==TRUE)==FALSE){
        q95max[R,C] = quantile(tmp1,probs=0.95,na.rm=TRUE)
        q95min[R,C] = quantile(tmp2,probs=0.95,na.rm=TRUE)
      } 
    message("Finished Calcs for R ",R," and C ",C)
      }
  }
  
 
  filesplit = strsplit(tmaxhistfilelist[i],"/")[[1]]
  newq95maxfile = paste(substr(filesplit[length(filesplit)],1,(nchar(filesplit[length(filesplit)])-3)),"_q95.nc",sep="")
  filesplit = strsplit(tminhistfilelist[i],"/")[[1]]
  newq95minfile = paste(substr(filesplit[length(filesplit)],1,(nchar(filesplit[length(filesplit)])-3)),"_q95.nc",sep="")
  
  londef = ncdim_def("lon",units="degrees E",vals=lon,longname="Longitude")
  latdef = ncdim_def("lat",units="degrees N",vals=lat,longname="Latitude")
  
  q95maxdef = ncvar_def("tmaxq95",units="degrees K",dim=list(londef,latdef),missval=1E20,longname="95th Percentile of Daily High Temperatures")
  q95mindef = ncvar_def("tminq95",units="degrees K",dim=list(londef,latdef),missval=1E20,longname="95th Percentile of Daily Low Temperatures")
  
  nc1 = nc_create(newq95maxfile,q95maxdef)
  nc2 = nc_create(newq95minfile,q95mindef)
  
  ncvar_put(nc1,q95maxdef,q95max)
  ncvar_put(nc2,q95mindef,q95min)
  
  nc_close(nc1)
  nc_close(nc2)
  
  gc()
 
}
