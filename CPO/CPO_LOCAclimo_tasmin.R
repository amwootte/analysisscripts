####################
#
# CPO Analysis
# Monthly Historical LOCA Climatology

library(ncdf4)
library(maps)
library(fields)
library(sp)

histperiod = c(1981,2005)
histdates = seq(as.Date("1950-01-01"),as.Date("2005-12-31"),by="day")
histmonths = seq(as.Date("1981-01-15"),as.Date("2005-12-15"),by="month")
latrange = c(26,40)
lonrange = c(251,270)

filepath = "/data4/data/DS_proj/LOCA/"

###########
#######
#  tasmin regular GCMs

filelisttmin1 = system(paste("ls ",filepath,"tasmin/historical/*historical.nc",sep=""),intern=TRUE)
split2 = do.call("rbind",strsplit(filelisttmin1,"_",fixed=TRUE))

TMINlist = list()
TMIN32list = list()

for(f in 1:length(filelisttmin1)){
  
  if(f!=20){
  nctest = nc_open(filelisttmin1[f])
  times = ncvar_get(nctest,"time")
  
  if(length(times)==length(histdates)) dates=histdates
  if(length(times)<length(histdates)) dates = histdates[-which(substr(histdates,6,10)=="02-29")]
  
  if(f==1){
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    if(lon[1]<0 & lonrange[1]>0) lonrange=lonrange-360
    if(lon[1]>0 & lonrange[1]>0) lonrange=lonrange-360;lon=lon-360;
    lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
    latidx = which(lat>= latrange[1] & lat<=latrange[2])
  }
  nc_close(nctest)
  TMIN = TMIN32 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
  
  for(h in 1:length(histmonths)){
    dateidx = which(substr(dates,1,7)==substr(histmonths[h],1,7))
    nctest = nc_open(filelisttmin1[f])
    tmin = ncvar_get(nctest,paste("tasmin_",split2[f,3],"_",split2[f,4],"_historical",sep=""),start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    if(min(tmin,na.rm=TRUE)>100){
      tmin = tmin-273.15
    } 
    nc_close(nctest)
    
    tmin32 = ifelse(tmin<=0,1,0)
    
    TMIN[,,h] = apply(tmin,c(1,2),mean,na.rm=TRUE)
    TMIN32[,,h] = apply(tmin32,c(1,2),sum,na.rm=TRUE)
    TMIN32[,,h] = ifelse(is.na(TMIN[,,h])==FALSE,TMIN32[,,h],NA)
    message("Finished calcs for histmonth ",histmonths[h]," and file ",f)
  }
  
  tminclimo = tmin32climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    tminclimo[,,month] = apply(TMIN[,,monidx],c(1,2),mean,na.rm=TRUE)
    tmin32climo[,,month] = apply(TMIN32[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and file ",f)
  }
  
  TMINlist[[f]] = tminclimo
  TMIN32list[[f]] = tmin32climo
  
  rm(TMIN)
  rm(TMIN32)
  rm(tminclimo)
  rm(tmin32climo)
  
  message("Finished calcs for file ",filelisttmin1[f])
  }
}

##########
#######
# tasmin odd GCMs

filelisttmin2 = system(paste("ls ",filepath,"tasmin/historical/*historical_1980.nc",sep=""),intern=TRUE)
split2o = do.call("rbind",strsplit(filelisttmin2,"_",fixed=TRUE))
split2o = split2o[,1:5]
split2o[,5] = "historical.nc"
startpoint = nrow(split2)
split2 = rbind(split2,split2o)
oddGCMs = split2o[,3]

for(i in 1:length(oddGCMs)){
  
  filelist = system(paste("ls ",filepath,"tasmin/historical/tasmin_",oddGCMs[i],"*historical_*.nc",sep=""),intern=TRUE)
  
  TMIN = TMIN32 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
  
  for(h in 1:length(histmonths)){
    
    fileidx = grep(substr(histmonths[h],1,4),filelist)
    dates = seq(as.Date(paste(substr(histmonths[h],1,4),"-01-01",sep="")),as.Date(paste(substr(histmonths[h],1,4),"-12-31",sep="")),by="day")
    
    nctest = nc_open(filelist[fileidx])
    times = ncvar_get(nctest,"time")
    if(length(times)<length(dates)) dates = dates[-which(substr(dates,6,10)=="02-29")]
    dateidx = which(substr(dates,6,7)==substr(histmonths[h],6,7))
    tmin = ncvar_get(nctest,paste("tasmin_",split2o[i,3],"_",split2o[i,4],"_historical",sep=""),start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    if(min(tmin,na.rm=TRUE)>100){
      tmin = tmin-273.15
    }
    nc_close(nctest)
    tmin32 = ifelse(tmin<=0,1,0)
    
    TMIN[,,h] = apply(tmin,c(1,2),mean,na.rm=TRUE)
    TMIN32[,,h] = apply(tmin32,c(1,2),sum,na.rm=TRUE)
    TMIN32[,,h] = ifelse(is.na(TMIN[,,h])==FALSE,TMIN32[,,h],NA)
    message("Finished calcs for histmonth ",histmonths[h])
  }
  
  tminclimo = tmin32climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    tminclimo[,,month] = apply(TMIN[,,monidx],c(1,2),mean,na.rm=TRUE)
    tmin32climo[,,month] = apply(TMIN32[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and file ",i)
  }
  
  TMINlist[[(startpoint+i)]] = tminclimo
  TMIN32list[[(startpoint+i)]] = tmin32climo
  
  rm(TMIN)
  rm(TMIN32)
  rm(tminclimo)
  rm(tmin32climo)
  message("Finished calcs for GCM ",oddGCMs[i])
}

###############
# netcdf creation

dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
dimY <- ncdim_def( "lat", "degrees_north", lat[latidx])
dimT <- ncdim_def("month","month",1:12)
# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

for(i in 1:nrow(split2)){
  
  if(i!=20){
  var1d <- ncvar_def("tminclimo","degrees_C", list(dimX,dimY,dimT), mv )
  var4d <- ncvar_def("tmin32climo","days", list(dimX,dimY,dimT), mv )
  
  #######
  # Create netcdf file
  
  nc <- nc_create(paste("/home/woot0002/tasmin_",split2[i,3],"_",split2[i,4],"_LOCA_histclimo.nc",sep="") ,  list(var1d,var4d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, TMINlist[[i]]) # no start or count: write all values\
  ncvar_put(nc, var4d, TMIN32list[[i]]) # no start or count: write all values\
  
  # close ncdf
  nc_close(nc)
  }
  
}

