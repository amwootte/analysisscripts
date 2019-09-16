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

########
# tasmax regular GCMs

filelisttmax1 = system(paste("ls ",filepath,"tasmax/historical/*historical.nc",sep=""),intern=TRUE)
split1 = do.call("rbind",strsplit(filelisttmax1,"_",fixed=TRUE))

TMAXlist = list()
TMAX95list = list()

for(f in 1:length(filelisttmax1)){
  nctest = nc_open(filelisttmax1[f])
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
  TMAX = TMAX95 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
  
  for(h in 1:length(histmonths)){
    dateidx = which(substr(dates,1,7)==substr(histmonths[h],1,7))
    nctest = nc_open(filelisttmax1[f])
    tmax = ncvar_get(nctest,paste("tasmax_",split1[f,3],"_",split1[f,4],"_historical",sep=""),start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    if(min(tmax,na.rm=TRUE)>100){
      tmax = tmax-273.15
    }
    nc_close(nctest)
    
    tmax95 = ifelse(tmax>=35,1,0)
    
    TMAX[,,h] = apply(tmax,c(1,2),mean,na.rm=TRUE)
    TMAX95[,,h] = apply(tmax95,c(1,2),sum,na.rm=TRUE)
    TMAX95[,,h] = ifelse(is.na(TMAX[,,h])==FALSE,TMAX95[,,h],NA)
    message("Finished calcs for histmonth ",histmonths[h]," and file ",f)
  }
  
  tmaxclimo = tmax95climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    tmaxclimo[,,month] = apply(TMAX[,,monidx],c(1,2),mean,na.rm=TRUE)
    tmax95climo[,,month] = apply(TMAX95[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and file ",f)
  }
  
  TMAXlist[[f]] = tmaxclimo
  TMAX95list[[f]] = tmax95climo
  
  rm(TMAX)
  rm(TMAX95)
  rm(tmaxclimo)
  rm(tmax95climo)
  
  message("Finished calcs for file ",filelisttmax1[f])
  
}

##########
#######
# tasmax odd GCMs

oddGCM = "MRI-CGCM3"
filelisttmax2 = system(paste("ls ",filepath,"tasmax/historical/*historical_*.nc",sep=""),intern=TRUE)
split1o = strsplit(filelisttmax2,"_",fixed=TRUE)
split1o = c(split1o[[1]][1],split1o[[1]][2],split1o[[1]][3],split1o[[1]][4],"historical.nc")
split1 = rbind(split1,split1o)

TMAX = TMAX95 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))

for(h in 1:length(histmonths)){
  
  fileidx = grep(substr(histmonths[h],1,4),filelisttmax2)
  dates = seq(as.Date(paste(substr(histmonths[h],1,4),"-01-01",sep="")),as.Date(paste(substr(histmonths[h],1,4),"-12-31",sep="")),by="day")
  
  nctest = nc_open(filelisttmax2[fileidx])
  times = ncvar_get(nctest,"time")
  if(length(times)<length(dates)) dates = dates[-which(substr(dates,6,10)=="02-29")]
  dateidx = which(substr(dates,6,7)==substr(histmonths[h],6,7))
  tmax = ncvar_get(nctest,paste("tasmax_",split1o[3],"_",split1o[4],"_historical",sep=""),start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
  if(min(tmax,na.rm=TRUE)>100){
    tmax = tmax-273.15
  }
  nc_close(nctest)
  tmax95 = ifelse(tmax>=35,1,0)
  
  TMAX[,,h] = apply(tmax,c(1,2),mean,na.rm=TRUE)
  TMAX95[,,h] = apply(tmax95,c(1,2),sum,na.rm=TRUE)
  TMAX95[,,h] = ifelse(is.na(TMAX[,,h])==FALSE,TMAX95[,,h],NA)
  message("Finished calcs for histmonth ",histmonths[h])
  }

tmaxclimo = tmax95climo = array(NA,dim=c(length(lonidx),length(latidx),12))
for(month in 1:12){
  monidx = which(as.numeric(substr(histmonths,6,7))==month)
  tmaxclimo[,,month] = apply(TMAX[,,monidx],c(1,2),mean,na.rm=TRUE)
  tmax95climo[,,month] = apply(TMAX95[,,monidx],c(1,2),mean,na.rm=TRUE)
  message("Finished climo calcs for month ",month," and file ",f)
}

TMAXlist[[nrow(split1)]] = tmaxclimo
TMAX95list[[nrow(split1)]] = tmax95climo

rm(TMAX)
rm(TMAX95)
rm(tmaxclimo)
rm(tmax95climo)

###########
# make tasmax netcdfs for LOCA

dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
dimY <- ncdim_def( "lat", "degrees_north", lat[latidx])
dimT <- ncdim_def("month","month",1:12)
# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

for(i in 1:nrow(split1)){

var1d <- ncvar_def("tmaxclimo","degrees_C", list(dimX,dimY,dimT), mv )
var4d <- ncvar_def("tmax95climo","degrees_C", list(dimX,dimY,dimT), mv )

#######
# Create netcdf file

nc <- nc_create(paste("/home/woot0002/tasmax_",split1[i,3],"_",split1[i,4],"_LOCA_histclimo.nc",sep="") ,  list(var1d,var4d) )

# Write some data to the file
ncvar_put(nc, var1d, TMAXlist[[i]]) # no start or count: write all values\
ncvar_put(nc, var4d, TMAX95list[[i]]) # no start or count: write all values\

# close ncdf
nc_close(nc)

}



###########
#######
#  tasmin regular GCMs

filelisttmin1 = system(paste("ls ",filepath,"tasmin/historical/*historical.nc",sep=""),intern=TRUE)
split2 = do.call("rbind",strsplit(filelisttmin1,"_",fixed=TRUE))

TMINlist = list()
TMIN32list = list()

for(f in 1:length(filelisttmin1)){
  
  #f=20
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

1:nrow(split2)
1:19

for(i in 21:nrow(split2)){
  
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



###########
#######
#  pr all GCMs

filelistpr1 = system(paste("ls ",filepath,"pr/historical/*historical.nc",sep=""),intern=TRUE)
split3 = do.call("rbind",strsplit(filelistpr1,"_",fixed=TRUE))

PRlist = list()
PR50list = list()

for(f in 1:length(filelistpr1)){
  nctest = nc_open(filelistpr1[f])
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
  PR = PR50 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
  
  for(h in 1:length(histmonths)){
    dateidx = which(substr(dates,1,7)==substr(histmonths[h],1,7))
    nctest = nc_open(filelistpr1[f])
    pr = ncvar_get(nctest,paste("pr_",split3[f,3],"_",split3[f,4],"_historical",sep=""),start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    if(max(pr,na.rm=TRUE)<1){
      pr = pr*86400
    }
    nc_close(nctest)
    
    pr50 = ifelse(pr>=50.4,1,0)
    
    PR[,,h] = apply(pr,c(1,2),sum,na.rm=TRUE)
    PR[,,h] = ifelse(is.na(pr[,,1])==FALSE,PR[,,h],NA)
    PR50[,,h] = apply(pr50,c(1,2),sum,na.rm=TRUE)
    PR50[,,h] = ifelse(is.na(pr[,,1])==FALSE,PR50[,,h],NA)
    message("Finished calcs for histmonth ",histmonths[h]," and file ",f)
  }
  
  prclimo = pr50climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    prclimo[,,month] = apply(PR[,,monidx],c(1,2),mean,na.rm=TRUE)
    pr50climo[,,month] = apply(PR50[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and file ",f)
  }
  
  PRlist[[f]] = prclimo
  PR50list[[f]] = pr50climo
  
  rm(PR)
  rm(PR50)
  rm(prclimo)
  rm(pr50climo)
  
  message("Finished calcs for file ",filelistpr1[f])
  
}

###############
# netcdf creation

dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
dimY <- ncdim_def( "lat", "degrees_north", lat[latidx])
dimT <- ncdim_def("month","month",1:12)
# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

for(i in 1:nrow(split3)){
  
  var1d <- ncvar_def("prclimo","mm", list(dimX,dimY,dimT), mv )
  var4d <- ncvar_def("pr50climo","days", list(dimX,dimY,dimT), mv )
  
  #######
  # Create netcdf file
  
  nc <- nc_create(paste("/home/woot0002/pr_",split3[i,3],"_",split3[i,4],"_LOCA_histclimo.nc",sep="") ,  list(var1d,var4d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, PRlist[[i]]) # no start or count: write all values\
  ncvar_put(nc, var4d, PR50list[[i]]) # no start or count: write all values\
  
  # close ncdf
  nc_close(nc)
  
}
