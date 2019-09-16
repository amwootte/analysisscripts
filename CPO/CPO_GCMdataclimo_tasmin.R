####################
#
# CPO Analysis
# Monthly Historical LOCA Climatology

library(ncdf4)
library(maps)
library(fields)
library(sp)

histperiod = c(1980,2005)
histdates = seq(as.Date("1950-01-01"),as.Date("2005-12-31"),by="day")
histmonths = seq(as.Date("1980-01-15"),as.Date("2005-12-15"),by="month")
latrange = c(26,40)
lonrange = c(251,270)

########
# tasmax regular GCMs

GCMdata = read.table("tasminCMIP5grab.csv",sep=",",header=TRUE)
GCMdata$model = as.character(GCMdata$model)
GCMdata$local_file = as.character(GCMdata$local_file)

GCMs = unique(GCMdata$model)

for(g in 17:length(GCMs)){
  
  g=18
  tmp = subset(GCMdata,model==GCMs[g])
  fileranges = do.call("rbind",strsplit(as.character(tmp$time),"-"))
  for(n in 1:ncol(fileranges)) fileranges[,n] = paste(substr(fileranges[,n],1,4),"-",substr(fileranges[,n],5,6),"-",substr(fileranges[,n],7,8),sep="")
  
  for(h in 1:length(histmonths)){
    
    #histyear = as.numeric(substr(histmonths[h],1,4))
    fileidx = which(as.Date(fileranges[,1])<= as.Date(histmonths[h]) & as.Date(fileranges[,2])>=as.Date(histmonths[h]))
    
    if(length(fileidx)>0){
    if(length(fileidx)>1) fileidx = fileidx[1]
    
    moddates = seq(as.Date(fileranges[fileidx,1]),as.Date(fileranges[fileidx,2]),by="day")
    
    nctest = nc_open(tmp$local_file[fileidx])
    times = ncvar_get(nctest,"time")
    
    if(length(times)==length(moddates)){
      dates=histdates
    } 
    if(length(times)<length(moddates)){
      moddates = moddates[-which(substr(moddates,6,10)=="02-29")]
      dates = histdates[-which(substr(histdates,6,10)=="02-29")]
    } 
    
    if(h==1){
      lon = ncvar_get(nctest,"lon")
      lat = ncvar_get(nctest,"lat")
      if(lon[1]<0 & lonrange[1]>0) lonrange=lonrange-360
      if(lon[1]>=0 & lonrange[1]>=0) lonrange=lonrange-360;lon=lon-360;
      lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
      latidx = which(lat>= latrange[1] & lat<=latrange[2])
    }
    
    dateidx = which(substr(moddates,1,7)==substr(histmonths[h],1,7))
    tmin = ncvar_get(nctest,"tasmin",start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))-273.15
    
    nc_close(nctest)
    tmin32 = ifelse(tmin<=0,1,0)
    }
    
    if(h==1) TMIN = TMIN32 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
    
    if(length(fileidx)>0){
    TMIN[,,h] = apply(tmin,c(1,2),mean,na.rm=TRUE)
    TMIN32[,,h] = apply(tmin32,c(1,2),sum,na.rm=TRUE)
    TMIN32[,,h] = ifelse(is.na(TMIN[,,h])==FALSE,TMIN32[,,h],NA)
    }
    message("Finished calcs for histmonth ",histmonths[h]," and GCM ",GCMs[g])
  }
  
  
  tminclimo = tmin32climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    tminclimo[,,month] = apply(TMIN[,,monidx],c(1,2),mean,na.rm=TRUE)
    tmin32climo[,,month] = apply(TMIN32[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and GCM ",GCMs[g])
  }
  
  dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
  dimY <- ncdim_def( "lat", "degrees_north",lat[latidx])
  dimT <- ncdim_def("month","month",1:12)
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def("tminclimo","degrees_C", list(dimX,dimY,dimT), mv )
  var4d <- ncvar_def("tmin32climo","days", list(dimX,dimY,dimT), mv )
  
  #######
  # Create netcdf file
  ens = as.character(unique(GCMdata[which(GCMdata$model==GCMs[g]),8]))
  nc <- nc_create(paste("/home/woot0002/tasmin_",GCMs[g],"_",ens,"_CMIP5_histclimo.nc",sep="") ,  list(var1d,var4d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, tminclimo) # no start or count: write all values\
  ncvar_put(nc, var4d, tmin32climo) # no start or count: write all values\
  
  # close ncdf
  nc_close(nc)
  
  message("Finished calcuations for GCM ",GCMs[g])
}

###########
# make tasmax netcdfs for LOCA
1:length(GCMs)

for(i in 21:length(GCMs)){

  dimX <- ncdim_def( "lon", "degrees_east", LONlist[[i]])
  dimY <- ncdim_def( "lat", "degrees_north",LATlist[[i]])
  dimT <- ncdim_def("month","month",1:12)
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  
var1d <- ncvar_def("tmaxclimo","degrees_C", list(dimX,dimY,dimT), mv )
var4d <- ncvar_def("tmax95climo","days", list(dimX,dimY,dimT), mv )

#######
# Create netcdf file
ens = as.character(unique(GCMdata[which(GCMdata$model==GCMs[i]),8]))
nc <- nc_create(paste("/home/woot0002/tasmin_",GCMs[i],"_",ens,"_CMIP5_histclimo.nc",sep="") ,  list(var1d,var4d) )

# Write some data to the file
ncvar_put(nc, var1d, TMAXlist[[i]]) # no start or count: write all values\
ncvar_put(nc, var4d, TMAX95list[[i]]) # no start or count: write all values\

# close ncdf
nc_close(nc)

}



