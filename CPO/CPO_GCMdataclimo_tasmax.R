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

GCMdata = read.table("tasmaxCMIP5grab.csv",sep=",",header=TRUE)
GCMdata$model = as.character(GCMdata$model)
GCMdata$local_file = as.character(GCMdata$local_file)

GCMs = unique(GCMdata$model)

TMAXlist = list()
TMAX95list = list()
LONlist = list()
LATlist = list()

for(g in 1:length(GCMs)){
  
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
      LONlist[[g]]=lon[lonidx]
      LATlist[[g]]=lat[latidx]
    }
    
    dateidx = which(substr(moddates,1,7)==substr(histmonths[h],1,7))
    tmax = ncvar_get(nctest,"tasmax",start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))-273.15
    
    nc_close(nctest)
    tmax95 = ifelse(tmax>=35,1,0)
    }
    
    if(h==1) TMAX = TMAX95 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
    
    if(length(fileidx)>0){
    TMAX[,,h] = apply(tmax,c(1,2),mean,na.rm=TRUE)
    TMAX95[,,h] = apply(tmax95,c(1,2),sum,na.rm=TRUE)
    TMAX95[,,h] = ifelse(is.na(TMAX[,,h])==FALSE,TMAX95[,,h],NA)
    }
    message("Finished calcs for histmonth ",histmonths[h]," and GCM ",GCMs[g])
  }
  
  
  tmaxclimo = tmax95climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    tmaxclimo[,,month] = apply(TMAX[,,monidx],c(1,2),mean,na.rm=TRUE)
    tmax95climo[,,month] = apply(TMAX95[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and GCM ",GCMs[g])
  }
  
  
  TMAXlist[[g]]=tmaxclimo
  TMAX95list[[g]]=tmax95climo
  
  message("Finished gathering data for GCM ",GCMs[g])
  
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
nc <- nc_create(paste("/home/woot0002/tasmax_",GCMs[i],"_",ens,"_CMIP5_histclimo.nc",sep="") ,  list(var1d,var4d) )

# Write some data to the file
ncvar_put(nc, var1d, TMAXlist[[i]]) # no start or count: write all values\
ncvar_put(nc, var4d, TMAX95list[[i]]) # no start or count: write all values\

# close ncdf
nc_close(nc)

}



