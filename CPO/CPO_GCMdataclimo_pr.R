####################
#
# CPO Analysis
# Monthly Historical LOCA Climatology

library(ncdf4)
library(maps)
library(fields)
library(sp)
library(PCICt)

histperiod = c(1980,2005)
histdates = seq(as.Date("1950-01-01"),as.Date("2005-12-31"),by="day")
histmonths = seq(as.Date("1980-01-15"),as.Date("2005-12-15"),by="month")
latrange = c(20,45)
lonrange = c(245,275)

########
# tasmax regular GCMs

GCMdata = read.table("/home/woot0002/old_csv/prCMIP5grab.csv",sep=",",header=TRUE)
GCMdata$model = as.character(GCMdata$model)
GCMdata$local_file = as.character(GCMdata$local_file)

GCMs = unique(GCMdata$model)

for(g in 1:length(GCMs)){
  
  #g=18
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
      if(length(which(substr(moddates,6,10)=="02-29"))>=1){
        moddates = moddates[-which(substr(moddates,6,10)=="02-29")]
        dates = histdates[-which(substr(histdates,6,10)=="02-29")]
      }
    } 
    if(length(times)<length(moddates)){
      modmonyear = unique(substr(moddates,1,7))
      histmonyear = unique(substr(histdates,1,7))
      daysin = as.character(1:30)
      daysin[1:9] = paste("0",daysin[1:9],sep="")
      
      moddates2 = c()
      for(i in 1:length(modmonyear)){
        moddates2 = c(moddates2,paste(modmonyear[i],daysin,sep="-"))
      }
      
      dates2 = c()
      for(i in 1:length(histmonyear)){
        dates2 = c(dates2,paste(histmonyear[i],daysin,sep="-"))
      }
      moddates = moddates2
      dates = dates2
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
    pr = ncvar_get(nctest,"pr",start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))*86400
    #pr = ncvar_get(nctest,"pr")*86400
    
    nc_close(nctest)
    pr50 = ifelse(pr>=50.4,1,0)
    }
    
    if(h==1) PR = PR50 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
    #if(h==1) PR = PR50 = array(NA,dim=c(length(lon),length(lat),length(histmonths)))
    
    if(length(fileidx)>0){
    PR[,,h] = apply(pr,c(1,2),sum,na.rm=TRUE)
    PR50[,,h] = apply(pr50,c(1,2),sum,na.rm=TRUE)
    PR50[,,h] = ifelse(is.na(PR[,,h])==FALSE,PR50[,,h],NA)
    }
    message("Finished calcs for histmonth ",histmonths[h]," and GCM ",GCMs[g])
  }
  
  
  prclimo = pr50climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  #prclimo = pr50climo = array(NA,dim=c(length(lon),length(lat),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    prclimo[,,month] = apply(PR[,,monidx],c(1,2),mean,na.rm=TRUE)
    pr50climo[,,month] = apply(PR50[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and GCM ",GCMs[g])
  }
  
  dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
  #dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north",lat[latidx])
  #dimY <- ncdim_def( "lat", "degrees_north",lat)
  dimT <- ncdim_def("month","month",1:12)
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def("prclimo","", list(dimX,dimY,dimT), mv )
  var4d <- ncvar_def("pr50climo","days", list(dimX,dimY,dimT), mv )
  
  #######
  # Create netcdf file
  ens = as.character(unique(GCMdata[which(GCMdata$model==GCMs[g]),8]))
  nc <- nc_create(paste("/home/woot0002/GCMs/pr_",GCMs[g],"_",ens,"_CMIP5_histclimo.nc",sep="") ,  list(var1d,var4d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, prclimo) # no start or count: write all values\
  ncvar_put(nc, var4d, pr50climo) # no start or count: write all values\
  
  # close ncdf
  nc_close(nc)
  
  message("Finished calcuations for GCM ",GCMs[g])
}



