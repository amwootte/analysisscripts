#######################
#
# Monthly climo calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = "pr"

files = system(paste("ls /data4/data/DS_proj/MACA_LIVNEH/",var,"/historical/*.nc",sep=""),intern=TRUE)

filelist = do.call("rbind",strsplit(files,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist[,8],"_",fixed=TRUE))
names(filelist2) = c("MACAv","var","GCM","exp","period","startyear","endyear","domain","ending")
filelist2[,6]=as.numeric(filelist2[,6])
filelist2[,7]=as.numeric(filelist2[,7])

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
monyear = dates[which(substr(dates,9,10)=="01")]
histmonths = seq(as.Date("1981-01-15"),as.Date("2005-12-15"),by="month")
startyear = 1981
endyear = 2005
years = startyear:endyear
months = 1:12
latrange = c(26,40)
lonrange = c(251,270)
if(var=="tasmax" | var=="tasmin") varin = "air_temperature"
if(var=="pr") varin = "precipitation"


GCMs = unique(filelist2[,3])
idxout = which(as.numeric(filelist2[,6])<startyear & as.numeric(filelist2[,7])<startyear)

filelist3 = filelist2[-idxout,]
filesin = files[-idxout]

for(i in 1:length(GCMs)){
  
  useidx = which(filelist3[,3]==GCMs[i])
  usefiles = filesin[useidx]
  usefileinfo = filelist3[useidx,]
  
  if(i==1){
  nctest = nc_open(usefiles[1])
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    if(lon[1]<0 & lonrange[1]>0) lonrange=lonrange-360
    if(lon[1]>0 & lonrange[1]>0) lonrange=lonrange-360;lon=lon-360;
    lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
    latidx = which(lat>= latrange[1] & lat<=latrange[2])
    nc_close(nctest)
  }
 
  datayears = as.numeric(usefileinfo[1,6]):as.numeric(usefileinfo[2,7])
  outputarray1 = outputarray2 = array(NA,dim=c(length(lonidx),length(latidx),length(datayears)*12))
  if(var=="pr"){
    outputarray3 = outputarray4 = outputarray5 = array(NA,dim=c(length(lonidx),length(latidx),length(datayears)*12))
  }
  
  for(t in 1:length(datayears)){
    
    if(datayears[t]<1990){
      fileidx = 1
    } else{
      fileidx = 2
    }
    
    nctest = nc_open(usefiles[fileidx])
    times = ncvar_get(nctest,"time")
    timeunits = nctest$dim$time$units
    filestartdate = as.Date(substr(timeunits,12,21))
    filedates = filestartdate+times
    
    dateidx = which(as.numeric(substr(filedates,1,4))==datayears[t])
    
    tmax = ncvar_get(nctest,varin,start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    nc_close(nctest)
    
    if(var=="tasmax"){
    if(min(tmax,na.rm=TRUE)>100){
      var1 = tmax-273.15
    }
    var2 = ifelse(var1>=35,1,0)
    }
    
    if(var=="tasmin"){
      if(min(tmax,na.rm=TRUE)>100){
        var1 = tmax-273.15
      }
      var2 = ifelse(var1<=0,1,0)
    }
    if(var=="pr"){
      if(max(tmax,na.rm=TRUE)<1){
        var1 = tmax*86400
      } else {
        var1 = tmax
      }
      var2 = ifelse(var1>=50.4,1,0)
      var3 = ifelse(var1>=1,1,0)
    }
   
    var1tmp = var2tmp = array(NA,dim=c(length(lonidx),length(latidx),12))
    if(var=="pr"){
      var3tmp = var4tmp = var5tmp = array(NA,dim=c(length(lonidx),length(latidx),12))
    }
    datesused = filedates[dateidx]
    for(m in 1:12){
      monidx = which(as.numeric(substr(datesused,6,7))==m)
      if(var=="tasmax" | var=="tasmin"){
        var1tmp[,,m] = apply(var1[,,monidx],c(1,2),mean,na.rm=TRUE)
        var2tmp[,,m] = apply(var2[,,monidx],c(1,2),sum,na.rm=TRUE)
        var2tmp[,,m] = ifelse(is.na(var1[,,1])==FALSE,var2tmp[,,m],NA)
      }
      if(var=="pr"){
        var1tmp[,,m] = apply(var1[,,monidx],c(1,2),sum,na.rm=TRUE)
        var2tmp[,,m] = apply(var2[,,monidx],c(1,2),sum,na.rm=TRUE)
        var2tmp[,,m] = ifelse(is.na(var1[,,1])==FALSE,var2tmp[,,m],NA)
        var3tmp[,,m] = apply(var3[,,monidx],c(1,2),sum,na.rm=TRUE)
        var3tmp[,,m] = ifelse(is.na(var1[,,1])==FALSE,var3tmp[,,m],NA)
        var4tmp[,,m] = apply(var1[,,monidx],c(1,2),max,na.rm=TRUE)
        var4tmp[,,m] = ifelse(is.na(var1[,,1])==FALSE,var4tmp[,,m],NA)
        var5tmp[,,m] = apply(var1[,,monidx],c(1,2),calcrollsum,size=5)
        var5tmp[,,m] = ifelse(is.na(var1[,,1])==FALSE,var5tmp[,,m],NA)
      }
    }
    
    idxfileout = 1:12+(12*(t-1))
    outputarray1[,,idxfileout]=var1tmp
    outputarray2[,,idxfileout]=var2tmp
    if(var=="pr"){
      outputarray3[,,idxfileout]=var3tmp
      outputarray4[,,idxfileout]=var4tmp
      outputarray5[,,idxfileout]=var5tmp
    }
    message("Finished calcs for ",datayears[t])
  }
  
  histmonthsused = seq(as.Date(paste(datayears[1],"-01-15",sep="")),as.Date("2005-12-15"),by="month")
  outidx = which(as.numeric(substr(histmonthsused,1,4))>=1981)
  outputarray1 = outputarray1[,,outidx]
  outputarray2 = outputarray2[,,outidx]
  if(var=="pr"){
    outputarray3 = outputarray3[,,outidx]
    outputarray4 = outputarray4[,,outidx]
    outputarray5 = outputarray5[,,outidx]
  }
  
  var1climo = var2climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  if(var=="pr"){
    var3climo = var4climo = var5climo = array(NA,dim=c(length(lonidx),length(latidx),12))
  }
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    var1climo[,,month] = apply(outputarray1[,,monidx],c(1,2),mean,na.rm=TRUE)
    var2climo[,,month] = apply(outputarray2[,,monidx],c(1,2),mean,na.rm=TRUE)
    if(var=="pr"){
      var3climo[,,month] = apply(outputarray3[,,monidx],c(1,2),mean,na.rm=TRUE)
      var4climo[,,month] = apply(outputarray4[,,monidx],c(1,2),mean,na.rm=TRUE)
      var5climo[,,month] = apply(outputarray5[,,monidx],c(1,2),mean,na.rm=TRUE)
    }
    
    message("Finished climo calcs for month ",month)
  }
  
  
  dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
  dimY <- ncdim_def( "lat", "degrees_north", lat[latidx])
  dimT <- ncdim_def("month","month",1:12)
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  if(var=="tasmax"){
    unitsout1 = "degrees_C"
    unitsout2 = "days"
    varname1 = "tmaxclimo"
    varname2 = "tmax95climo"
  }
  if(var=="tasmin"){
    unitsout1 = "degrees_C"
    unitsout2 = "days"
    varname1 = "tminclimo"
    varname2 = "tmin32climo"
  }
  if(var=="pr"){
    unitsout1 = "mm"
    unitsout2 = "days"
    unitsout3 = "days"
    unitsout4 = "mm"
    unitsout5 = "mm"
    
    varname1 = "prclimo"
    varname2 = "pr50climo"
    varname3 = "r1mmclimo"
    varname4 = "rx1dayclimo"
    varname5 = "rx5dayclimo"
  }
  
  var1d <- ncvar_def(varname1,unitsout1, list(dimX,dimY,dimT), mv )
  var2d <- ncvar_def(varname2,unitsout2, list(dimX,dimY,dimT), mv )
  if(var=="pr"){
    var3d <- ncvar_def(varname3,unitsout3, list(dimX,dimY,dimT), mv )
    var4d <- ncvar_def(varname4,unitsout4, list(dimX,dimY,dimT), mv )
    var5d <- ncvar_def(varname5,unitsout5, list(dimX,dimY,dimT), mv )
  }
  
  #######
  # Create netcdf file
  if(var!="pr"){
    nc <- nc_create(paste("/home/woot0002/MACA_LIVNEH/",usefileinfo[1,1],"_",usefileinfo[1,2],"_",usefileinfo[1,3],"_",usefileinfo[1,4],"_histclimo.nc",sep="") ,  list(var1d,var2d) )
  } else {
    nc <- nc_create(paste("/home/woot0002/MACA_LIVNEH/",usefileinfo[1,1],"_",usefileinfo[1,2],"_",usefileinfo[1,3],"_",usefileinfo[1,4],"_histclimo.nc",sep="") ,  list(var1d,var2d,var3d,var4d,var5d) )
  }
  # Write some data to the file
  ncvar_put(nc, var1d, var1climo) # no start or count: write all values\
  ncvar_put(nc, var2d, var2climo) # no start or count: write all values\
  if(var=="pr"){
    ncvar_put(nc, var3d, var3climo) # no start or count: write all values\
    ncvar_put(nc, var4d, var4climo) # no start or count: write all values\
    ncvar_put(nc, var5d, var5climo) # no start or count: write all values\
  }
  
  # close ncdf
  nc_close(nc)
  message("Finished calcs for GCMs ",GCMs[i])
}
     
    
    

  


