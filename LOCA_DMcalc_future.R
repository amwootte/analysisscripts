library(iterators)
#library(parallel)
#library(foreach)
#library(doParallel)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)
library(SPEI)

DM = 0

if(DM==0){
  boundsin = c(-0.75,-0.5)
}
if(DM==1){
  boundsin = c(-1.25,-0.75)
}
if(DM==2){
  boundsin = c(-1.55,-1.25)
}
if(DM==3){
  boundsin = c(-1.95,-1.55)
}
if(DM==4){
  boundsin = c(-1000,-1.95)
}

filesin = system("ls /data4/data/DS_proj/LOCA/SPI1/future/*.nc",intern=TRUE)

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesin[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  varin = paste(filesplit2[1],filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
  nctest = nc_open(filesin[i])
  
  if(i==1){
    dataunits = nctest$var[[1]]$units
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    time=ncvar_get(nctest,"time")
    timeunits = nctest$var[[1]]$dim[[3]]$units
    latunits = nctest$var[[1]]$dim[[2]]$units
    lonunits = nctest$var[[1]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    timeoutmon= time
    
    yearmon = unique(substr(dates,1,4))
  }
  
  testdat = ncvar_get(nctest,varin)
  
  start_time <- Sys.time()
  DMout = ifelse(testdat>=boundsin[1] & testdat<boundsin[2],1,0)
  DMout = ifelse(is.na(testdat)==FALSE,DMout,NA)
  end_time <- Sys.time()
  end_time - start_time
  
  nc_close(nctest)
  
  varout=paste("DM",DM,sep="")
  filenamepart = paste(varout,filesplit2[2],filesplit2[3],filesplit2[4],sep="_")
  filenameout = paste("/data4/data/DS_proj/LOCA/DM",DM,"/future/",filenamepart,sep="")
  
  varout2 = paste(varout,filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits, lat)
  dimT <- ncdim_def("time",timeunits,timeoutmon)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varout2,dataunits,longname=varout2, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, DMout) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(DMout)
  gc()
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(filesin)," Time:",ptmend[3]," secs")
    
}

