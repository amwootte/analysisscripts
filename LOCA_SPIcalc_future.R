source("/data2/3to5/I35/scripts/analysisfunctions.R")
source("/data2/3to5/I35/scripts/springpheno_v0.5.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(multiApply)
require(SPEI)

tempres=1
filesin = system("ls /data4/data/DS_proj/LOCA/prmon/future/*.nc",intern=TRUE)

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesin[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  varin = paste(filesplit2[1],filesplit2[3],filesplit2[4],substr(filesplit2[5],1,nchar(filesplit2[5])-3),sep="_")
  
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
    dates=seq(as.Date("2006-01-01"),as.Date("2100-12-01"),by="month")
    timeoutmon= time
    
    yearmon = unique(substr(dates,1,7))
  }
  
  testdat = ncvar_get(nctest,varin)
  
  start_time <- Sys.time()
  SPIout = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
  #SPIout[r,c,]= 
  for(r in 1:length(lon)){
    for(c in 1:length(lat)){
      if(all(is.na(testdat[r,c,])==TRUE)==FALSE){
        SPIout[r,c,] = spi(testdat[r,c,],scale=tempres,distribution = 'Gamma',na.rm=TRUE)$fitted
        SPIout[r,c,] = ifelse(SPIout[r,c,]== -Inf | SPIout[r,c,]== Inf,NA,SPIout[r,c,])
      }
    }
  }
  end_time <- Sys.time()
  end_time - start_time
  
  
  nc_close(nctest)
  
  varout=paste("SPI",tempres,sep="")
  filenamepart = paste(varout,filesplit2[3],filesplit2[4],filesplit2[5],sep="_")
  filenameout = paste("/data4/data/DS_proj/LOCA/SPI",tempres,"/future/",filenamepart,sep="")
  
  varout2 = paste(varout,filesplit2[3],filesplit2[4],substr(filesplit2[5],1,nchar(filesplit2[5])-3),sep="_")
  
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
  ncvar_put(nc, var1d, SPIout) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(SPIout)
  gc()
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(filesin)," Time:",ptmend[3]," secs")
    
  }
  
 
 

