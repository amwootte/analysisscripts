source("/data2/3to5/I35/scripts/analysisfunctions.R")
source("/data2/3to5/I35/scripts/springpheno_v0.5.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(multiApply)

filesin = system("ls /data4/data/DS_proj/LOCA/pr/historical/*.nc",intern=TRUE)

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
    timeout= time
    
    yearmon = unique(substr(dates,1,7))
    timeoutmon = timeout[which(substr(dates,9,10)=="01")]
    
  }
  
  tmp = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
  for(y in 1:length(yearmon)){
    idxmon = which(substr(dates,1,7)==yearmon[y])
    testdat = ncvar_get(nctest,varin,start=c(1,1,idxmon[1]),count=c(-1,-1,length(idxmon)))
    if(mean(testdat,na.rm=TRUE)<1){
      testdat = testdat*86400
    }
    
      tmp1 = apply(testdat,c(1,2),sum,na.rm=TRUE)  
    tmp[,,y]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
    rm(testdat)
    gc()
    message("Finished with yearmon ",y," / ",length(yearmon))
  }
  nc_close(nctest)
  testdat=tmp
  
  filenamepart = paste(filesplit2[1],"mon",filesplit2[2],filesplit2[3],filesplit2[4],sep="_")
  filenameout = paste("/data4/data/DS_proj/LOCA/prmon/historical/",filenamepart,sep="")
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits, lat)
  dimT <- ncdim_def("time",timeunits,timeoutmon)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varin,dataunits,longname=varin, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, testdat) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(testdat)
  gc()
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(filesin))
    
  }
  
 
 

