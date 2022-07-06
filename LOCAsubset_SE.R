source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(RPushbullet)

var = "pr"
tempres = "monthly"
varout = "rx5day"

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",var,"/historical/*.nc",sep=""),intern=TRUE)

latbnds = c(14,53)
lonbnds = c(-125,-67)

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesin[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  varin = paste(filesplit2[1],filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  varout2 = paste(varout,filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
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
    timesin = which(as.numeric(substr(dates,1,4))>=1976)
    dates = dates[timesin]
    timestart = timesin[1]
    timecount = length(timesin)
    timeout= time[timesin]
    
    if(tempres=="monthly"){
      yearmon = unique(substr(dates,1,7))
      timeout = timeout[which(substr(dates,9,10)=="01")]
    } 
    
    LON = rep(lon,each=length(lat))
    LAT = rep(lat,length(lon))
    R = rep(1:length(lon),each=length(lat))
    C = rep(1:length(lat),length(lon))
    modelgrid = data.frame(R,C,LON,LAT)
    names(modelgrid) = c("R","C","lon","lat")
    if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
    
    idxuse = which(modelgrid$lat>=latbnds[1] & modelgrid$lat<=latbnds[2] & modelgrid$lon>=lonbnds[1] & modelgrid$lon<=lonbnds[2])
    modelgriduse = modelgrid[idxuse,]
    
    lonstart = min(modelgriduse$R)
    latstart = min(modelgriduse$C)
    loncount = length(unique(modelgriduse$R))
    latcount = length(unique(modelgriduse$C))
    
    lonout = unique(modelgriduse$lon)
    latout = unique(modelgriduse$lat)
    
  }
  
  #if(var=="pr"){
  #  testdat=testdat*86400
  #}
  
  
  if(tempres=="monthly"){
    tmp = array(NA,dim=c(loncount,latcount,length(yearmon)))
    for(y in 1:length(yearmon)){
      idxmon = which(substr(dates,1,7)==yearmon[y])
      
      testdat = ncvar_get(nctest,varin,start=c(lonstart,latstart,idxmon[1]),count=c(loncount,latcount,length(idxmon)))
      if(var=="tasmax"){
        if(testdat[320,144,1]>200){
          testdat = testdat-273.15
        }
      }
      
      if(mean(testdat,na.rm=TRUE)<1 & (var=="pr")){
        testdat = testdat*86400
      }
      
      if(varout=="pr"){
        tmp1 = apply(testdat,c(1,2),sum,na.rm=TRUE)  
      } 
      if(varout=="rx1day"){
        tmp1 = apply(testdat,c(1,2),max,na.rm=TRUE)
      }
      if(varout=="rx5day"){
        tmp1 = apply(testdat,c(1,2),calcrollsum,size=5)
      }
      if(varout=="tasmax" | varout=="tasmin"){
        tmp1 = apply(testdat,c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="tmax100"){
        td = ifelse(testdat >=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      
      tmp[,,y]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
      message("Finished with yearmon ",y," / ",length(yearmon))
    }
    nc_close(nctest)
    testdat=tmp
   
  }
  
  #####
  if(tempres=="daily"){
    filenameout = paste("/home/woot0002/RCMES/data/LOCA/",varout2,".nc",sep="")
  }
  if(tempres=="monthly"){
    filenameout = paste("/home/woot0002/RCMES/data/LOCA/",varout2,"_mon_CONUS.nc",sep="")
  }
  
  dimX <- ncdim_def( "lon", lonunits, lonout)
  dimY <- ncdim_def( "lat", latunits, latout)
  dimT <- ncdim_def("time",timeunits,timeout)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varout,dataunits,longname=varout, list(dimX,dimY,dimT), mv ,compression=9)

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
  
  #pbPost(type="note",title="LOCA trimmed file written out",body=paste("Trimmed file written for ",tempres," ",varin," it took ",ptmend[3]," secs.",sep=""),email="amwootte@ou.edu")
  
}

