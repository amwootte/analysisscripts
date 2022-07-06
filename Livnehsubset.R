source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(RPushbullet)

var = "rx5day"
tempres = "monthly"
varfile = "Prec"

filesin = system(paste("ls /data4/data/OBS/LIVNEH/*.nc",sep=""),intern=TRUE)

filesplit = do.call(rbind,strsplit(filesin,split="/",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,6],split=".",fixed=TRUE))
filesplit2 = data.frame(filesplit2)
filesplit2$year = as.numeric(substr(filesplit2[,2],1,4))

years = c(1976,2005)
latbnds = c(14,53)
lonbnds = c(-125,-67)

fileidx = which(filesplit2$year>=years[1] & filesplit2$year<=years[2])

filesin = filesin[fileidx]
yearmon = as.character(unique(filesplit2[fileidx,2]))

timeout = c()

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  nctest = nc_open(filesin[i])
  
  if(i==1){
    varnames = names(nctest$var)
    dataunits = nctest$var[[which(varnames==varfile)]]$units
    
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    timeunits = nctest$var[[which(varnames==varfile)]]$dim[[3]]$units
    latunits = nctest$var[[which(varnames==varfile)]]$dim[[2]]$units
    lonunits = nctest$var[[which(varnames==varfile)]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    indexdate = as.Date(substr(timeunits,12,21))
    
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
    
    dataout = array(NA,dim=c(loncount,latcount,length(yearmon)))
    
  }
  
  time=ncvar_get(nctest,"time")
  dates=indexdate+time
  
  timeout = c(timeout,time[1])
  
  testdat = ncvar_get(nctest,varfile,start=c(lonstart,latstart,1),count=c(loncount,latcount,-1))
  nc_close(nctest)
  
  if(tempres=="monthly"){
      if(var=="pr"){
        tmp1 = apply(testdat,c(1,2),sum,na.rm=TRUE)  
      } 
     if(var=="rx1day"){
       tmp1 = apply(testdat,c(1,2),max,na.rm=TRUE)
      }
      if(var=="rx5day"){
        tmp1 = apply(testdat,c(1,2),calcrollsum,size=5)
      }
      if(var=="tasmax" | var=="tasmin"){
        tmp1 = apply(testdat,c(1,2),mean,na.rm=TRUE)
      }
      if(var=="tmax100"){
        td = ifelse(testdat>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      dataout[,,i]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
    }
  
}

if(tempres=="daily"){
  filenameout = paste("/home/woot0002/RCMES/data/OBS/",var,"_LIVNEH.nc",sep="")
}
if(tempres=="monthly"){
  filenameout = paste("/home/woot0002/RCMES/data/OBS/",var,"_LIVNEH_mon_CONUS.nc",sep="")
}

dimX <- ncdim_def( "lon", lonunits, lonout)
dimY <- ncdim_def( "lat", latunits, latout)
#timeout = timeout[-c(1:2)]
dimT <- ncdim_def("time",timeunits,timeout)

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

var1d <- ncvar_def(var,dataunits,longname=var, list(dimX,dimY,dimT), mv ,compression=9)

#######
# Create netcdf file

nc <- nc_create(filenameout ,  var1d )
# Write some data to the file
ncvar_put(nc, var1d, dataout) # no start or count: write all values\
# close ncdf
nc_close(nc)
rm(testdat)
gc()
ptmend = proc.time()-ptm

#pbPost(type="note",title="LIVNEH trimmed file written out",body=paste("Trimmed file written for ",tempres," ",var," it took ",ptmend[3]," secs.",sep=""),email="amwootte@ou.edu")
