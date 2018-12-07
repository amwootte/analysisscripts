library(ncdf4)
library(maps)
library(fields)
library(sp)

varname = "tasmin"

futfiles = system(paste("ls /data2/3to5/I35/",varname,"/DeltaSD/*rcp85*.nc",sep=""),intern=TRUE)
files = system(paste("ls /data2/3to5/I35/",varname,"/EDQM/*historical*.nc",sep=""),intern=TRUE)

filesplit1 = do.call("rbind",strsplit(files,"/",fixed=TRUE))
filesplit2 = do.call("rbind",strsplit(filesplit1[,7],"-",fixed=TRUE))
filesplit3 = do.call("rbind",strsplit(filesplit2[,3],"_",fixed=TRUE))

filesplit = cbind(filesplit2[,1],filesplit3)
filesplit = cbind(filesplit,filesplit2[,4])

filesplit1 = do.call("rbind",strsplit(futfiles,"/",fixed=TRUE))
filesplit2 = do.call("rbind",strsplit(filesplit1[,7],"-",fixed=TRUE))
filesplit3 = do.call("rbind",strsplit(filesplit2[,3],"_",fixed=TRUE))

futfilesplit = cbind(filesplit2[,1],filesplit3)
futfilesplit = cbind(futfilesplit,filesplit2[,4])


###
# prep a mask file

nctest = nc_open(files[1])
mask = ncvar_get(nctest,varname,start=c(1,1,1),count=c(-1,-1,1))
nc_close(nctest)
mask = ifelse(is.na(mask)==TRUE,0,1)

######

for(i in 1:nrow(filesplit)){
  
  obsletter = substr(filesplit[i,2],4,4)
  if(obsletter=="D"){
    obsfile = system(paste("ls /home/woot0002/",varname,"*daymet*.nc",sep=""),intern=TRUE)
  }
  if(obsletter=="L"){
    obsfile = system(paste("ls /home/woot0002/",varname,"*livneh*.nc",sep=""),intern=TRUE)
  }
  if(obsletter=="P"){
    obsfile = system(paste("ls /home/woot0002/",varname,"*prism*.nc",sep=""),intern=TRUE)
  }
  
  nctest = nc_open(obsfile)
  lon = ncvar_get(nctest,"lon")
  lat = ncvar_get(nctest,"lat")
  times = ncvar_get(nctest,"time")
  
  latunits = nctest$dim[[1]]$units
  lonunits = nctest$dim[[3]]$units
  timeunits = nctest$dim[[4]]$units
  
  dataset = ncvar_get(nctest,varname)
  if(varname=="pr") dataunits = "kg m-2 s-1"
  if(varname=="tasmax" | varname=="tasmin") dataunits = "K"
  nc_close(nctest)
  
  mv = 1E20
  
  for(t in 1:length(times)){
    dataset[,,t] = ifelse(mask==1,dataset[,,t],NA)
  }
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits,lat)
  dimT <- ncdim_def("time",timeunits,times)
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varname,dataunits, list(dimX,dimY,dimT), mv )
  
  #######
  # Create netcdf file
  filename = paste("/home/woot0002/",futfilesplit[i,1],"-DeltaSD-",filesplit[i,2],"_",filesplit[i,3],"_",filesplit[i,4],"_",filesplit[i,5],"_",filesplit[i,6],"-",filesplit[i,7],sep="")
  nc <- nc_create(filename ,  list(var1d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, dataset) # no start or count: write all values\
   
  # close ncdf
  nc_close(nc)
  
  message("Finished historical write to match file: ",futfiles[i])
}





