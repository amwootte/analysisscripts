source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rgdal)
#library(RPushbullet)

var = "tasmax"
varfile = "tmax"

unitsout="C"

filesin = system(paste("ls /data4/data/OBS/DAYMETv4/na/",var,"/*.nc",sep=""),intern=TRUE)

filesplit = do.call(rbind,strsplit(filesin,split="/",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,ncol(filesplit)],split="_",fixed=TRUE))
filesplit2 = data.frame(filesplit2)
filesplit2$year = as.numeric(substr(filesplit2[,ncol(filesplit2)],1,4))

years = c(1980,2020)
latbnds = c(14,53)
lonbnds = c(-125,-67)

fileidx = which(filesplit2$year>=years[1] & filesplit2$year<=years[2])

filesin = filesin[fileidx]
yearmon = as.character(unique(filesplit2[fileidx,2]))

for(i in 7:length(filesin)){

  filesplit = strsplit(filesin[i],"/",fixed=TRUE)
  filenamein = filesplit[length(filesplit)]
  
nctest = nc_open(filesin[i])
varnames = names(nctest$var)
times = ncvar_get(nctest,"time")
timeunits = nctest$var[[which(varnames==varfile)]]$dim[[3]]$units
#timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
indexdate = as.Date(substr(timeunits,12,21))
nc_close(nctest)

for(t in 1:length(times)){

  ptm=proc.time()
  
pft <- raster(filesin[i],band=t,varname=varfile)

if(t==1){
bnds1 = c(-2000000,2500000)
bnds2 = c(-1900000,1000000)
bnds = extent(bnds1,bnds2)
bounds = extent(lonbnds,latbnds)
}

pftcropped = crop(pft,bnds)
test = projectRaster(pftcropped,crs=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

rm(pft)
rm(pftcropped)
#testcropped = crop(test,bounds)

testcroppedmat = raster::as.matrix(test)

if(t==1){
  coords = coordinates(test)
  lon = unique(coords[,1])
  lat = unique(coords[,2])
  lat = rev(lat)
  
  filenameout1 = paste("/home/woot0002/RCMES/data/OBS/DAYMETv4/",filenamein[[1]][length(filenamein[[1]])],sep="")
  
  dimX <- ncdim_def( "lon", "degrees W", lon)
  dimY <- ncdim_def( "lat", "degrees N", lat)
  dimT <- ncdim_def("time",timeunits,times)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varfile,unitsout,longname=varfile, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout1 ,  var1d )
} else {
  nc = nc_open(filenameout1,write=TRUE)
}

outarray = t(testcroppedmat)[,length(lat):1]

#sourcename = test@file@name
#sourcesplit = strsplit(sourcename,"/",fixed=TRUE)[[1]]
#system(paste("rm /",sourcesplit[2],"/",sourcesplit[3],"/",sourcesplit[4],"/*",sep=""))


# Write some data to the file
ncvar_put(nc, var1d, outarray,start=c(1,1,t),count=c(-1,-1,1)) # no start or count: write all values\
# close ncdf
nc_close(nc)

ptmend = proc.time()-ptm
message("Finished reprojection and regrid for time ",t," / ",length(times))
message("time ",ptmend[3]," secs ")


}
message("Finished reprojection and regrid for file ",i," / ",length(filesin))
}
