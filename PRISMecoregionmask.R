
library(rgdal)
library(raster)
library(ncdf4)
library(maps)
library(fields)
library(mapdata)
library(rasterVis)
library(maptools)

setwd("shapefiles")
mapvals = readShapePoly("caribbean_usgs_country_bound_geo84")

test = readShapePoly("PR_lifezones_2007coast_10jan2012")
projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

##################
# WRF basics

test2 = nc_open("pr_ppt.01.nc")
lonname = test2$var[[1]]$dim[[1]]$name
latname = test2$var[[1]]$dim[[2]]$name
lonwrf = ncvar_get(test2,lonname)
latwrf = ncvar_get(test2,latname)
spunits = test2$var[[1]]$dim[[2]]$units
WRFdata = ncvar_get(test2,"pr_ppt")
nc_close(test2)

lonidx = which(lonwrf<= -65.585 & lonwrf>= -67.27)
latidx = which(latwrf<= 18.53 & latwrf>= 17.925)

#plot(WRF.sub1)
#plot(test,add=TRUE)

#par(mfrow=c(1,1))

#plot(r.sub)
#plot(test, add = TRUE)

testsfc = list(x=lonwrf,y=rev(latwrf),z=WRFdata[,length(latwrf):1])
surface(testsfc,type="I",xlim=c(-67.27,-65.585),ylim=c(17.925,18.53))
plot(mapvals,add=TRUE)

#WRFras = raster(t(WRFdata[lonidx,latidx])[length(latidx):1,])
WRFras = raster(t(WRFdata[lonidx,latidx]))
extent(WRFras) = c(-67.27,-65.585,17.925,18.53)
plot(WRFras)
plot(mapvals,add=TRUE)

WRFras2 = WRFras
extent(WRFras2)=extent(test)

WRFcoordst = coordinates(WRFras2)
xtrans = unique(WRFcoordst[,1])
ytrans = unique(WRFcoordst[,2])

WRFcoordst = data.frame(WRFcoordst)
WRFcoordst$R = rep(lonidx,length(latidx))
WRFcoordst$C = rep(latidx,each=length(lonidx))
WRFcoordst$NAME = paste(WRFcoordst$x,WRFcoordst$y,sep="_")

WRFeco = matrix(NA,nrow=length(lonidx),ncol=length(latidx))

for(i in 1:6){
  
  # Subset to Eco-region
  nestates <- c(i)
  test.sub <- test[as.character(test@data$code) %in% nestates, ]
  
  # Crop elevation data by extent of state subset
  WRF.sub <- crop(WRFras2, extent(test.sub))
  WRF.sub <- mask(WRFras2, test.sub)
  
  vals = getValues(WRF.sub)
  valsin = ifelse(is.na(vals)==FALSE,i,NA)
  
  WRFecotemp = matrix(valsin,nrow=length(lonidx),ncol=length(latidx))
  WRFecotemp = WRFecotemp[,length(latidx):1]
  WRFeco[which(is.na(WRFecotemp)==FALSE,arr.ind=TRUE)]=WRFecotemp[which(is.na(WRFecotemp)==FALSE,arr.ind=TRUE)]
  
}

testsfc = list(x=lonwrf[lonidx],y=rev(latwrf[latidx]),z=WRFeco)
surface(testsfc,type="I")
plot(mapvals,add=TRUE)

WRFecoout = matrix(NA,nrow=length(lonwrf),ncol=length(latwrf))
WRFecoout[lonidx,rev(latidx)]=WRFeco

WRFras = raster(t(WRFecoout[lonidx,latidx]))
extent(WRFras) = c(-67.27,-65.585,17.925,18.53)
plot(WRFras)
plot(mapvals,add=TRUE)

testsfc = list(x=lonwrf,y=rev(latwrf),z=WRFecoout[,length(latwrf):1])
surface(testsfc,type="I")
plot(mapvals,add=TRUE)

#latwrf = rev(latwrf)

###
# make mask into netcdf format

dimX <- ncdim_def( "lon", "degrees_east", lonwrf)
dimY <- ncdim_def( "lat", "degrees_north", latwrf )

# Make varables of various dimensionality, for illustration purposes
mv <- -999 # missing value to use

var1d <- ncvar_def("code", "",longname="Eco-region Codes", list(dimX,dimY), mv )

# Create the test file
nc <- nc_create("PRISMecoregionmask.nc" , var1d )

# Write some data to the file
ncvar_put( nc, var1d, WRFecoout ) # no start or count: write all values\

# close ncdf
nc_close(nc)


