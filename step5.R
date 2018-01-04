##################
#
# Step 5: File Format Change 

library(ncdf4)
library(rgdal)
library(raster)
library(rasterVis)
library(maps)
library(maptools)

# must pass in the following: step2_filename, varname, difftype, futureperiod,tempperiod,outfiletype
# configured currently only for annual tempperiod and GeoTIFFs

test = nc_open(step2_filename)
projdiff_rcp26= ncvar_get(test,"projmeandiff_rcp26")
projdiff_rcp45= ncvar_get(test,"projmeandiff_rcp45")
projdiff_rcp85= ncvar_get(test,"projmeandiff_rcp85")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

if(outfiletype="GTiff"){

dataras = raster(t(projdiff_rcp26[,length(lat):1]))
extent(dataras) = c(lon[1]-360,lon[length(lon)]-360,lat[1],lat[length(lat)])
crs(dataras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- writeRaster(dataras, filename=paste("/data2/3to5/I35/GeoTIFFs/",var,"_rcp26_meanchange_",difftype,"_",futureperiod[1],"-",futureperiod[2],".tif",sep=""), format="GTiff", overwrite=TRUE)

dataras = raster(t(projdiff_rcp45[,length(lat):1]))
extent(dataras) = c(lon[1]-360,lon[length(lon)]-360,lat[1],lat[length(lat)])
crs(dataras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- writeRaster(dataras, filename=paste("/data2/3to5/I35/GeoTIFFs/",var,"_rcp45_meanchange_",difftype,"_",futureperiod[1],"-",futureperiod[2],".tif",sep=""), format="GTiff", overwrite=TRUE)

dataras = raster(t(projdiff_rcp85[,length(lat):1]))
extent(dataras) = c(lon[1]-360,lon[length(lon)]-360,lat[1],lat[length(lat)])
crs(dataras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- writeRaster(dataras, filename=paste("/data2/3to5/I35/GeoTIFFs/",var,"_rcp85_meanchange_",difftype,"_",futureperiod[1],"-",futureperiod[2],".tif",sep=""), format="GTiff", overwrite=TRUE)

}