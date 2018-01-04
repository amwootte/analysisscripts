##################
#
# GeoTIFF write test with 3^5 output

library(ncdf4)
library(rgdal)
library(raster)
library(rasterVis)
library(maps)
library(maptools)

var = "heatwaves"
scen = "rcp26"

test = nc_open(paste("/data2/3to5/I35/ens_means/",var,"_ensmean_absolute_2041-2070.nc",sep=""))
vardata= ncvar_get(test,paste("projmeandiff_",scen,sep=""))
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

dataras = raster(t(vardata[,length(lat):1]))
extent(dataras) = c(lon[1]-360,lon[length(lon)]-360,lat[1],lat[length(lat)])
plot(dataras)
map("state",add=TRUE)

crs(dataras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

rf <- writeRaster(dataras, filename=paste(var,"_",scen,"_meanchange_2041-2070.tif",sep=""), format="GTiff", overwrite=TRUE)

