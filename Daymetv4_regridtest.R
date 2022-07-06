source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rgdal)
#library(RPushbullet)

var = "prcp"
tempres = "monthly"
varfile = "prcp"

filesin = system(paste("ls /data4/data/OBS/DAYMETv4/na/pr/*.nc",sep=""),intern=TRUE)

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

#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

i=1

nctest = nc_open(filesin[i])
varnames = names(nctest$var)
times = ncvar_get(nctest,"time")
timeunits = nctest$var[[which(varnames==varfile)]]$dim[[3]]$units
#timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
indexdate = as.Date(substr(timeunits,12,21))
nc_close(nctest)

dates=indexdate+times


t=1

pft <- raster(filesin[i],band=t,varname="prcp")
test = projectRaster(pft,crs=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
system("rm /tmp/Rtmpcmd34W/raster/*")
#dates=seq(as.Date("1980-01-01"),as.Date("1980-12-31"),by="day")
if(t==1){
  bounds = extent(lonbnds,latbnds)
}

testcropped = crop(test,bounds)

testcroppedmat = raster::as.matrix(testcropped)

if(t==1){
  coords = coordinates(testcropped)
  lon = unique(coords[,1])
  lat = unique(coords[,2])
  lat = rev(lat)
}

testcroppedmat2 = t(testcroppedmat)[,length(lat):1]

if(t==1){
  
}

testsfc = list(x=lon,y=lat,z=testcroppedmat2)
surface(testsfc,type="I")
map("world",add=TRUE)
map("state",add=TRUE)
