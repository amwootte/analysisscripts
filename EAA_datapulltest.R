
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)
library(rgdal)
source("/data2/3to5/I35/scripts/analysisfunctions.R")

varname = "tasmax"
projfile = "/data2/3to5/I35/tasmax/PARM/tasmax_day_I35txdetrp1-PARM-B38D01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc"
histfile = "/data2/3to5/I35/tasmax/PARM/tasmax_day_I35txdetrp1-PARM-B30D01K00_historical_r1i1p1_I35Land_19810101-20051231.nc"
obsfile = "/data2/3to5/I35/tasmax/Daymet/tasmax_day_daymet_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

nctest = nc_open(obsfile)
times = ncvar_get(nctest,"time")
tmpV = ncvar_get(nctest,varname,start=c(1,1,1),count=c(-1,-1,1))
lon = ncvar_get(nctest,"lon")
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

#shapefile1 = "/home/woot0002/shapefiles/SanAntonioAirport"
#shapefile2 = "/home/woot0002/shapefiles/DurnellRanch"
#shapefile3 = "/home/woot0002/shapefiles/Jwells"

message("opening shapefile 1")
test1 = readOGR(dsn="/home/woot0002/shapefiles/",layer="SanAntonioAirport")
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
test_trans1 = spTransform(test1, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

message("opening shapefile 2")
test2 = readOGR(dsn="/home/woot0002/shapefiles/",layer="DurnellRanch")
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
test_trans2 = test2#spTransform(test2, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

message("opening shapefile 3")
test3 = readOGR(dsn="/home/woot0002/shapefiles/",layer="Jwells")
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
test_trans3 = test3 #spTransform(test3, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")


modrasV = raster(t(tmpV)[length(lat):1,])
#rm(tmpV)
if(all(lon>0)){
  extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
} else {
  extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
}

plot(modrasV,xlim=c(-100.05,-98),ylim=c(29,30),main=paste(varname," Jan 1, 1981",sep=""))
map("state",add=TRUE)
plot(test_trans1,add=TRUE)
plot(test_trans2,add=TRUE)
plot(test_trans3,add=TRUE)

####

coords1 = data.frame(coordinates(test_trans1))
coords2 = data.frame(coordinates(test_trans2))
coords3 = data.frame(coordinates(test_trans3))

names(coords1) = names(coords2) = names(coords3) = c("loc_lon","loc_lat")

coordsall = rbind(coords1,coords2)
coordsall = rbind(coordsall,coords3)

coordsall$location = c("SanAntonio1","SanAntonio2","SanAntonio3","DurnellRanch","Jwells1","Jwells2")

####

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

###
# get cells to use
R = c()
C = c()
modlon = c()
modlat = c()

for(n in 1:nrow(coordsall)){
  
  loc_lon = as.numeric(coordsall$loc_lon[n])
  loc_lat = as.numeric(coordsall$loc_lat[n])
  if(loc_lon>0) loc_lon=loc_lon-360
  
  pointarea = distfunc(loc_lon,loc_lat,modelgrid)
  R[n] = pointarea$R[1]
  C[n] = pointarea$C[1]
  modlon[n] = pointarea$lon[1]
  modlat[n] = pointarea$lat[1]
  
}

coordsall$R = R
coordsall$C = C
coordsall$modlon = modlon
coordsall$modlat = modlat


for(n in 1:nrow(coordsall)){
  
  ###
  # obs out
  
  nctest = nc_open(obsfile)
  times = ncvar_get(nctest,"time")
  tmpV = ncvar_get(nctest,varname,start=c(coordsall$R[n],coordsall$C[n],1),count=c(1,1,-1))
  if(varname=="pr") tmpV = tmpV*86400
  nc_close(nctest)
  
  histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  obsfileout = paste("/home/woot0002/EAA/",varname,"_1981-2005_Daymet_",coordsall$location[n],".csv",sep="")
  obsdataout = data.frame(histdates,tmpV)
  names(obsdataout) = c("dates",varname)
  write.table(obsdataout,obsfileout,sep=",",row.names=FALSE)
  
  ###
  # hist simulation out
  
  nctest = nc_open(histfile)
  times = ncvar_get(nctest,"time")
  tmpV = ncvar_get(nctest,varname,start=c(coordsall$R[n],coordsall$C[n],1),count=c(1,1,-1))
  if(varname=="pr") tmpV = tmpV*86400
  nc_close(nctest)
  
  histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  histfileout = paste("/home/woot0002/EAA/",varname,"_1981-2005_MPI-ESM-LR_PARM_Daymet_",coordsall$location[n],".csv",sep="")
  histdataout = data.frame(histdates,tmpV)
  names(histdataout) = c("dates",varname)
  write.table(histdataout,histfileout,sep=",",row.names=FALSE)
  
  ###
  # proj simulation out
  
  nctest = nc_open(projfile)
  times = ncvar_get(nctest,"time")
  tmpV = ncvar_get(nctest,varname,start=c(coordsall$R[n],coordsall$C[n],1),count=c(1,1,-1))
  if(varname=="pr") tmpV = tmpV*86400
  nc_close(nctest)
  
  projdates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
  projfileout = paste("/home/woot0002/EAA/",varname,"_2006-2099_RCP85_MPI-ESM-LR_PARM_Daymet_",coordsall$location[n],".csv",sep="")
  projdataout = data.frame(projdates,tmpV)
  names(projdataout) = c("dates",varname)
  write.table(projdataout,projfileout,sep=",",row.names=FALSE)
  
}

write.table(coordsall,"/home/woot0002/EAA/coordinates.csv",sep=",",row.names=FALSE)


