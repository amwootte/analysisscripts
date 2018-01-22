#################################
#
# Re-gridding Script
#
# Adrienne Wootten 8/19/2015
#
#################################
#
# NOTE: File will auto regrid for all timeslices available
#
####################################

#/data/osu/prism/daily/combo

#args <- commandArgs(trailingOnly = TRUE)

#args = c("wicci","prcp","twentyc3m")

#dataset = args[1]
#var = args[2]
#period = args[3]

#files = system(paste("ls /data/osu/prism/daily/combo/",period,"*",var,"*.nc",sep=""),intern=TRUE)
files = c()
for(y in 1980:2014){
files = c(files,system(paste("ls NLDAS2/",y,"/*",sep=""),intern=TRUE))
}


#files
#####################

library(ncdf)
library(sp)
library(fields)

newlat=c(25.13511, 25.27025, 25.40538, 25.54052, 25.67565, 25.81079, 25.94592, 26.08106, 26.21619, 26.35133, 26.48646, 26.62160, 26.75673, 26.89187, 27.02700, 27.16214, 27.29727, 27.43241, 27.56754, 27.70268, 27.83781, 27.97295, 28.10808, 28.24322, 28.37835, 28.51349, 28.64862, 28.78376, 28.91889, 29.05403, 29.18916, 29.32430, 29.45943, 29.59457, 29.72970, 29.86484, 29.99997, 30.13511, 30.27024, 30.40538, 30.54051, 30.67565, 30.81078, 30.94592, 31.08105, 31.21619, 31.35132, 31.48646, 31.62159, 31.75673, 31.89186, 32.02700, 32.16213, 32.29727, 32.43240, 32.56754, 32.70267, 32.83781, 32.97294, 33.10808, 33.24321, 33.37835, 33.51348, 33.64862, 33.78375, 33.91889, 34.05402, 34.18916, 34.32429, 34.45943, 34.59456, 34.72970, 34.86483, 34.99997, 35.13510, 35.27024, 35.40537, 35.54051, 35.67564, 35.81078, 35.94591, 36.08105, 36.21618, 36.35132, 36.48645, 36.62159, 36.75672, 36.89186, 37.02699, 37.16213, 37.29726, 37.43240, 37.56753, 37.70267, 37.83780, 37.97294, 38.10807, 38.24321, 38.37834, 38.51348, 38.64861, 38.78375, 38.91888, 39.05402, 39.18915, 39.32429, 39.45942, 39.59456, 39.72969, 39.86483, 39.99996, 40.13510, 40.27023)

newlon=c(-102.97288, -102.83774, -102.70261, -102.56747, -102.43234, -102.29720, -102.16207, -102.02693, -101.89180, -101.75666, -101.62153, -101.48639, -101.35126, -101.21612, -101.08099, -100.94585, -100.81072, -100.67558, -100.54045, -100.40531, -100.27018, -100.13504,  -99.99991,  -99.86477, -99.72964, -99.59450, -99.45937, -99.32423, -99.18910, -99.05396, -98.91883, -98.78369, -98.64856, -98.51342, -98.37829, -98.24315, -98.10802, -97.97288, -97.83775, -97.70261, -97.56748, -97.43234, -97.29721, -97.16207, -97.02694, -96.89180, -96.75667, -96.62153, -96.48640, -96.35126, -96.21613, -96.08099, -95.94586, -95.81072, -95.67559, -95.54045, -95.40532, -95.27018, -95.13505, -94.99991, -94.86478, -94.72964, -94.59451, -94.45937, -94.32424, -94.18910, -94.05397, -93.91883, -93.78370, -93.64856, -93.51343, -93.37829, -93.24316, -93.10802, -92.97289, -92.83775, -92.70262, -92.56748, -92.43235, -92.29721, -92.16208, -92.02694, -91.89181, -91.75667, -91.62154, -91.48640, -91.35127, -91.21613, -91.08100, -90.94586, -90.81073, -90.67559, -90.54046, -90.40532, -90.27019, -90.13505, -89.99992, -89.86478, -89.72965, -89.59451, -89.45938, -89.32424, -89.18911, -89.05397, -88.91884, -88.78370, -88.64857, -88.51343, -88.37830, -88.24316, -88.10803, -87.97289, -87.83776, -87.70262, -87.56749, -87.43235, -87.29722, -87.16208, -87.02695, -86.89181, -86.75668, -86.62154, -86.48641, -86.35127, -86.21614, -86.08100, -85.94587, -85.81073, -85.67560, -85.54046, -85.40533, -85.27019, -85.13506, -84.99992, -84.86479, -84.72965, -84.59452, -84.45938, -84.32425, -84.18911, -84.05398, -83.91884, -83.78371, -83.64857, -83.51344, -83.37830, -83.24317, -83.10803, -82.97290, -82.83776, -82.70263, -82.56749, -82.43236, -82.29722, -82.16209, -82.02695, -81.89182, -81.75668, -81.62155, -81.48641, -81.35128, -81.21614, -81.08101, -80.94587, -80.81074, -80.67560, -80.54047, -80.40533, -80.27020, -80.13506, -79.99993, -79.86479, -79.72966, -79.59452, -79.45939, -79.32425, -79.18912, -79.05398, -78.91885, -78.78371, -78.64858, -78.51344, -78.37831, -78.24317, -78.10804, -77.97290, -77.83777, -77.70263, -77.56750, -77.43236, -77.29723, -77.16209, -77.02696, -76.89182, -76.75669, -76.62155, -76.48642, -76.35128, -76.21615, -76.08101, -75.94588, -75.81074, -75.67561, -75.54047, -75.40534, -75.27020, -75.13507, -74.99993, -74.86480, -74.72966, -74.59453, -74.45939, -74.32426, -74.18912, -74.05399)

#newlat=c(33.51348, 33.64862, 33.78375, 33.91889, 34.05402, 34.18916, 34.32429, 34.45943, 34.59456, 34.72970, 34.86483, 34.99997, 35.13510, 35.27024, 35.40537, 35.54051, 35.67564, 35.81078, 35.94591, 36.08105, 36.21618, 36.35132, 36.48645, 36.62159)
#newlon = c(-79.99993, -79.86479, -79.72966, -79.59452, -79.45939, -79.32425, -79.18912, -79.05398, -78.91885, -78.78371, -78.64858, -78.51344, -78.37831, -78.24317, -78.10804, -77.97290, -77.83777, -77.70263, -77.56750, -77.43236, -77.29723, -77.16209, -77.02696, -76.89182, -76.75669, -76.62155, -76.48642, -76.35128, -76.21615, -76.08101, -75.94588, -75.81074, -75.67561, -75.54047, -75.40534, -75.27020, -75.13507)

newrownum = length(newlon) # number of cells in x direction
newcolnum = length(newlat)  # number of cells in y direction

times = 1:length(files)-1

# Go through all the files to do regridding for each time

for(f in 1:length(files)){
ptm1 = proc.time()

message("File regrid began at ",Sys.time()," for file: ",files[f])

# read file
test = open.ncdf(files[f])

lonname = test$var[[1]]$dim[[1]]$name
latname = test$var[[1]]$dim[[2]]$name

lon = get.var.ncdf(test,lonname)
lat = get.var.ncdf(test,latname)

sravgname = test$var[[1]]$name
srtotname = test$var[[10]]$name

presname = test$var[[2]]$name
spechumname = test$var[[3]]$name

tavgname = test$var[[4]]$name
tmaxname = test$var[[8]]$name
tminname = test$var[[7]]$name

ugrdname = test$var[[5]]$name
vgrdname = test$var[[6]]$name

prcpname = test$var[[9]]$name

#timelength = length(test$var[[1]]$dim[[3]]$vals)

times = get.var.ncdf(test,test$var[[1]]$dim[[3]]$name)
timeunits = test$var[[1]]$dim[[3]]$units
spunits = test$var[[1]]$dim[[2]]$units
missingval = test$var[[1]]$missval

tempunits = test$var[[4]]$units
prcpunits = test$var[[9]]$units
windunits = test$var[[5]]$units
srunits = test$var[[1]]$units
humunits = test$var[[3]]$units
presunits = test$var[[2]]$units

tavg = get.var.ncdf(test,tavgname)
tmax = get.var.ncdf(test,tmaxname)
tmin = get.var.ncdf(test,tminname)

prcp = get.var.ncdf(test,prcpname)

ugrd = get.var.ncdf(test,ugrdname)
vgrd = get.var.ncdf(test,vgrdname)

pres = get.var.ncdf(test,presname)
spechum = get.var.ncdf(test,spechumname)

sravg = get.var.ncdf(test,sravgname)
srtot = get.var.ncdf(test,srtotname)

close.ncdf(test)

#####
# Regrid everything

#ptm = proc.time()

###
# temperatures K

test2tavg = interp.surface.grid(list(x=lon,y=lat,z=tavg),grid.list=list(x=newlon,y=newlat))
newtavg = matrix(test2tavg[[3]],nrow=length(newlon),ncol=length(newlat))
test2tmax = interp.surface.grid(list(x=lon,y=lat,z=tmax),grid.list=list(x=newlon,y=newlat))
newtmax = matrix(test2tmax[[3]],nrow=length(newlon),ncol=length(newlat))
test2tmin = interp.surface.grid(list(x=lon,y=lat,z=tmin),grid.list=list(x=newlon,y=newlat))
newtmin = matrix(test2tmin[[3]],nrow=length(newlon),ncol=length(newlat))

####
# precip mm

test2prcp = interp.surface.grid(list(x=lon,y=lat,z=prcp),grid.list=list(x=newlon,y=newlat))
newprcp = matrix(test2prcp[[3]],nrow=length(newlon),ncol=length(newlat))

####
# sr W/m^2

test2sravg = interp.surface.grid(list(x=lon,y=lat,z=sravg),grid.list=list(x=newlon,y=newlat))
newsravg = matrix(test2sravg[[3]],nrow=length(newlon),ncol=length(newlat))
test2srtot = interp.surface.grid(list(x=lon,y=lat,z=srtot),grid.list=list(x=newlon,y=newlat))
newsrtot = matrix(test2srtot[[3]],nrow=length(newlon),ncol=length(newlat))

####
# pressure Pa

test2pres = interp.surface.grid(list(x=lon,y=lat,z=pres),grid.list=list(x=newlon,y=newlat))
newpres = matrix(test2pres[[3]],nrow=length(newlon),ncol=length(newlat))

####
# specific humidity

test2shum = interp.surface.grid(list(x=lon,y=lat,z=spechum),grid.list=list(x=newlon,y=newlat))
newshum = matrix(test2shum[[3]],nrow=length(newlon),ncol=length(newlat))

####
# winds m/s

test2ugrd = interp.surface.grid(list(x=lon,y=lat,z=ugrd),grid.list=list(x=newlon,y=newlat))
newugrd = matrix(test2ugrd[[3]],nrow=length(newlon),ncol=length(newlat))
test2vgrd = interp.surface.grid(list(x=lon,y=lat,z=vgrd),grid.list=list(x=newlon,y=newlat))
newvgrd = matrix(test2vgrd[[3]],nrow=length(newlon),ncol=length(newlat))


#message("Regridded Time slice ",j," / ",timelength)
#ptmend = proc.time()-ptm
#message("Time = ",ptmend[3]," secs")

###########
# Write new netcdf for Regridded data

dimX <- dim.def.ncdf( "lon", "degrees_east", newlon)
dimY <- dim.def.ncdf( "lat", "degrees_north", newlat )

dimT <- dim.def.ncdf( "time",timeunits,times)
# Make varables of various dimensionality, for illustration purposes
mv <- missingval # missing value to use

var1d <- var.def.ncdf(tmaxname, "K", list(dimX,dimY,dimT), mv )
var2d <- var.def.ncdf(tminname, "K", list(dimX,dimY,dimT), mv )
var3d <- var.def.ncdf(tavgname, "K", list(dimX,dimY,dimT), mv )

var4d <- var.def.ncdf(prcpname, "mm", list(dimX,dimY,dimT), mv )

var5d <- var.def.ncdf(ugrdname, "m/s", list(dimX,dimY,dimT), mv )
var6d <- var.def.ncdf(vgrdname, "m/s", list(dimX,dimY,dimT), mv )

var7d <- var.def.ncdf(sravgname, "W/m^2", list(dimX,dimY,dimT), mv )
var8d <- var.def.ncdf(srtotname, "W/m^2", list(dimX,dimY,dimT), mv )

var9d <- var.def.ncdf(presname, "Pa", list(dimX,dimY,dimT), mv )
var10d <- var.def.ncdf(spechumname, "", list(dimX,dimY,dimT), mv )

filesplit = strsplit(files[f],"/")[[1]]
fileout = paste(filesplit[1],"/res15km/SE/",filesplit[2],"/",filesplit[3],sep="")

# Create the file
nc <- create.ncdf(fileout , list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d) )

# Write some data to the file
put.var.ncdf( nc, var1d, newtmax ) # no start or count: write all values\
put.var.ncdf( nc, var2d, newtmin )
put.var.ncdf( nc, var3d, newtavg )

put.var.ncdf( nc, var4d, newprcp )

put.var.ncdf( nc, var5d, newugrd )
put.var.ncdf( nc, var6d, newvgrd )

put.var.ncdf( nc, var7d, newsravg )
put.var.ncdf( nc, var8d, newsrtot )

put.var.ncdf( nc, var9d, newpres )
put.var.ncdf( nc, var10d, newshum )


# close ncdf
close.ncdf(nc)

message("Finished regridding file ",f," / ",length(files))
ptm1end = proc.time()-ptm1
message("Time = ",ptm1end[3]," secs")
message("File regrid complete at ",Sys.time())

}



