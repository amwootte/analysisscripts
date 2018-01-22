#################################
#
# Downscaled Data Re-gridding Script version 2 - MACA
#
# Adrienne Wootten 12/08/2016
#
# This version accesses new functions built to address handling 
# missing values when interpolating
#
#################################
#
# Required Arguments
#
# args[1] - Dataset name
#    	Examples here include Hostetler, and Maurer
# args[2] - GCM name
#	Example here for GFDL
# args[3] - desired output resolution in km
#	using 50 as an example here
# args[4] - Parameter of interest
#	must be specific to the dataset in question.
#
#
# NOTE: File will auto regrid for all timeslices available
#
####################################

#args <- commandArgs(trailingOnly = TRUE)

args = c("maca","pr","rcp85")

dataset = args[1]
var = args[2]
period = args[3]

#files = system(paste("ls /home/users/amwootte/PR_downscaling/DSdata/rawoutput/",dataset,"/macav2livneh_",var,"*",period,"*.nc",sep=""),intern=TRUE)
files = system(paste("ls maca/macav2livneh_",var,"*",period,"*.nc",sep=""),intern=TRUE)
files
#####################

library(ncdf4)
library(sp)
library(fields)
source("regridfunctions.R")

# SERAP grid
newlat = c(33.5625, 33.6875, 33.8125, 33.9375, 34.0625, 34.1875, 34.3125, 34.4375, 34.5625, 34.6875, 34.8125, 34.9375, 35.0625, 35.1875, 35.3125, 35.4375, 35.5625, 35.6875, 35.8125, 35.9375, 36.0625, 36.1875, 36.3125,
    36.4375, 36.5625, 36.6875)

newlon = c( -79.9375, -79.8125, -79.6875, -79.5625, -79.4375, -79.3125, -79.1875, -79.0625, -78.9375, -78.8125, -78.6875, -78.5625, -78.4375, -78.3125,-78.1875, -78.0625, -77.9375, -77.8125, -77.6875, -77.5625, -77.4375,-77.3125, -77.1875, -77.0625, -76.9375, -76.8125, -76.6875, -76.5625,  -76.4375, -76.3125, -76.1875, -76.0625, -75.9375, -75.8125, -75.6875,-75.5625, -75.4375, -75.3125, -75.1875, -75.0625)

# 15km grid for the SE
#newlat=c(25.13511, 25.27025, 25.40538, 25.54052, 25.67565, 25.81079, 25.94592, 26.08106, 26.21619, 26.35133, 26.48646, 26.62160, 26.75673, 26.89187, 27.02700, 27.16214, 27.29727, 27.43241, 27.56754, 27.70268, 27.83781, 27.97295, 28.10808, 28.24322, 28.37835, 28.51349, 28.64862, 28.78376, 28.91889, 29.05403, 29.18916, 29.32430, 29.45943, 29.59457, 29.72970, 29.86484, 29.99997, 30.13511, 30.27024, 30.40538, 30.54051, 30.67565, 30.81078, 30.94592, 31.08105, 31.21619, 31.35132, 31.48646, 31.62159, 31.75673, 31.89186, 32.02700, 32.16213, 32.29727, 32.43240, 32.56754, 32.70267, 32.83781, 32.97294, 33.10808, 33.24321, 33.37835, 33.51348, 33.64862, 33.78375, 33.91889, 34.05402, 34.18916, 34.32429, 34.45943, 34.59456, 34.72970, 34.86483, 34.99997, 35.13510, 35.27024, 35.40537, 35.54051, 35.67564, 35.81078, 35.94591, 36.08105, 36.21618, 36.35132, 36.48645, 36.62159, 36.75672, 36.89186, 37.02699, 37.16213, 37.29726, 37.43240, 37.56753, 37.70267, 37.83780, 37.97294, 38.10807, 38.24321, 38.37834, 38.51348, 38.64861, 38.78375, 38.91888, 39.05402, 39.18915, 39.32429, 39.45942, 39.59456, 39.72969, 39.86483, 39.99996, 40.13510, 40.27023)

#newlon=c(-102.97288, -102.83774, -102.70261, -102.56747, -102.43234, -102.29720, -102.16207, -102.02693, -101.89180, -101.75666, -101.62153, -101.48639, -101.35126, -101.21612, -101.08099, -100.94585, -100.81072, -100.67558, -100.54045, -100.40531, -100.27018, -100.13504,  -99.99991,  -99.86477, -99.72964, -99.59450, -99.45937, -99.32423, -99.18910, -99.05396, -98.91883, -98.78369, -98.64856, -98.51342, -98.37829, -98.24315, -98.10802, -97.97288, -97.83775, -97.70261, -97.56748, -97.43234, -97.29721, -97.16207, -97.02694, -96.89180, -96.75667, -96.62153, -96.48640, -96.35126, -96.21613, -96.08099, -95.94586, -95.81072, -95.67559, -95.54045, -95.40532, -95.27018, -95.13505, -94.99991, -94.86478, -94.72964, -94.59451, -94.45937, -94.32424, -94.18910, -94.05397, -93.91883, -93.78370, -93.64856, -93.51343, -93.37829, -93.24316, -93.10802, -92.97289, -92.83775, -92.70262, -92.56748, -92.43235, -92.29721, -92.16208, -92.02694, -91.89181, -91.75667, -91.62154, -91.48640, -91.35127, -91.21613, -91.08100, -90.94586, -90.81073, -90.67559, -90.54046, -90.40532, -90.27019, -90.13505, -89.99992, -89.86478, -89.72965, -89.59451, -89.45938, -89.32424, -89.18911, -89.05397, -88.91884, -88.78370, -88.64857, -88.51343, -88.37830, -88.24316, -88.10803, -87.97289, -87.83776, -87.70262, -87.56749, -87.43235, -87.29722, -87.16208, -87.02695, -86.89181, -86.75668, -86.62154, -86.48641, -86.35127, -86.21614, -86.08100, -85.94587, -85.81073, -85.67560, -85.54046, -85.40533, -85.27019, -85.13506, -84.99992, -84.86479, -84.72965, -84.59452, -84.45938, -84.32425, -84.18911, -84.05398, -83.91884, -83.78371, -83.64857, -83.51344, -83.37830, -83.24317, -83.10803, -82.97290, -82.83776, -82.70263, -82.56749, -82.43236, -82.29722, -82.16209, -82.02695, -81.89182, -81.75668, -81.62155, -81.48641, -81.35128, -81.21614, -81.08101, -80.94587, -80.81074, -80.67560, -80.54047, -80.40533, -80.27020, -80.13506, -79.99993, -79.86479, -79.72966, -79.59452, -79.45939, -79.32425, -79.18912, -79.05398, -78.91885, -78.78371, -78.64858, -78.51344, -78.37831, -78.24317, -78.10804, -77.97290, -77.83777, -77.70263, -77.56750, -77.43236, -77.29723, -77.16209, -77.02696, -76.89182, -76.75669, -76.62155, -76.48642, -76.35128, -76.21615, -76.08101, -75.94588, -75.81074, -75.67561, -75.54047, -75.40534, -75.27020, -75.13507, -74.99993, -74.86480, -74.72966, -74.59453, -74.45939, -74.32426, -74.18912, -74.05399)
#newlat=c(33.51348, 33.64862, 33.78375, 33.91889, 34.05402, 34.18916, 34.32429, 34.45943, 34.59456, 34.72970, 34.86483, 34.99997, 35.13510, 35.27024, 35.40537, 35.54051, 35.67564, 35.81078, 35.94591, 36.08105, 36.21618, 36.35132, 36.48645, 36.62159)
#newlon = c(-79.99993, -79.86479, -79.72966, -79.59452, -79.45939, -79.32425, -79.18912, -79.05398, -78.91885, -78.78371, -78.64858, -78.51344, -78.37831, -78.24317, -78.10804, -77.97290, -77.83777, -77.70263, -77.56750, -77.43236, -77.29723, -77.16209, -77.02696, -76.89182, -76.75669, -76.62155, -76.48642, -76.35128, -76.21615, -76.08101, -75.94588, -75.81074, -75.67561, -75.54047, -75.40534, -75.27020, -75.13507)

newrownum = length(newlon) # number of cells in x direction
newcolnum = length(newlat)  # number of cells in y direction

for(f in 1:length(files)){
ptm1 = proc.time()

message("File regrid began at ",Sys.time()," for file: ",files[f])

# read file
test = nc_open(files[f])

#test = nc_open("/projdata/usda/pinemap/maca/past/test_hada.nc")

lonname = test$var[[1]]$dim[[1]]$name
latname = test$var[[1]]$dim[[2]]$name

modelname = test$var[[1]]$name
missingval = test$var[[1]]$missval

dataunits = test$var[[1]]$units

timelength = length(test$var[[1]]$dim[[3]]$vals)
times = ncvar_get(test,test$var[[1]]$dim[[3]]$name)
timeunits = test$var[[1]]$dim[[3]]$units

spunits = test$var[[1]]$dim[[2]]$units

newmatrix = array(data=NA,dim=c(newrownum,newcolnum,timelength))

##
# lat bound 1-252
# lon bound 342-815

lon = ncvar_get(test,lonname)#,start=342,count=473)-360
lat = ncvar_get(test,latname)#,start=1,count=252)

testvar = ncvar_get(test,modelname)#,start=c(342,1,1),count=c(473,252,-1))
nc_close(test)

#lat2 = seq(min(lat),max(lat),by=0.1)
#lon2 = seq(min(lon),max(lon),by=0.1)

lat2= lat
lon2 = lon

for(j in 1:timelength){

ptm = proc.time()

test2 = interp.surface.gridfix(list(x=lon2,y=lat2,z=testvar[,,j]),grid.list=list(x=newlon,y=newlat))

newvar = matrix(test2[[3]],nrow=length(newlon),ncol=length(newlat))
newvar2 = newvar

newmatrix[,,j] = newvar2

message("Regridded Time slice ",j," / ",timelength)
ptmend = proc.time()-ptm
message("Time = ",ptmend[3]," secs")

}

###########
# Write new netcdf for Regridded data

library(ncdf)

if(length(grep("degrees",spunits))>0){
dimX <- dim.def.ncdf( "lon", "degrees_north", newlon)
dimY <- dim.def.ncdf( "lat", "degrees_east", newlat )
} else {
dimX <- dim.def.ncdf( "x", spunits, newlon)
dimY <- dim.def.ncdf( "y", spunits, newlat )
}

dimT <- dim.def.ncdf( "time",timeunits,times)
# Make varables of various dimensionality, for illustration purposes
mv <- missingval # missing value to use

var1d <- var.def.ncdf(modelname, dataunits, list(dimX,dimY,dimT), mv )

# Create the test file
fileout = strsplit(files[f],"/",fixed=TRUE)[[1]][2]

nc <- create.ncdf(paste("commongrid/projection/",fileout,sep="") , var1d )

# Write some data to the file
put.var.ncdf( nc, var1d, newmatrix ) # no start or count: write all values\

# close ncdf
close.ncdf(nc)

message("Finished regridding file ",f," / ",length(files))
ptm1end = proc.time()-ptm1
message("Time = ",ptm1end[3]," secs")
message("File regrid complete at ",Sys.time())

}
