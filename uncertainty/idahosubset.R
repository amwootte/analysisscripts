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
#for(y in 1980:2014){
#/projdata/usda/pinemap/idaho
#files = c(files,system(paste("ls NLDAS2/",y,"/*",sep=""),intern=TRUE))
#}

variable ="vs"

files = system(paste("ls /projdata/usda/pinemap/idaho_v2/",variable,"*.nc",sep=""),intern=TRUE)

#files
#####################

library(ncdf4)
library(sp)
library(fields)

latrange = c(33.51348,36.62159)
lonrange = c(-80,-75)

#times = 1:length(files)-1
# Go through all the files to do regridding for each time

years = do.call(rbind,strsplit(files,"_",fixed=TRUE))
years = as.numeric(substr(years[,3],1,4))

times = seq(as.Date(paste(years[1],"-01-01",sep="")),as.Date(paste(years[length(years)],"-12-31",sep="")),by="day")

timein = c()

newmatrix = array(NA,dim=c(120,75,length(times)))

for(f in 1:length(files)){
ptm1 = proc.time()

message("File subset began at ",Sys.time()," for file: ",files[f])

# read file
test = nc_open(files[f])

lonname = test$var[[1]]$dim[[2]]$name
latname = test$var[[1]]$dim[[3]]$name

lon = ncvar_get(test,lonname)
lat = ncvar_get(test,latname)

varname = test$var[[1]]$name

#timelength = length(test$var[[1]]$dim[[3]]$vals)

timein = c(timein,ncvar_get(test,test$var[[1]]$dim[[1]]$name))
timeunits = test$var[[1]]$dim[[1]]$units
spunits = test$var[[1]]$dim[[2]]$units
missingval = test$var[[1]]$missval

varunits = test$var[[1]]$units

vardata = ncvar_get(test,varname)

nc_close(test)

vardata = vardata[,,length(lat):1]
lat = rev(lat)

lonidx = which(lon<=lonrange[2] & lon>=lonrange[1])
latidx = which(lat<=latrange[2] & lat>=latrange[1])

timeidx = 1:dim(vardata)[1]

if(f==1){
timeidxin = timeidx
maxidx = max(timeidxin)
} else {
 timeidxin = timeidx+maxidx
maxidx = max(timeidxin)
}

for(i in 1:length(timeidx)) newmatrix[,,timeidxin[i]] = vardata[timeidx[i],lonidx,latidx]

message("Finished subseting file ",f," / ",length(files))
ptm1end = proc.time()-ptm1
message("Time = ",ptm1end[3]," secs")
message("File subset complete at ",Sys.time())

}



#message("Regridded Time slice ",j," / ",timelength)
#ptmend = proc.time()-ptm
#message("Time = ",ptmend[3]," secs")

###########
# Write new netcdf for Regridded data

library(ncdf)

dimX <- dim.def.ncdf( "lon", "degrees_east", lon[lonidx])
dimY <- dim.def.ncdf( "lat", "degrees_north", lat[latidx])

dimT <- dim.def.ncdf( "time",timeunits,timein)
# Make varables of various dimensionality, for illustration purposes
mv <- missingval # missing value to use

var1d <- var.def.ncdf(varname, varunits, list(dimX,dimY,dimT), mv )

fileout = paste("idaho_v2/",variable,".nc",sep="")

# Create the file
nc <- create.ncdf(fileout , var1d)

# Write some data to the file
put.var.ncdf( nc, var1d, newmatrix ) # no start or count: write all values\

# close ncdf
close.ncdf(nc)

message("File subset complete")



