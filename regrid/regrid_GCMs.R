#################################
#
# Downscaled Data Re-gridding Script - GCMs
#
# Adrienne Wootten 11/29/2018
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

args = c("GCMs","pr")

dataset = args[1]
var = args[2]

files = system(paste("ls /home/woot0002/GCMs/",var,"*.nc",sep=""),intern=TRUE)

#files
#####################

library(ncdf4)
library(sp)
library(fields)

####
# get 3^5 lons and lats for new grid
gridfile = "/data2/3to5/I35/pr/EDQM/pr_day_I35prp1-QDM-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
nctest = nc_open(gridfile)
# Select correct file based on the above checks
newlat=ncvar_get(nctest,"lat")
newlon=ncvar_get(nctest,"lon")
newrownum = length(newlon) # number of cells in x direction
newcolnum = length(newlat)  # number of cells in y direction
nc_close(nctest)

for(f in 1:length(files)){
  #for(f in 1:10){
  ptm1 = proc.time()
  
  message("File regrid began at ",Sys.time()," for file: ",files[f])
  
  # read file
  test = nc_open(files[f])
  lonname = test$var[[1]]$dim[[1]]$name
  latname = test$var[[1]]$dim[[2]]$name
  missingval = test$var[[1]]$missval
  
  if(var=="pr"){
    varname1 = "prclimo"
    varname2 = "pr50climo"
    units1 = "mm"
    units2 = "days"
  }
  
  if(var=="tasmax"){
    varname1 = "tmaxclimo"
    varname2 = "tmax95climo"
    units1 = "degrees_C"
    units2 = "days"
  }
  
  if(var=="tasmin"){
    varname1 = "tminclimo"
    varname2 = "tmin32climo"
    units1 = "degrees_C"
    units2 = "days"
  }
  
  timelength = length(test$var[[1]]$dim[[3]]$vals)
  times = ncvar_get(test,test$var[[1]]$dim[[3]]$name)
  timeunits = test$var[[1]]$dim[[3]]$units
  
  spunits = test$var[[1]]$dim[[2]]$units
  
  lon = ncvar_get(test,lonname)
  lat = ncvar_get(test,latname)
  
  if(f==1 & all(lon<0)==TRUE){
    newlon = newlon-360
  }
  
  testvar1 = ncvar_get(test,varname1)
  testvar2 = ncvar_get(test,varname2)
  nc_close(test)
  
  #lat2= lat
  #lon2 =lon
  
  newmatrix1 = newmatrix2 = array(data=NA,dim=c(newrownum,newcolnum,timelength))
  
  for(j in 1:timelength){
    
    ptm = proc.time()
    
    int1 = interp.surface.grid(list(x=lon,y=lat,z=testvar1[,,j]),grid.list=list(x=newlon,y=newlat))
    int2 = interp.surface.grid(list(x=lon,y=lat,z=testvar2[,,j]),grid.list=list(x=newlon,y=newlat))
    
    newvar = matrix(int1[[3]],nrow=length(newlon),ncol=length(newlat))
    if(varname1=="pr"){
      newvar = ifelse(newvar<0,0,newvar)
    }
    newmatrix1[,,j] = newvar
    
    newvar = matrix(int2[[3]],nrow=length(newlon),ncol=length(newlat))
    newvar = ifelse(newvar<0,0,newvar)
    newmatrix2[,,j] = newvar
    
    message("Regridded Time slice ",j," / ",timelength)
    ptmend = proc.time()-ptm
    message("Time = ",ptmend[3]," secs")
    
  }
  
  #testsfc = list(x=lon,y=lat,z=testvar1[,,1])
  #testsfc2 = list(x=newlon,y=newlat,z=newmatrix1[,,1])
  
  #surface(testsfc,type="I")
  #map("state",add=TRUE)
  
  #surface(testsfc2,type="I")
  #map("state",add=TRUE)
  
  ###########
  # Write new netcdf for Regridded data
  
  if(length(grep("degrees",spunits))>0){
    dimX <- ncdim_def( "lon", "degrees_west",newlon)
    dimY <- ncdim_def( "lat", "degrees_north", newlat)
  } else {
    dimX <- ncdim_def( "x", spunits, newlon)
    dimY <- ncdim_def( "y", spunits, newlat )
  }
  
  dimT <- ncdim_def( "time",timeunits,times)
  # Make varables of various dimensionality, for illustration purposes
  mv <- missingval # missing value to use
  
  var1d <- ncvar_def(varname1,units1, list(dimX,dimY,dimT), mv )
  var4d <- ncvar_def(varname2,units2, list(dimX,dimY,dimT), mv )
  
  # Create the test file
  fileout = strsplit(files[f],"/",fixed=TRUE)[[1]][5]
  
  nc <- nc_create(paste("/home/woot0002/GCMs/regrid/",fileout,sep="") ,  list(var1d,var4d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, newmatrix1) # no start or count: write all values\
  ncvar_put(nc, var4d, newmatrix2) # no start or count: write all values\
  
  # close ncdf
  nc_close(nc)
  
  message("Finished regridding file ",f," / ",length(files))
  ptm1end = proc.time()-ptm1
  message("Time = ",ptm1end[3]," secs")
  message("File regrid complete at ",Sys.time())
  
}
