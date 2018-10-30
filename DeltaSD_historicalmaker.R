###########################
# DeltaSD historical file creation
#
# Creating some faux DeltaSD historical data based on the Daymet, Livneh, and PRISM data 
#

library(ncdf4)
library(sp)
library(fields)
load("/home/woot0002/3to5/Mask3to5.Rdata")

filepath = "/data2/3to5/I35/pr/DeltaSD/"
filenames = c("pr_day_I35prp1-DeltaSD-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A10L01K00_historical_r6i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A20D01K00_historical_r1i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A20L01K00_historical_r1i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A20P01K00_historical_r1i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A30D01K00_historical_r1i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A30L01K00_historical_r1i1p1_I35Land_19810101-20051231.nc",
"pr_day_I35prp1-DeltaSD-A30P01K00_historical_r1i1p1_I35Land_19810101-20051231.nc")

varname = "pr"
longvarname = "Precipitation"
varunits = "kg m-2 s-1"
mv = 1.0000000e+20
lonunits = "degrees_east"
latunits = "degrees_north"
timeunits = "days since 1961-01-01 00:00"

for(i in 1:length(filenames)){
  
  if(length(grep("D01",filenames[i]))==1){
    obsfile = "/home/woot0002/3to5/pr_day_daymet_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
    history = "Created to provide filename consistent Daymet data for DeltaSD"
  }
  if(length(grep("L01",filenames[i]))==1){
    obsfile = "/home/woot0002/3to5/pr_day_livneh_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
    history = "Created to provide filename consistent Livneh data for Livneh"
  }
  if(length(grep("P01",filenames[i]))==1){
    obsfile = "/home/woot0002/3to5/pr_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
    history = "Created to provide filename consistent PRISM data for PRISM"
  }
  
  test = nc_open(obsfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  
  outarray = array(NA,dim=c(length(lon),length(lat),length(times)))
  for(j in 1:length(times)){
    vardat = ncvar_get(test,varname,start=c(1,1,j),count=c(-1,-1,1))
    outarray[,,j] = ifelse(mask==1,vardat,NA)
  }
  nc_close(test)
  
  dimX <- ncdim_def("lon", lonunits, lon)
  dimY <- ncdim_def("lat", latunits, lat)
  dimT <- ncdim_def("times",timeunits,times)
  
  var1d <- ncvar_def(varname,varunits,longname=longvarname, list(dimX,dimY,dimT), mv )
  nc <- nc_create(paste(filepath,filenames[i],sep="") ,  list(var1d) )
  ncvar_put(nc, var1d, outarray)
  ncatt_put(nc,0,"history",history)
  nc_close(nc)
  message("Historical file creation complete for filenames: ",filenames[i])
}



