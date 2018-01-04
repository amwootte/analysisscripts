library(ncdf4)
library(maps)
library(fields)

setwd("/home/woot0002")

# Red River
test = nc_open("/data/static_web/RedRiver/Downscaled/MPI-ESM-LR/historical/pr_day_RRprp1-BCQM-A30aaL01K00_historical_r1i1p1_RR_19610101-20051231.nc")
tempdata1 = ncvar_get(test,"pr")
  lat1 = ncvar_get(test,"lat")
  lon1 = ncvar_get(test,"lon")-360
  times1 = ncvar_get(test,"time")
  startdate1 = as.Date(substr(test$var[[5]]$dim[[3]]$units,12,21),"%Y-%m-%d")
  timeunits1 = test$var[[5]]$dim[[3]]$units
  dates1 = startdate1+times1
  domainmask1 = ifelse(is.na(tempdata1[,,1])==FALSE,1,0)
nc_close(test)

# 3^5 domain
test = nc_open("pr_day_I35prp1-EDQM-A38L01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc")
tempdata2 = ncvar_get(test,"pr")
lat2 = ncvar_get(test,"lat")
lon2 = ncvar_get(test,"lon")-360
times2 = ncvar_get(test,"time")
startdate2 = as.Date(substr(test$var[[5]]$dim[[3]]$units,12,21),"%Y-%m-%d")
timeunits2 = test$var[[5]]$dim[[3]]$units
dates2 = startdate2+times2
domainmask2 = ifelse(is.na(tempdata2[,,1])==FALSE,1,0)
nc_close(test)

#############

testsfcRR = list(x=lon1,y=lat1,z=domainmask1)
surface(testsfcRR,type="I")
map("state",add=TRUE)

testsfc = list(x=lon2,y=lat2,z=domainmask2)
surface(testsfc,type="I")
map("state",add=TRUE)

#############
# Create mask netcdf

dimX <- ncdim_def("lon","degrees_east",lon1)
dimY <- ncdim_def("lat","degrees_north",lat1)

var1 <- ncvar_def("RRmask","",dim=list(dimX,dimY),missval=999999)

nc <- nc_create("RRdomainmask.nc",var1)

ncvar_put(nc,var1,domainmask1)
nc_close(nc)
