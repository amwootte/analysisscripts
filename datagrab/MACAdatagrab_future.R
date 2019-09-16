#############################
#
# MACA monthly location grab
# historical GCM only
#
#############################
########
# Load library and functions

library(ncdf4)

distfunc = function(lon,lat,modelgrid){
  #model grid must have the lats and lons in the same format as the location lat and lon. modelgrid should also have r and c per grid point
  dist=3963 * acos(sin(lat/57.2958)*sin(modelgrid$lat/57.2958)+cos(lat/57.2958)*cos(modelgrid$lat/57.2958)*cos((modelgrid$lon/57.2958)-(lon/57.2958)))
  minidx = which(dist==min(dist,na.rm=TRUE))
  modelgrid[minidx,]
}

########
# get model grid information

test = nc_open("/data/static_web/SPTC_Climate_Projections_Data_and_Images/Data/Average_temperatures/Gridpoint_netcdf_all_domain/CMIP5_MACA/CCSM4/historical/Max_Min_temperatures_MACA_CCSM4_historical_1950.nc")
lonmod = ncvar_get(test,"lon")
latmod = ncvar_get(test,"lat")
nc_close(test)

LON = rep(lonmod,each=length(latmod))
LAT = rep(latmod,length(lonmod))
R = rep(1:length(lonmod),each=length(latmod))
C = rep(1:length(latmod),length(lonmod))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")

#######
# Get city locations and set GCM

locations = read.table("/home/woot0002/REUproject/ticklocationdata.csv",header=TRUE,sep=",")
GCM = "CNRM-CM5"

#######
# Grab temperature and precipitation data for each location.

for(i in 1:nrow(locations)){

lat = locations$Latitude[i] # pull out latitude, longitude, and location name
lon = locations$Longitude[i]
locationname = locations$Location[i]

pointarea = distfunc(lon,lat,modelgrid) # find the right model grid point for each location

monthyear = seq(as.Date("2006-01-15"),as.Date("2099-12-15"),by="month") # set time dimension, separate columns for monthyear, month, and year
monthyear = format(monthyear,"%Y-%m")
year = as.numeric(substr(monthyear,1,4))
month = as.numeric(substr(monthyear,6,7))

minTmin=maxTmax=Tmean=Precip = c() # blank vectors for each temperature and precipitation value

for(y in 2006:2099){
  # determine the right netcdf file names
  filename1 = paste("/data/static_web/SPTC_Climate_Projections_Data_and_Images/Data/Average_temperatures/Gridpoint_netcdf_all_domain/CMIP5_MACA/",GCM,"/future/high_emissions/Max_Min_temperatures_MACA_",GCM,"_rcp85_",y,".nc",sep="") 
  filename2 = paste("/data/static_web/SPTC_Climate_Projections_Data_and_Images/Data/Average_temperatures/Gridpoint_netcdf_all_domain/CMIP5_MACA/",GCM,"/future/high_emissions/monthly_average_temperature_SC_US_",GCM,"_rcp85_",y,".nc",sep="")
  filename3 = paste("/data/static_web/SPTC_Climate_Projections_Data_and_Images/Data/Average_precipitation/Gridpoint_netcdf_all_domain/CMIP5_MACA/",GCM,"/future_high/monthly_precip_accum_SC_US_",GCM,"_rcp85_",y,".nc",sep="")
  
  # gather monthly maximum and monthly minimum temperature
  test = nc_open(filename1)
  minTmin = c(minTmin,ncvar_get(test,"monthly_minimum_Tmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1)))
  maxTmax = c(maxTmax,ncvar_get(test,"monthly_maximum_Tmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1)))
  nc_close(test)
  
  # gather monthly mean temperature
  test = nc_open(filename2)
  Tmean = c(Tmean,ncvar_get(test,"average_temperature",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))*9/5+32)
  nc_close(test)
  
  # gather monthly total precipitation
  test = nc_open(filename3)
  Precip = c(Precip,ncvar_get(test,"precip_accumulation",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1)))
  nc_close(test)
  
  message("Completed data grab for year ",y)
}
  
outputdat = data.frame(monthyear,year,month,maxTmax,Tmean,minTmin,Precip) # create data frame from the vectors created above
write.table(outputdat,paste("/home/woot0002/REUproject/",GCM,"_future_climate_",locationname,".csv",sep=""),row.names=FALSE,sep=",") # write out file for location
  message("completed grab for location: ",locationname," ",i," / ",nrow(locations))
}
  