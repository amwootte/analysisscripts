###############
# Time Capsule Point grab for OKC

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##############
# User supplied inputs

varname = "tasmin" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname # don't change this

distfunc = function(coords = c(lat,lon),dataset){
  
  lat = coords[1]
  lon = coords[2]
  
  dist = 3963 * acos(sin( lat / 57.2958 ) * sin( dataset$LATS / 57.2958 ) + cos( lat / 57.2958 ) * cos( dataset$LATS /57.2958 ) * cos( ( dataset$LONS / 57.2958) - ( lon / 57.2958 ) ))
  
  dataset[which(dist==min(dist)),]
  
}

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

###########
# 1. Data Gather and conversion

filename = paste("/data2/3to5/I35/",varin,"/Livneh/",varin,"_day_livneh_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc",sep="")

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
dateidx = which(substr(dates,6,10)=="02-29")
  
  test = nc_open(filename)
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    rows = 1:length(lon) # row number is the index for longitudes
    cols = 1:length(lat) # columns number is the index for latitude
    
    R = rep(rows,each=length(cols)) # This whole section gets lat, lon, rows, and cols, from a vector format to a table with the grids together
    C = rep(cols,length(rows))
    LONS = rep(lon,each=length(cols))
    LATS = rep(lat,length(rows))
    
    gridmeta = data.frame(R,C,LONS,LATS) 
    citylon = -97.5164 # this happens to be coordinates for Oklahoma City, OK
    citylat = 35.4676  # longitude values need to be in degrees W
    
    closest_point = distfunc(coords=c(citylat,citylon),gridmeta) # the result is a table with one row with the indices, lat, and lon for the closest grid point match. Included is the distance in miles.
    r=closest_point$R[1]
    c=closest_point$C[1]
  
  temp = ncvar_get(test,varin,start=c(r,c,1),count=c(1,1,-1))
  nc_close(test)
  
  histlist = temp
  
################
# Turn everything into data frames

dateshist = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

histlist = data.frame(histlist)
names(histlist)= "Livneh"
histlist = cbind(dateshist,histlist)

#######
# save output

save(list=c("varname","histlist"),file=paste(varname,"_Livneh_output.Rdata",sep=""))



