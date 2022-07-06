library(ncdf4)
library(maps)
library(fields)
library(sp)

###
# Daymet check

filename = "/data2/3to5/I35/tasmax/Daymet/tasmax_day_daymet_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

test = nc_open(filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
times = ncvar_get(test,"time")
tasmax = ncvar_get(test,"tasmax")
nc_close(test)

latunits = "degrees_north"
latlongname = "latitude"
lonunits = "degrees_east"
lonlongname = "longitude"
timeunits = "days since 1961-01-01 00:00"
tasmaxunits = "K"
tasmaxlongname = "Daily Maximum Near-Surface Air Temperature"


dimX <- ncdim_def( "lon",lonunits, lon)
dimY <- ncdim_def( "lat",latunits, lat)
dimH <- ncdim_def("time",timeunits,times)

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

var1d <- ncvar_def("tasmax",tasmaxunits,longname=tasmaxlongname, list(dimX,dimY,dimH), mv )

#######
# Create netcdf file

nc <- nc_create("/home/woot0002/Daymettasmax_forcedncdf4.nc" ,list(var1d),force_v4=TRUE )

# Write some data to the file
ncvar_put(nc, var1d, tasmax) # no start or count: write all values\

# close ncdf
nc_close(nc)






#nc <- create.nc("/home/woot0002/Daymettasmax_classicformat.nc",format="classic")
#dim.def.nc(nc, "lon", length(lon))
#dim.def.nc(nc, "lat", length(lat))
#dim.def.nc(nc, "time", length(times), unlim=TRUE)
#dim.def.nc(nc, "max_string_length", 32)

## Create three variables, one as coordinate variable
#var.def.nc(nc, "time", "NC_DOUBLE", "time")
#var.def.nc(nc, "lon", "NC_DOUBLE", "lon")
#var.def.nc(nc, "lat", "NC_DOUBLE", "lat")
#var.def.nc(nc, "tasmax", "NC_FLOAT", c("lon","lat","time"))

#att.put.nc(nc,"lon","long_name","NC_CHAR",lonlongname)
#att.put.nc(nc,"lat","long_name","NC_CHAR",latlongname)
#att.put.nc(nc,"lon","units","NC_CHAR",lonunits)
#att.put.nc(nc,"lat","units","NC_CHAR",latunits)

#att.put.nc(nc,"time","long_name","NC_CHAR","time")
#att.put.nc(nc,"time","calendar","NC_CHAR","julian")
#att.put.nc(nc,"time","units","NC_CHAR",timeunits)

#att.put.nc(nc,"tasmax","long_name","NC_CHAR",tasmaxlongname)
#att.put.nc(nc,"tasmax","units","NC_CHAR",tasmaxunits)
#att.put.nc(nc,"tasmax","standard_name","NC_CHAR","air_temperature")
#att.put.nc(nc, "tasmax", "_FillValue", "NC_FLOAT", 1e+20)
#att.put.nc(nc, "tasmax", "missing_value", "NC_FLOAT", 1e+20)

## Put subsets of the data:
#var.put.nc(nc, "time", times,start=1,count=length(times))
#var.put.nc(nc, "lon", lon,start=1,count=length(lon))
#var.put.nc(nc, "lat", lat,start=1,count=length(lat))
#var.put.nc(nc, "tasmax", tasmax, start=c(1,1,1),count=c(length(lon),length(lat),length(times)))

#close.nc(nc)


