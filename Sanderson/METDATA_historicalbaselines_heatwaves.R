####################
#
# METDATA historical baseline to 3^5 domain
# 
#####################
###
# Arguments set

datasetname = "METDATA" # common name of dataset to start with

###
# load libraries and master list files
library(ncdf4)
library(maps) # these two just to check plotting issues if need be
library(fields)
source("/home/woot0002/scripts/analysisfunctions.R")
latrange = c(26,40)
lonrange = c(251,270)

###
# setting variable names

###########
# Start calculating yearly averages for each file

years = 1981:2005

for(y in 1:length(years)){
    test = nc_open(paste("/data4/data/OBS/METDATA/tasmax/tmmx_",years[y],".nc",sep=""))
    time = ncvar_get(test,"day")
    timeunits = test$var[[1]]$dim[[3]]$units
    startdate = as.Date(substr(timeunits,12,21))
    times = startdate+time
      if(y==1){
        lon = ncvar_get(test,"lon")
        lat = ncvar_get(test,"lat")
        #lat = rev(lat)
        if(lon[1]<0 & lonrange[1]>0) lonrange=lonrange-360
        #if(lon[1]>0 & lonrange[1]>0) lonrange=lonrange-360;lon=lon-360;
        lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
        latidx = which(lat>= latrange[1] & lat<=latrange[2])
        yeararray = array(NA,dim=c(length(lon[lonidx]),length(lat[latidx]),length(years)))
      }
    tmaxarray = tminarray = array(NA,dim=c(length(lonidx),length(latidx),length(times)))
   
    nc_close(test)
    
      for(d in 1:length(times)){
        test1 = nc_open(paste("/data4/data/OBS/METDATA/tasmax/tmmx_",years[y],".nc",sep=""))
        tmp = ncvar_get(test1,"air_temperature",start=c(latidx[1],lonidx[1],d),count=c(length(latidx),length(lonidx),1))
        tmp2 = t(tmp)
        tmp3 = tmp2[,length(latidx):1]
        tmaxarray[,,d] = tmp3
        nc_close(test1)
        
        test2 = nc_open(paste("/data4/data/OBS/METDATA/tasmin/tmmn_",years[y],".nc",sep=""))
        tmp = ncvar_get(test2,"air_temperature",start=c(latidx[1],lonidx[1],d),count=c(length(latidx),length(lonidx),1))
        tmp2 = t(tmp)
        tmp3 = tmp2[,length(latidx):1]
        tminarray[,,d] = tmp3
        nc_close(test2)
        
        rm(tmp)
        rm(tmp2)
        rm(tmp3)
        message("Finished data gather for day ",d," / ",length(times))
      }
      
      
      test = nc_open("/data2/3to5/I35/METDATA/tmmx_q95.nc")
      q95tmax = ncvar_get(test,"tmaxq95")
      nc_close(test)
      test = nc_open("/data2/3to5/I35/METDATA/tmmn_q95.nc")
      q95tmin = ncvar_get(test,"tminq95")
      nc_close(test)
      for(r in 1:length(lonidx)){
        for(c in 1:length(latidx)){
          yeararray[r,c,y] = heatwave.calc(tmaxarray[r,c,],tminarray[r,c,],q95tmax[r,c],q95tmin[r,c])
          message("Calcs complete for r ",r," and c ",c)
        }
      }  
   
    #testsfc = list(x=lon,y=lat,z=temp3[,length(lat):1])
    #surface(testsfc,type="I")
    
    rm(vardata1)
    rm(vardata2)
    rm(tmaxarray)
    rm(tminarray)
    message("Finished data grab for year: ",years[y])
}
    
climovalues = apply(yeararray,c(1,2),mean,na.rm=TRUE)
rm(yeararray)

#############
# Regrid climatology

test=nc_open("/home/woot0002/uncertainty/anomdat/tasmax_MPI-ESM-LR_EDQMP0_rcp85_anom.nc")
newlon = ncvar_get(test,"lon")
newlat = ncvar_get(test,"lat")
temp = ncvar_get(test,"tasmax",start=c(1,1,1),count=c(-1,-1,1))
nc_close(test)

dmask = ifelse(is.na(temp)==FALSE,1,0)

if(all(newlon>0)==TRUE) lon=lon+360

  test2tavg = interp.surface.grid(list(x=lon[lonidx],y=rev(lat[latidx]),z=climovalues),grid.list=list(x=newlon,y=newlat))
  newmatrix = matrix(test2tavg[[3]],nrow=length(newlon),ncol=length(newlat))
  newmatrix = ifelse(dmask==1,newmatrix,NA)
  
  testsfc = list(x=newlon,y=newlat,z=newmatrix)
  surface(testsfc,type="I")
  
lon = newlon
lat = newlat

#############
# Write data with new file and variable names to uncertainty folder

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

fileout = paste("/data2/3to5/I35/METDATA/heatwaves_histclimo.nc",sep="")

var1d <- ncvar_def("heatwaves_climo", "events_per_year", list(dimX,dimY), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , var1d)
ncvar_put( nc1, var1d, newmatrix) # no start or count: write all values\

nc_close(nc1)





















 
