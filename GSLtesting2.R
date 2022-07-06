####################
#
# GSL comparison OKC
# 1) Get tasmax / tasmin from Daymet for OKC
# 2) Calculate GSL using GSLcalc
# 3) calculate GSL from last spring freeze to first fall freeze (three variants with tasmin, 32, 28, and 25)
# 4) redo 2 calculated separately from function GSLcalc
#

library(ncdf4)
library(maps)
library(fields)
library(sp)
library(mailR)
source("/data2/3to5/I35/scripts/analysisfunctions.R")

locationtable = read.table(file="/data2/3to5/I35/scripts/location_name.csv",header=TRUE,sep=",",colClasses=c("character","numeric","numeric"))
locationtable[,3] = as.numeric(locationtable[,3])
locationtable[,2] = as.numeric(locationtable[,2])

tmaxfile = "/data2/3to5/I35/tasmax/DeltaSD/tasmax_day_I35txdetrp1-DeltaSD-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
tminfile = "/data2/3to5/I35/tasmin/DeltaSD/tasmin_day_I35tndetrp1-DeltaSD-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"

tmaxseries = tminseries = c()

years = 1981:2005

  test = nc_open(tmaxfile)
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    LON = rep(lon,each=length(lat))
    LAT = rep(lat,length(lon))
    R = rep(1:length(lon),each=length(lat))
    C = rep(1:length(lat),length(lon))
    modelgrid = data.frame(R,C,LON,LAT)
    names(modelgrid) = c("R","C","lon","lat")
    if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
    
    loc_lon = as.numeric(locationtable[61,3])
    loc_lat = as.numeric(locationtable[61,2])
    if(loc_lon>0) loc_lon=loc_lon-360
    pointarea = distfunc(loc_lon,loc_lat,modelgrid)
    locstart = c(pointarea$R[1]-1,pointarea$C[1]-1)
    locend = c(pointarea$R[1]+1,pointarea$C[1]+1)
  
  varname=test$var[[1]]$name
  dset = ncvar_get(test,varname)
  nc_close(test)
  
  #tmp = dset[,,1]
  #tmp = t(tmp)
  #tmp = tmp[,length(lat):1]
  
  #testsfc = list(x=lon,y=rev(lat),z=tmp)
  #surface(testsfc,type="I")
  ##map("state",add=TRUE)
  #points(x=pointarea$lon[1],y=pointarea$lat[1],pch=19)
  #points(x=lon[pointarea$R[1]],y=rev(lat)[pointarea$C[1]],pch=19,col="pink")
  
  #tmp[pointarea$R[1],pointarea$C[1]]
  #dset[pointarea$C[1],pointarea$R[1],1]
  
  tmaxseries = c(tmaxseries,dset[pointarea$C[1],pointarea$R[1],])
  
  test = nc_open(tminfile)
  varname=test$var[[1]]$name
  dset = ncvar_get(test,varname)
  nc_close(test)
  
  tminseries = c(tminseries,dset[pointarea$C[1],pointarea$R[1],])

histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
if(length(histdates)>length(tmaxseries)){
  datesoutidx = which(substr(histdates,6,10)=="02-29")
  histdates = histdates[-datesoutidx]
}

tempdata = data.frame(histdates,tmaxseries,tminseries)
names(tempdata) = c("date","tmax","tmin")

tempdata$year = as.numeric(substr(tempdata$date,1,4))

####
# read in station data

statdat = read.table("/home/woot0002/okcstationdata.csv",sep=",",header=TRUE)
statdat$TMAX = ((statdat$TMAX-32)*(5/9))+273.15
statdat$TMIN = ((statdat$TMIN-32)*(5/9))+273.15
statdat$year = as.numeric(substr(statdat$DATE,1,4))

###
# GSLcalc run with METDATA and station data

GSL_MD = GSL_stat = c()
tmax_MD = tmin_MD = tmax_stat = tmin_stat=c()
for(y in 1:length(years)){
  
  yearidx = which(tempdata$year==years[y])
  GSL_MD[y] = GSLcalc(tempdata$tmax[yearidx],tempdata$tmin[yearidx],tempdata$date[yearidx],tempdata$date[yearidx[1]])
  tmax_MD[y] = mean(tempdata$tmax[yearidx])
  tmin_MD[y] = mean(tempdata$tmin[yearidx])
  
  yearidx2 = which(statdat$year==years[y])
  GSL_stat[y] = GSLcalc(statdat$TMAX[yearidx2],statdat$TMIN[yearidx2],as.Date(statdat$DATE[yearidx2]),as.Date(statdat$DATE[yearidx2[1]]))
  tmax_stat[y] = mean(statdat$TMAX[yearidx2])
  tmin_stat[y] = mean(statdat$TMIN[yearidx2])
}


