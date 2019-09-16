library(ncdf4)
source("/data2/3to5/I35/scripts/analysisfunctions.R")

varname = "pr"
#loc_name = "Deming"
#loc_lon = 360-107.7586
#loc_lat = 32.2687

#loc_name = "SanAntonio"
#loc_lon = 360-106.8659
#loc_lat = 33.9178

#loc_name = "Albuquerque"
#loc_lon = 360-106.6504
#loc_lat = 35.0844


obsfile = "/data2/3to5/I35/pr/Daymet/pr_day_daymet_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

blocksize = 3

varin = varname

#############
# 1a- historical calcs

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
years = 1981:2005
 
  test = nc_open(obsfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  
    ###
    # create model grid
    
    LON = rep(lon,each=length(lat))
    LAT = rep(lat,length(lon))
    R = rep(1:length(lon),each=length(lat))
    C = rep(1:length(lat),length(lon))
    modelgrid = data.frame(R,C,LON,LAT)
    names(modelgrid) = c("R","C","lon","lat")
    if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
    
    ###
    # get cells to use
    loc_lon = as.numeric(loc_lon)
    loc_lat = as.numeric(loc_lat)
    if(loc_lon>0) loc_lon=loc_lon-360
    
    pointarea = distfunc(loc_lon,loc_lat,modelgrid)
    locstart = c(pointarea$R-floor(blocksize/2),pointarea$C-floor(blocksize/2))
    locend = c(pointarea$R+floor(blocksize/2),pointarea$C+floor(blocksize/2))
    
  test = nc_open(obsfile)
    tempdata = ncvar_get(test,varin,start=c(locstart[1],locstart[2],1),count=c(blocksize,blocksize,-1))
    nc_close(test)
    tmp = apply(tempdata,3,mean,na.rm=TRUE)*86400
    
    #### dates to use
    if(length(dates)>length(tmp)){
      rmidx = which(substr(dates,6,10)=="02-29")
      histdates = dates[-rmidx]
    } 
    
    if(length(dates)==length(tmp)){
      rmidx = which(substr(dates,6,10)=="02-29")
      histdates = dates[-rmidx]
      tmp = tmp[-rmidx]
    }

    histyearframe = array(NA,dim=c(365,ncol=length(years)))
    
    for(y in 1:length(years)){
      yearidx = which(as.numeric(substr(histdates,1,4))==years[y])
      histyearframe[1:365,y]= tmp[yearidx]
    }


#####

histyearmean = apply(histyearframe,2,mean,na.rm=TRUE) 
histclimo = mean(histyearmean,na.rm=TRUE) 
histyears = 1981:2005
histyearanom = histyearframe-histclimo

histcumanom2 = histyearanom
for(i in 1:length(histcumanom)){
  if(i==1){
    histcumanom[1] = histyearanom[1]
  } else {
    histcumanom[i] = histyearanom[i]+histcumanom[(i-1)]
  }
}

histcumanom = histyearanom
for(y in 1:length(years)){
for(i in 1:365){
  if(i==1){
    histcumanom[1,y] = histyearanom[1,y]
  } else {
    histcumanom[i,y] = histyearanom[i,y]+histcumanom[(i-1),y]
  }
}
}


#plot(histyearframe[,1],type="l")
#plot(histcumanom[,1],type="l")
#plot(histcumanom[,2],type="l")

##########

histdailymean = apply(histyearframe,1,mean,na.rm=TRUE) 
histdailyanom = histdailymean-histclimo

histdailycumanom = c()
for(i in 1:length(histdailyanom)){
  if(i==1){
    histdailycumanom[1]=histdailyanom[1]
  } else {
    histdailycumanom[i] = histdailyanom[i]+histdailycumanom[(i-1)]
  }
}

plot(histdailycumanom,type="l",main="PRISM")

save(list=c("histdailycumanom"),file=paste("wetseasontesting_Daymet_",loc_name,".Rdata",sep=""))






