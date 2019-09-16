source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

setwd("/home/woot0002/")

PARMdtrfile = "dtr_day_I35dtrdetrp1-PARM-B10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"

EDQMtminfile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc" 
EDQMtmaxfile = "/data2/3to5/I35/tasmax/EDQM/tasmax_day_I35txp1-EDQM-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"

PRISMtminfile = "/data2/3to5/I35/tasmin/PRISM/tasmin_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
PRISMtmaxfile = "/data2/3to5/I35/tasmax/PRISM/tasmax_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

test = nc_open(PARMdtrfile)
PARMdtr = ncvar_get(test,"dtr")
lat = ncvar_get(test,"lat")
lon = ncvar_get(test,"lon")
nc_close(test)

PARMdtr = ifelse(PARMdtr<=0,NA,PARMdtr)


histdates= seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

if(length(histdates)>dim(PARMdtr)[3]){
  histdates2 = histdates[-which(substr(histdates,6,10)=="02-29")]
}

test = nc_open(EDQMtmaxfile)
EDQMtmax = ncvar_get(test,"tasmax")
nc_close(test)

test = nc_open(EDQMtminfile)
EDQMtmin = ncvar_get(test,"tasmin")
nc_close(test)

EDQMdtr = EDQMtmax-EDQMtmin
EDQMdtr = ifelse(EDQMdtr<=0,NA,EDQMdtr)

test = nc_open(PRISMtmaxfile)
PRISMtmax = ncvar_get(test,"tasmax")
nc_close(test)

test = nc_open(PRISMtminfile)
PRISMtmin = ncvar_get(test,"tasmin")
nc_close(test)

PRISMdtr = PRISMtmax-PRISMtmin
PRISMdtr = ifelse(PRISMdtr<=0,NA,PRISMdtr)

library(e1071)

PARMskew = EDQMskew = PRISMskew = matrix(NA,nrow=length(lon),ncol=length(lat))

for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    
    if(is.na(PARMdtr[r,c,1])==FALSE){
      PARMskew[r,c] = skewness(PARMdtr[r,c,],na.rm=TRUE)
      EDQMskew[r,c] = skewness(EDQMdtr[r,c,],na.rm=TRUE)
      PRISMskew[r,c] = skewness(PRISMdtr[r,c,],na.rm=TRUE)
    }
  }
  message("Finished calcs for row ",r," / ",length(lon))
}


diffcolorbar = colorramp(c(PARMskew,EDQMskew,PRISMskew),colorchoice="bluetored",type="difference",Blimit=50)

testsfc1 = list(x=(lon-360),y=lat,z=PARMskew)
surface(testsfc1,type="I",main="PARM dtr skew",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMskew)
surface(testsfc1,type="I",main="EDQM dtr skew",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMskew)
surface(testsfc1,type="I",main="PRISM dtr skew",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

