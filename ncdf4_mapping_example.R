source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

EDQMpastfile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
EDQMfuturefile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A18P01K00_rcp85_r6i1p1_I35Land_20060101-20991231.nc"

con = nc_open(EDQMpastfile)

con

names(con)
names(con$var)
names(con$dim)

names(con$var$tasmin)
names(con$var[[6]])
names(con$var[["tasmin"]])

con$var$tasmin$units

tminhist = ncvar_get(con,"tasmin")
lat = ncvar_get(con,"lat")
lon = ncvar_get(con,"lon")
nc_close(test)

con = nc_open(EDQMfuturefile)
tminfut = ncvar_get(con,"tasmin")
nc_close(test)

histdates= seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
futdates= seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

if(length(histdates)>dim(tminhist)[3]){ # date correction if the model data doesn't have leap days
  histdates2 = histdates[-which(substr(histdates,6,10)=="02-29")]
}

if(length(futdates)>dim(tminfut)[3]){
  futdates2 = futdates[-which(substr(futdates,6,10)=="02-29")]
}

mask = ifelse(is.na(tminhist[,,1])==FALSE,1,0) # sometimes with aggregation you'll end up with NAs being replaced with zeros. This will make a mask you can use to correct that.
earlyidx = which(as.numeric(substr(futdates2,1,4))>=2006 & as.numeric(substr(futdates2,1,4))<=2039)
mididx = which(as.numeric(substr(futdates2,1,4))>=2040 & as.numeric(substr(futdates2,1,4))<=2069)
lateidx = which(as.numeric(substr(futdates2,1,4))>=2070 & as.numeric(substr(futdates2,1,4))<=2099)

histclimo = apply(tminhist,c(1,2),mean,na.rm=TRUE)
earlyclimo = apply(tminfut[,,earlyidx],c(1,2),mean,na.rm=TRUE)
midclimo = apply(tminfut[,,mididx],c(1,2),mean,na.rm=TRUE)
lateclimo = apply(tminfut[,,lateidx],c(1,2),mean,na.rm=TRUE)

earlychange = earlyclimo - histclimo
midchange = midclimo - histclimo
latechange = lateclimo - histclimo

zlims = c(255,290)
brks = seq(255,290,by=5)
colbar = colorRampPalette(c("darkblue","blue","green","yellow","orange","red","darkred"))(length(brks)-1)
sfc = list(x=(lon-360),y=lat,z=tminhist[,,1])
surface(sfc,type="I",main="Tmin Jan 1, 1981",zlim=zlims,col=colbar,breaks=brks,xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

rawcolorbar = colorramp(c(histclimo),type="raw",colorchoice="yellowtored",Blimit=30)
testsfc1 = list(x=(lon-360),y=lat,z=histclimo)
surface(testsfc1,type="I",main="historical Tmin (1981-2005)",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

diffcolorbar = colorramp(c(earlychange,midchange,latechange),type="difference",colorchoice="bluetored",Blimit=30) # This is a custom made function I built to generate the colorbars. Otherwise R defaults to 64 colors in a rainbow ramp.
testsfc1 = list(x=(lon-360),y=lat,z=earlychange)
surface(testsfc1,type="I",main="Projected Change Tmin (2006-2039 minus 1981-2005)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=midchange)
surface(testsfc1,type="I",main="Projected Change Tmin (2040-2069 minus 1981-2005)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=latechange)
surface(testsfc1,type="I",main="Projected Change Tmin (2070-2099 minus 1981-2005)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

