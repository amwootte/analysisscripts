##############
# URG check

library(ncdf4)
library(maps)
library(fields)
library(sp)
library(mailR)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

datacheck = read.table("/data2/3to5/I35/URG/tasmax/EDQM/future/tasmax_daily_MPI-ESM-LR_EDQM_Daymet_rcp45_2006-2099.txt",sep=" ")

cols = c(772, 779, 781, 783, 785, 794, 824, 826, 910, 911)

datacheck = datacheck[,cols]

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
datesidx = which(substr(dates,1,7)=="2099-12")

datacheck = datacheck[datesidx,]
days= 1:31

datadiff= datacheck[-31,]

for(i in 1:ncol(datacheck)){
  datadiff[,i]=diff(datacheck[,i])
  if(i ==1){
    plot(datacheck[,i]~days,type="l",ylim=range(datacheck))
  } else {
    lines(datacheck[,i]~days,type="l")
  }
}

for(i in 1:ncol(datadiff)){
  if(i ==1){
    plot(diff(datadiff[,i]),type="l",ylim=c(-10,10))
  } else {
    lines(diff(datadiff[,i]),type="l")
  }
}

#####
# netcdf mapping

nctest = nc_open("/data2/3to5/I35/tasmax/EDQM/tasmax_day_I35txp1-EDQM-A34D01K00_rcp45_r1i1p1_I35Land_20060101-20991231.nc")
lon = ncvar_get(nctest,"lon")
lat = ncvar_get(nctest,"lat")
times = ncvar_get(nctest,"time")
tasmax = ncvar_get(nctest,"tasmax",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))
nc_close(nctest)

nctest = nc_open("/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A34D01K00_rcp45_r1i1p1_I35Land_20060101-20991231.nc")
tasmin = ncvar_get(nctest,"tasmin",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))
nc_close(nctest)

nctest = nc_open("/data2/3to5/I35/pr/EDQM/pr_day_I35prp1-QDM-A34D01K00_rcp45_r1i1p1_I35Land_20060101-20991231.nc")
pr = ncvar_get(nctest,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
nc_close(nctest)

tasmaxdiff = tasmax
for(i in 1:31){
  if(i==1){
    tasmaxdiff[,,i]=NA
  } else {
    tasmaxdiff[,,i]=tasmax[,,i]-tasmax[,,(i-1)]
  }
}


tempcolorbar = colorramp(tasmax,colorchoice="yellowtored",Blimit=20,use_fixed_scale = FALSE,fixed_scale=NA,type="raw")
prcolorbar = colorramp(pr,colorchoice="whitetogreen",Blimit=30,use_fixed_scale = TRUE,fixed_scale=c(0,125),type="raw")
tempdiffcolorbar = colorramp(tasmaxdiff,colorchoice="bluetored",Blimit=30,use_fixed_scale = FALSE,fixed_scale=NA)

datesin = dates[datesidx]

pdf("/home/woot0002/Dec2099checks.pdf",width=35,height=15)
par(mfrow=c(3,7))

for(i in 14:20){
  
  testsfc = list(x=lon-360,y=lat,z=tasmax[,,i])
  surface(testsfc,type="I",main=paste("High Temperature (K): ",datesin[i],sep=""),zlim=tempcolorbar[[1]],col=tempcolorbar[[3]],breaks=tempcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(tasmax[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tasmax[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(tasmax[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
}

frame()

for(i in 15:20){
  
  testsfc = list(x=lon-360,y=lat,z=tasmaxdiff[,,i])
  surface(testsfc,type="I",main=paste("Temp Change: ",datesin[(i-1)],"-",datesin[i],sep=""),zlim=tempdiffcolorbar[[1]],col=tempdiffcolorbar[[3]],breaks=tempdiffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(tasmaxdiff[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tasmaxdiff[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(tasmaxdiff[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
}

for(i in 14:20){
  
  testsfc = list(x=lon-360,y=lat,z=pr[,,i])
  surface(testsfc,type="I",main=paste("Precip (mm): ",datesin[i],sep=""),zlim=prcolorbar[[1]],col=prcolorbar[[3]],breaks=prcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(pr[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(pr[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(pr[,,i],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
}
dev.off()