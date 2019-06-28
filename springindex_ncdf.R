####
# Spring Phenology R test code - netcdf
# AMW 6/14/2019
#
# This script demonstrates appropriate use for calculations of first leaf, first bloom, damage index, and false springs using the functions
# in springpheno.R. Specific to vector analyses.
#
#  Questions on code should be emailed to amwootte@ou.edu

source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
setwd("/home/woot0002") # edit or remove paths for scripts and files as needed
source("scripts/springpheno.R")

tasmaxfile = "/data2/3to5/I35/tasmax/EDQM/tasmax_day_I35txp1-EDQM-A38L01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc"
tasminfile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A38L01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc"

dates= seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
years = 2006:2099


# EDQM Future files
# tasmax_day_I35txp1-EDQM-A38L01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc
# tasmin_day_I35tnp1-EDQM-A38L01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc

# EDQM historical files
# tasmax_day_I35txp1-EDQM-A30L01K00_historical_r1i1p1_I35Land_19810101-20051231.nc
# tasmin_day_I35tnp1-EDQM-A30L01K00_historical_r1i1p1_I35Land_19810101-20051231.nc

# Livneh files
# tasmax_day_livneh_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc
# tasmin_day_livneh_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc

test = nc_open(tasmaxfile) # getting lat, lon, and creating a mask to run only points where there are data.
mask = ncvar_get(test,"tasmax",start=c(1,1,1),count=c(-1,-1,1))
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)
pointstorun = which(is.na(mask)==FALSE,arr.ind=TRUE)


FLIout = FBIout = DMGout = array(NA,dim=c(length(lon),length(lat),length(years),4))
EFSout = LFSout = lastfreezeout = array(NA,dim=c(length(lon),length(lat),length(years)))
EFSEIout = LFSEIout = array(NA,dim=c(length(lon),length(lat)))

for(i in 1:nrow(pointstorun)){
    
    R = pointstorun[i,1]
    C = pointstorun[i,2]
   
    test = nc_open(tasmaxfile)
    tasmax = ncvar_get(test,"tasmax",start=c(R,C,1),count=c(1,1,-1))
    tasmax = (tasmax-273.15)*9/5+32
    nc_close(test)
    
    test = nc_open(tasminfile)
    tasmin = ncvar_get(test,"tasmin",start=c(R,C,1),count=c(1,1,-1))
    tasmin = (tasmin-273.15)*9/5+32
    nc_close(test)
    
      TMAX = TMIN = matrix(NA,nrow=366,ncol=length(years))
      for(y in 1:length(years)){
        idx = which(as.numeric(substr(dates,1,4))==years[y])
        if(length(idx)==366){
          TMAX[,y] = tasmax[idx]
          TMIN[,y] = tasmin[idx]
        }
        if(length(idx)==365){
          TMAX[c(1:59,61:366),y] = tasmax[idx]
          TMIN[c(1:59,61:366),y] = tasmin[idx]
        }
      }
      
      RESULTS = calc_si(TMAX,TMIN,lat[C])
      
      EFSout[R,C,] = RESULTS$FSmat[,1]
      LFSout[R,C,] = RESULTS$FSmat[,2]
      lastfreezeout[R,C,] = RESULTS$lastfreeze
      
      EFSEIout[R,C] = RESULTS$FSEImat[1]
      LFSEIout[R,C] = RESULTS$FSEImat[2]
      
      FLIout[R,C,,1] = RESULTS$FLImat[,1]
      FLIout[R,C,,2] = RESULTS$FLImat[,2]
      FLIout[R,C,,3] = RESULTS$FLImat[,3]
      FLIout[R,C,,4] = RESULTS$FLImat[,4]
      
      FBIout[R,C,,1] = RESULTS$FBImat[,1]
      FBIout[R,C,,2] = RESULTS$FBImat[,2]
      FBIout[R,C,,3] = RESULTS$FBImat[,3]
      FBIout[R,C,,4] = RESULTS$FBImat[,4]
      
      DMGout[R,C,,1] = RESULTS$DMGmat[,1]
      DMGout[R,C,,2] = RESULTS$DMGmat[,2]
      DMGout[R,C,,3] = RESULTS$DMGmat[,3]
      DMGout[R,C,,4] = RESULTS$DMGmat[,4]
    
    message("Calculations complete for point ",i," / ",nrow(pointstorun))
}

plot(FLIout[95,70,,1]~years,type="b",lwd=2,pch=19,ylim=range(FLIout[95,70,,]))
lines(FLIout[95,70,,2]~years,col="blue")
lines(FLIout[95,70,,3]~years,col="purple")
lines(FLIout[95,70,,4]~years,col="magenta")

plot(FBIout[95,70,,1]~years,type="b",lwd=2,pch=19,ylim=range(FBIout[95,70,,]))
lines(FBIout[95,70,,2]~years,col="blue")
lines(FBIout[95,70,,3]~years,col="purple")
lines(FBIout[95,70,,4]~years,col="magenta")

plot(DMGout[95,70,,1]~years,type="b",lwd=2,pch=19,ylim=range(DMGout[95,70,,]))
lines(DMGout[95,70,,2]~years,col="blue")
lines(DMGout[95,70,,3]~years,col="purple")
lines(DMGout[95,70,,4]~years,col="magenta")

plot(lastfreezeout[95,70,]~years,type="b",lwd=2,pch=19,ylim=range(lastfreezeout[95,70,]))

FLIclimo = apply(FLIout,c(1,2,4),mean,na.rm=TRUE)
FBIclimo = apply(FBIout,c(1,2,4),mean,na.rm=TRUE)
DMGclimo = apply(DMGout,c(1,2,4),mean,na.rm=TRUE)

par(mfrow=c(2,2))
colorbar = colorramp(as.vector(FLIclimo),colorchoice="yellowtored",type="raw",Blimit=30,use_fixed_scale = FALSE,fixed_scale=NA)
testsfc1 = list(x=lon-360,y=lat,z=FLIclimo[,,1])
surface(testsfc1,type="I",main="Mean FLI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=FLIclimo[,,2])
surface(testsfc1,type="I",main="Plant 1 FLI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=FLIclimo[,,3])
surface(testsfc1,type="I",main="Plant 2 FLI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=FLIclimo[,,4])
surface(testsfc1,type="I",main="Plant 3 FLI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)



par(mfrow=c(2,2))
colorbar = colorramp(as.vector(FBIclimo),colorchoice="yellowtored",type="raw",Blimit=30,use_fixed_scale = FALSE,fixed_scale=NA)
testsfc1 = list(x=lon-360,y=lat,z=FBIclimo[,,1])
surface(testsfc1,type="I",main="Mean FBI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=FBIclimo[,,2])
surface(testsfc1,type="I",main="Plant 1 FBI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=FBIclimo[,,3])
surface(testsfc1,type="I",main="Plant 2 FBI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=FBIclimo[,,4])
surface(testsfc1,type="I",main="Plant 3 FBI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)


par(mfrow=c(2,2))
colorbar = colorramp(as.vector(DMGclimo),colorchoice="yellowtored",type="raw",Blimit=30,use_fixed_scale = FALSE,fixed_scale=NA)
testsfc1 = list(x=lon-360,y=lat,z=DMGclimo[,,1])
surface(testsfc1,type="I",main="Mean DMG - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=DMGclimo[,,2])
surface(testsfc1,type="I",main="Plant 1 DMG - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=DMGclimo[,,3])
surface(testsfc1,type="I",main="Plant 2 DMG - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=DMGclimo[,,4])
surface(testsfc1,type="I",main="Plant 3 DMG - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)


par(mfrow=c(1,2))
colorbar = colorramp(c(EFSEIout,LFSEIout),colorchoice="yellowtored",type="raw",Blimit=20,use_fixed_scale = FALSE,fixed_scale=NA)
testsfc1 = list(x=lon-360,y=lat,z=EFSEIout)
surface(testsfc1,type="I",main="Early FSEI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=LFSEIout)
surface(testsfc1,type="I",main="Late FSEI - Livneh based",zlim=colorbar[[1]],col=colorbar[[3]],breaks=colorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

save(list=c("LFSEIout","EFSEIout","FLIout","FBIout","DMGout","EFSout","LFSout","lastfreezeout","years","lat","lon"),file="EDQMF_SIxresults.Rdata")


####










