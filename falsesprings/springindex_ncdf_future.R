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

tasmaxfile = "/data2/3to5/I35/tasmax/EDQM/tasmax_day_I35txp1-EDQM-A28L01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc"
tasminfile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A28L01K00_rcp85_r1i1p1_I35Land_20060101-20991231.nc"

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

save(list=c("LFSEIout","EFSEIout","FLIout","FBIout","DMGout","EFSout","LFSout","lastfreezeout","years","lat","lon"),file="EDQM-A28L01K00_rcp85_SIxresults.Rdata")


####










