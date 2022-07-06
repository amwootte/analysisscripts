library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("/data2/3to5/I35/scripts/colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maptools)
library(ggplot2)
library(zoo)
library(lars)

streamflowID = gaugenum = "08151500"
futureperiod = c(2070,2099)

####

load(file=paste("/home/woot0002/streamflowreg_3daymovingmean_",streamflowID,".Rdata",sep=""))
datesin = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

######

  diff = fulldata2$streamflowmod-fulldata2$streamflow_mean
  diffsq = diff^2
  
  RMSE = sqrt(mean(diffsq,na.rm=TRUE))

RMSE

####

histdattabyearly = aggregate(fulldata2[,c(which(names(fulldata2)=="streamflow_mean"),which(names(fulldata2)=="streamflowmod"))],by=list(year=fulldata2$YEAR),mean,na.rm=TRUE)

  diff = histdattabyearly$streamflowmod-histdattabyearly$streamflow_mean
  diffsq = diff^2
  
  RMSEyearly = sqrt(mean(diffsq,na.rm=TRUE))

  RMSEyearly
  
####

# climatology of yearly average streamflow ensemble RMSE
climotab = apply(histdattabyearly[,2:ncol(histdattabyearly)],2,mean,na.rm=TRUE)

error = climotab[2]-climotab[1]

error

sqrt(mean(error^2))


