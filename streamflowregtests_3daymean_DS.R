#####
#
# streamflow testing - regression with C-PrEP

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

gaugenum = "08151500"
loclon = -98.669444 # gauging station location 
loclat = 30.751111
futureperiod = c(2070,2099)
runcalcs=FALSE
# load DS files

prfile_hist = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_*historical*.nc",intern=TRUE))

tasmaxfile_hist = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/tasmax/PARM/tasmax_*historical*.nc",intern=TRUE))

tasminfile_hist = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_*historical*.nc",intern=TRUE))

prfile_proj = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_day*rcp*.nc",intern=TRUE))

tasmaxfile_proj = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_day*rcp*.nc",intern=TRUE))

tasminfile_proj = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_day*rcp*.nc",intern=TRUE))

#####
# Create file breakdown tables

filebreakdown = do.call(rbind,strsplit(prfile_proj,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(prfile_hist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

if(runcalcs==TRUE){

####
# find data point to use in all.

nctest = nc_open(prfile_hist[1])
lon = ncvar_get(nctest,"lon")-360
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

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

pointarea = distfunc(loclon,loclat,modelgrid)

######
# Start a loop through historical files.
histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
histresults = list()
for(i in 1:length(prfile_hist)){
  
  ###
  # Get PRCP, TMAX, and TMIN data
  
  nctest = nc_open(prfile_hist[i])
  PRCP = ncvar_get(nctest,"pr",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))*86400
  nc_close(nctest)
  
  nctest = nc_open(tasmaxfile_hist[i])
  TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  nctest = nc_open(tasminfile_hist[i])
  TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  if(length(PRCP)<length(histdates) | length(TMAX)<length(histdates) | length(TMIN)<length(histdates)){
    idxout = which(substr(histdates,6,10)=="02-29")
    datesin = histdates[-idxout]
    if(length(PRCP)==length(histdates)){
      PRCP = PRCP[-idxout]
    }
    if(length(TMAX)==length(histdates)){
      TMAX = TMAX[-idxout]
    }
    if(length(TMIN)==length(histdates)){
      TMIN = TMIN[-idxout]
    }
  } else {
    datesin = histdates
  }
  
  #####
  # unit conversion
  
  PRCP = PRCP/25.4
  TMAX = TMAX-273.15
  TMAX = (TMAX*9/5)+32
  TMIN = TMIN-273.15
  TMIN = (TMIN*9/5)+32
  DATE = datesin
  climate = data.frame(DATE,PRCP,TMAX,TMIN)
  
  #####
  # combine data
  
  fulldata = climate
  
  fulldata$PRCP=as.numeric(fulldata$PRCP)
  fulldata$TMAX=as.numeric(fulldata$TMAX)
  fulldata$TMIN=as.numeric(fulldata$TMIN)
  
  fulldata$TAVG=(as.numeric(fulldata$TMAX)+as.numeric(fulldata$TMIN))/2
  
  fulldata$HOT1S = ifelse(fulldata$TMAX>mean(fulldata$TMAX,na.rm=TRUE)+sd(fulldata$TMAX,na.rm=TRUE),1,0) # hotdays 1 sd, 2 sd, 3 sd, and above 90F
  fulldata$HOT2S = ifelse(fulldata$TMAX>mean(fulldata$TMAX,na.rm=TRUE)+2*sd(fulldata$TMAX,na.rm=TRUE),1,0)
  fulldata$HOT3S = ifelse(fulldata$TMAX>mean(fulldata$TMAX,na.rm=TRUE)+3*sd(fulldata$TMAX,na.rm=TRUE),1,0)
  fulldata$HOT90 = ifelse(fulldata$TMAX>90,1,0)
  
  fulldata$COLD1S = ifelse(fulldata$TMIN<mean(fulldata$TMIN,na.rm=TRUE)+sd(fulldata$TMIN,na.rm=TRUE),1,0) # hotdays 1 sd, 2 sd, 3 sd, and above 90F
  fulldata$COLD2S = ifelse(fulldata$TMIN<mean(fulldata$TMIN,na.rm=TRUE)+2*sd(fulldata$TMIN,na.rm=TRUE),1,0)
  fulldata$COLD3S = ifelse(fulldata$TMIN<mean(fulldata$TMIN,na.rm=TRUE)+3*sd(fulldata$TMIN,na.rm=TRUE),1,0)
  fulldata$COLD32 = ifelse(fulldata$TMIN<32,1,0)
  
  fulldata$TAVG_3D = c(rep(NA,2),rollapply(fulldata$TAVG,width=3,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_1W = c(rep(NA,6),rollapply(fulldata$TAVG,width=7,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_2W = c(rep(NA,13),rollapply(fulldata$TAVG,width=14,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_3W = c(rep(NA,20),rollapply(fulldata$TAVG,width=21,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_4W = c(rep(NA,27),rollapply(fulldata$TAVG,width=28,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_1M = c(rep(NA,29),rollapply(fulldata$TAVG,width=30,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_2M = c(rep(NA,59),rollapply(fulldata$TAVG,width=60,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_3M = c(rep(NA,89),rollapply(fulldata$TAVG,width=90,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_4M = c(rep(NA,119),rollapply(fulldata$TAVG,width=120,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_5M = c(rep(NA,149),rollapply(fulldata$TAVG,width=150,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_6M = c(rep(NA,179),rollapply(fulldata$TAVG,width=180,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_12M = c(rep(NA,364),rollapply(fulldata$TAVG,width=365,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_24M = c(rep(NA,729),rollapply(fulldata$TAVG,width=730,FUN=mean,na.rm=TRUE))
  
  fulldata$PRCPAVG_3D = c(rep(NA,2),rollapply(fulldata$PRCP,width=3,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_1W = c(rep(NA,6),rollapply(fulldata$PRCP,width=7,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_2W = c(rep(NA,13),rollapply(fulldata$PRCP,width=14,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_3W = c(rep(NA,20),rollapply(fulldata$PRCP,width=21,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_4W = c(rep(NA,27),rollapply(fulldata$PRCP,width=28,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_1M = c(rep(NA,29),rollapply(fulldata$PRCP,width=30,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_2M = c(rep(NA,59),rollapply(fulldata$PRCP,width=60,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_3M = c(rep(NA,89),rollapply(fulldata$PRCP,width=90,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_4M = c(rep(NA,119),rollapply(fulldata$PRCP,width=120,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_5M = c(rep(NA,149),rollapply(fulldata$PRCP,width=150,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_6M = c(rep(NA,179),rollapply(fulldata$PRCP,width=180,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_12M = c(rep(NA,364),rollapply(fulldata$PRCP,width=365,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_24M = c(rep(NA,729),rollapply(fulldata$PRCP,width=730,FUN=mean,na.rm=TRUE))
  
  fulldata$PRCP_3D = c(rep(NA,2),rollapply(fulldata$PRCP,width=3,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_1W = c(rep(NA,6),rollapply(fulldata$PRCP,width=7,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_2W = c(rep(NA,13),rollapply(fulldata$PRCP,width=14,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_3W = c(rep(NA,20),rollapply(fulldata$PRCP,width=21,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_4W = c(rep(NA,27),rollapply(fulldata$PRCP,width=28,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_1M = c(rep(NA,29),rollapply(fulldata$PRCP,width=30,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_2M = c(rep(NA,59),rollapply(fulldata$PRCP,width=60,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_3M = c(rep(NA,89),rollapply(fulldata$PRCP,width=90,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_4M = c(rep(NA,119),rollapply(fulldata$PRCP,width=120,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_5M = c(rep(NA,149),rollapply(fulldata$PRCP,width=150,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_6M = c(rep(NA,179),rollapply(fulldata$PRCP,width=180,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_12M = c(rep(NA,364),rollapply(fulldata$PRCP,width=365,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_24M = c(rep(NA,729),rollapply(fulldata$PRCP,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$MAXPRCP1W_2W = c(rep(NA,13),rollapply(fulldata$PRCP_1W,width=14,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_3W = c(rep(NA,20),rollapply(fulldata$PRCP_1W,width=21,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_4W = c(rep(NA,27),rollapply(fulldata$PRCP_1W,width=28,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_1M = c(rep(NA,29),rollapply(fulldata$PRCP_1W,width=30,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_2M = c(rep(NA,59),rollapply(fulldata$PRCP_1W,width=60,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_3M = c(rep(NA,89),rollapply(fulldata$PRCP_1W,width=90,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_4M = c(rep(NA,119),rollapply(fulldata$PRCP_1W,width=120,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_5M = c(rep(NA,149),rollapply(fulldata$PRCP_1W,width=150,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_6M = c(rep(NA,179),rollapply(fulldata$PRCP_1W,width=180,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_12M = c(rep(NA,364),rollapply(fulldata$PRCP_1W,width=365,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_24M = c(rep(NA,729),rollapply(fulldata$PRCP_1W,width=730,FUN=max,na.rm=TRUE))
  
  fulldata$MDRN = ifelse(fulldata$PRCP<0.01,0,1)
  fulldata$MDRN_1W = c(rep(NA,6),rollapply(fulldata$MDRN,width=7,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_2W = c(rep(NA,13),rollapply(fulldata$MDRN,width=14,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_3W = c(rep(NA,20),rollapply(fulldata$MDRN,width=21,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_4W = c(rep(NA,27),rollapply(fulldata$MDRN,width=28,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_1M = c(rep(NA,29),rollapply(fulldata$MDRN,width=30,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_2M = c(rep(NA,59),rollapply(fulldata$MDRN,width=60,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_3M = c(rep(NA,89),rollapply(fulldata$MDRN,width=90,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_4M = c(rep(NA,119),rollapply(fulldata$MDRN,width=120,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_5M = c(rep(NA,149),rollapply(fulldata$MDRN,width=150,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_6M = c(rep(NA,179),rollapply(fulldata$MDRN,width=180,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_12M = c(rep(NA,364),rollapply(fulldata$MDRN,width=365,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_24M = c(rep(NA,729),rollapply(fulldata$MDRN,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$TMAXAVG_1W = c(rep(NA,6),rollapply(fulldata$TMAX,width=7,FUN=mean,na.rm=TRUE))
  
  fulldata$MAXTMAX1W_2W = c(rep(NA,13),rollapply(fulldata$TMAXAVG_1W,width=14,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_3W = c(rep(NA,20),rollapply(fulldata$TMAXAVG_1W,width=21,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_4W = c(rep(NA,27),rollapply(fulldata$TMAXAVG_1W,width=28,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_1M = c(rep(NA,29),rollapply(fulldata$TMAXAVG_1W,width=30,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_2M = c(rep(NA,59),rollapply(fulldata$TMAXAVG_1W,width=60,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_3M = c(rep(NA,89),rollapply(fulldata$TMAXAVG_1W,width=90,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_4M = c(rep(NA,119),rollapply(fulldata$TMAXAVG_1W,width=120,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_5M = c(rep(NA,149),rollapply(fulldata$TMAXAVG_1W,width=150,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_6M = c(rep(NA,179),rollapply(fulldata$TMAXAVG_1W,width=180,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_12M = c(rep(NA,364),rollapply(fulldata$TMAXAVG_1W,width=365,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_24M = c(rep(NA,729),rollapply(fulldata$TMAXAVG_1W,width=730,FUN=max,na.rm=TRUE))
  
  fulldata$TMINAVG_1W = c(rep(NA,6),rollapply(fulldata$TMAX,width=7,FUN=mean,na.rm=TRUE))
  
  fulldata$MINTMIN1W_2W = c(rep(NA,13),rollapply(fulldata$TMINAVG_1W,width=14,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_3W = c(rep(NA,20),rollapply(fulldata$TMINAVG_1W,width=21,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_4W = c(rep(NA,27),rollapply(fulldata$TMINAVG_1W,width=28,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_1M = c(rep(NA,29),rollapply(fulldata$TMINAVG_1W,width=30,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_2M = c(rep(NA,59),rollapply(fulldata$TMINAVG_1W,width=60,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_3M = c(rep(NA,89),rollapply(fulldata$TMINAVG_1W,width=90,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_4M = c(rep(NA,119),rollapply(fulldata$TMINAVG_1W,width=120,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_5M = c(rep(NA,149),rollapply(fulldata$TMINAVG_1W,width=150,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_6M = c(rep(NA,179),rollapply(fulldata$TMINAVG_1W,width=180,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_12M = c(rep(NA,364),rollapply(fulldata$TMINAVG_1W,width=365,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_24M = c(rep(NA,729),rollapply(fulldata$TMINAVG_1W,width=730,FUN=min,na.rm=TRUE))
  
  fulldata$HOT1S_1W = c(rep(NA,6),rollapply(fulldata$HOT1S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_2W = c(rep(NA,13),rollapply(fulldata$HOT1S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_3W = c(rep(NA,20),rollapply(fulldata$HOT1S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_4W = c(rep(NA,27),rollapply(fulldata$HOT1S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_1M = c(rep(NA,29),rollapply(fulldata$HOT1S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_2M = c(rep(NA,59),rollapply(fulldata$HOT1S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_3M = c(rep(NA,89),rollapply(fulldata$HOT1S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_4M = c(rep(NA,119),rollapply(fulldata$HOT1S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_5M = c(rep(NA,149),rollapply(fulldata$HOT1S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_6M = c(rep(NA,179),rollapply(fulldata$HOT1S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_12M = c(rep(NA,364),rollapply(fulldata$HOT1S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_24M = c(rep(NA,729),rollapply(fulldata$HOT1S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$HOT2S_1W = c(rep(NA,6),rollapply(fulldata$HOT2S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_2W = c(rep(NA,13),rollapply(fulldata$HOT2S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_3W = c(rep(NA,20),rollapply(fulldata$HOT2S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_4W = c(rep(NA,27),rollapply(fulldata$HOT2S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_1M = c(rep(NA,29),rollapply(fulldata$HOT2S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_2M = c(rep(NA,59),rollapply(fulldata$HOT2S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_3M = c(rep(NA,89),rollapply(fulldata$HOT2S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_4M = c(rep(NA,119),rollapply(fulldata$HOT2S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_5M = c(rep(NA,149),rollapply(fulldata$HOT2S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_6M = c(rep(NA,179),rollapply(fulldata$HOT2S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_12M = c(rep(NA,364),rollapply(fulldata$HOT2S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_24M = c(rep(NA,729),rollapply(fulldata$HOT2S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$HOT3S_1W = c(rep(NA,6),rollapply(fulldata$HOT3S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_2W = c(rep(NA,13),rollapply(fulldata$HOT3S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_3W = c(rep(NA,20),rollapply(fulldata$HOT3S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_4W = c(rep(NA,27),rollapply(fulldata$HOT3S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_1M = c(rep(NA,29),rollapply(fulldata$HOT3S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_2M = c(rep(NA,59),rollapply(fulldata$HOT3S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_3M = c(rep(NA,89),rollapply(fulldata$HOT3S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_4M = c(rep(NA,119),rollapply(fulldata$HOT3S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_5M = c(rep(NA,149),rollapply(fulldata$HOT3S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_6M = c(rep(NA,179),rollapply(fulldata$HOT3S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_12M = c(rep(NA,364),rollapply(fulldata$HOT3S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_24M = c(rep(NA,729),rollapply(fulldata$HOT3S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$HOT90_1W = c(rep(NA,6),rollapply(fulldata$HOT90,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_2W = c(rep(NA,13),rollapply(fulldata$HOT90,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_3W = c(rep(NA,20),rollapply(fulldata$HOT90,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_4W = c(rep(NA,27),rollapply(fulldata$HOT90,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_1M = c(rep(NA,29),rollapply(fulldata$HOT90,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_2M = c(rep(NA,59),rollapply(fulldata$HOT90,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_3M = c(rep(NA,89),rollapply(fulldata$HOT90,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_4M = c(rep(NA,119),rollapply(fulldata$HOT90,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_5M = c(rep(NA,149),rollapply(fulldata$HOT90,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_6M = c(rep(NA,179),rollapply(fulldata$HOT90,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_12M = c(rep(NA,364),rollapply(fulldata$HOT90,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_24M = c(rep(NA,729),rollapply(fulldata$HOT90,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD1S_1W = c(rep(NA,6),rollapply(fulldata$COLD1S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_2W = c(rep(NA,13),rollapply(fulldata$COLD1S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_3W = c(rep(NA,20),rollapply(fulldata$COLD1S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_4W = c(rep(NA,27),rollapply(fulldata$COLD1S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_1M = c(rep(NA,29),rollapply(fulldata$COLD1S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_2M = c(rep(NA,59),rollapply(fulldata$COLD1S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_3M = c(rep(NA,89),rollapply(fulldata$COLD1S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_4M = c(rep(NA,119),rollapply(fulldata$COLD1S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_5M = c(rep(NA,149),rollapply(fulldata$COLD1S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_6M = c(rep(NA,179),rollapply(fulldata$COLD1S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_12M = c(rep(NA,364),rollapply(fulldata$COLD1S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_24M = c(rep(NA,729),rollapply(fulldata$COLD1S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD2S_1W = c(rep(NA,6),rollapply(fulldata$COLD2S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_2W = c(rep(NA,13),rollapply(fulldata$COLD2S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_3W = c(rep(NA,20),rollapply(fulldata$COLD2S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_4W = c(rep(NA,27),rollapply(fulldata$COLD2S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_1M = c(rep(NA,29),rollapply(fulldata$COLD2S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_2M = c(rep(NA,59),rollapply(fulldata$COLD2S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_3M = c(rep(NA,89),rollapply(fulldata$COLD2S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_4M = c(rep(NA,119),rollapply(fulldata$COLD2S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_5M = c(rep(NA,149),rollapply(fulldata$COLD2S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_6M = c(rep(NA,179),rollapply(fulldata$COLD2S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_12M = c(rep(NA,364),rollapply(fulldata$COLD2S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_24M = c(rep(NA,729),rollapply(fulldata$COLD2S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD3S_1W = c(rep(NA,6),rollapply(fulldata$COLD3S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_2W = c(rep(NA,13),rollapply(fulldata$COLD3S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_3W = c(rep(NA,20),rollapply(fulldata$COLD3S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_4W = c(rep(NA,27),rollapply(fulldata$COLD3S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_1M = c(rep(NA,29),rollapply(fulldata$COLD3S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_2M = c(rep(NA,59),rollapply(fulldata$COLD3S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_3M = c(rep(NA,89),rollapply(fulldata$COLD3S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_4M = c(rep(NA,119),rollapply(fulldata$COLD3S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_5M = c(rep(NA,149),rollapply(fulldata$COLD3S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_6M = c(rep(NA,179),rollapply(fulldata$COLD3S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_12M = c(rep(NA,364),rollapply(fulldata$COLD3S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_24M = c(rep(NA,729),rollapply(fulldata$COLD3S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD32_1W = c(rep(NA,6),rollapply(fulldata$COLD32,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_2W = c(rep(NA,13),rollapply(fulldata$COLD32,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_3W = c(rep(NA,20),rollapply(fulldata$COLD32,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_4W = c(rep(NA,27),rollapply(fulldata$COLD32,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_1M = c(rep(NA,29),rollapply(fulldata$COLD32,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_2M = c(rep(NA,59),rollapply(fulldata$COLD32,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_3M = c(rep(NA,89),rollapply(fulldata$COLD32,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_4M = c(rep(NA,119),rollapply(fulldata$COLD32,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_5M = c(rep(NA,149),rollapply(fulldata$COLD32,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_6M = c(rep(NA,179),rollapply(fulldata$COLD32,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_12M = c(rep(NA,364),rollapply(fulldata$COLD32,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_24M = c(rep(NA,729),rollapply(fulldata$COLD32,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$YEAR = as.numeric(substr(fulldata$DATE,1,4))
  years = unique(fulldata$YEAR)
  
  fulldata2 = fulldata
  
  #####
  
  load(paste("/home/woot0002/streamflowtests_3daymovingmean_",gaugenum,".Rdata",sep=""))
  #load("/home/woot0002/streamflowtests_3daymovingmean_08146000.Rdata")
  
  corsmaxidx = c()
  corsmaxval = c()
  for(I in 1:length(corsusedall)){
    
    tmp = do.call("c",corsusedall[[I]])
    corsmaxidx[I] = which(tmp==max(tmp))
    corsmaxval[I] = max(tmp)
    
  }
  
  maxidx = which(corsmaxval==max(corsmaxval))
  varnamesin = varsusedall[[maxidx]][[corsmaxidx[maxidx]]]
  
  allnames = names(fulldata2)
  usevaridx = which(allnames %in% varnamesin)
  
  #####
  # All vars
  
  xmat = as.matrix(fulldata2[,usevaridx])

  load(paste("/home/woot0002/",histfilebreakdown$obs[i],"streamflowfit_",gaugenum,".Rdata",sep=""))
  #load(paste("/home/woot0002/",histfilebreakdown$obs[i],"streamflowfit_08146000.Rdata",sep=""))
  best_step <- fit$df[which.min(fit$RSS)]
  predictions <- predict(fit, xmat, s=best_step, type="fit")$fit
  predstreamflow = ifelse(predictions<0,0,predictions)
  
  histresults[[i]] = list(predstreamflow)
  message("Finished gathering results for file ",i," / ",length(prfile_hist))
}

######
# Start a loop through historical files.
projdates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
projresults = list()
for(i in 1:length(prfile_proj)){
  
  ###
  # Get PRCP, TMAX, and TMIN data
  
  nctest = nc_open(prfile_proj[i])
  PRCP = ncvar_get(nctest,"pr",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))*86400
  nc_close(nctest)
  
  nctest = nc_open(tasmaxfile_proj[i])
  TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  nctest = nc_open(tasminfile_proj[i])
  TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  if(length(PRCP)<length(projdates) | length(TMAX)<length(projdates) | length(TMIN)<length(projdates)){
    idxout = which(substr(projdates,6,10)=="02-29")
    datesin = projdates[-idxout]
    if(length(PRCP)==length(projdates)){
      PRCP = PRCP[-idxout]
    }
    if(length(TMAX)==length(projdates)){
      TMAX = TMAX[-idxout]
    }
    if(length(TMIN)==length(projdates)){
      TMIN = TMIN[-idxout]
    }
  } else {
    datesin = projdates
  }
  
  #####
  # unit conversion
  
  PRCP = PRCP/25.4
  TMAX = TMAX-273.15
  TMAX = (TMAX*9/5)+32
  TMIN = TMIN-273.15
  TMIN = (TMIN*9/5)+32
  DATE = datesin
  climate = data.frame(DATE,PRCP,TMAX,TMIN)
  
  #####
  # combine data
  
  fulldata = climate
  
  fulldata$PRCP=as.numeric(fulldata$PRCP)
  fulldata$TMAX=as.numeric(fulldata$TMAX)
  fulldata$TMIN=as.numeric(fulldata$TMIN)
  
  fulldata$TAVG=(as.numeric(fulldata$TMAX)+as.numeric(fulldata$TMIN))/2
  
  fulldata$HOT1S = ifelse(fulldata$TMAX>mean(fulldata$TMAX,na.rm=TRUE)+sd(fulldata$TMAX,na.rm=TRUE),1,0) # hotdays 1 sd, 2 sd, 3 sd, and above 90F
  fulldata$HOT2S = ifelse(fulldata$TMAX>mean(fulldata$TMAX,na.rm=TRUE)+2*sd(fulldata$TMAX,na.rm=TRUE),1,0)
  fulldata$HOT3S = ifelse(fulldata$TMAX>mean(fulldata$TMAX,na.rm=TRUE)+3*sd(fulldata$TMAX,na.rm=TRUE),1,0)
  fulldata$HOT90 = ifelse(fulldata$TMAX>90,1,0)
  
  fulldata$COLD1S = ifelse(fulldata$TMIN<mean(fulldata$TMIN,na.rm=TRUE)+sd(fulldata$TMIN,na.rm=TRUE),1,0) # hotdays 1 sd, 2 sd, 3 sd, and above 90F
  fulldata$COLD2S = ifelse(fulldata$TMIN<mean(fulldata$TMIN,na.rm=TRUE)+2*sd(fulldata$TMIN,na.rm=TRUE),1,0)
  fulldata$COLD3S = ifelse(fulldata$TMIN<mean(fulldata$TMIN,na.rm=TRUE)+3*sd(fulldata$TMIN,na.rm=TRUE),1,0)
  fulldata$COLD32 = ifelse(fulldata$TMIN<32,1,0)
  
  fulldata$TAVG_3D = c(rep(NA,2),rollapply(fulldata$TAVG,width=3,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_1W = c(rep(NA,6),rollapply(fulldata$TAVG,width=7,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_2W = c(rep(NA,13),rollapply(fulldata$TAVG,width=14,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_3W = c(rep(NA,20),rollapply(fulldata$TAVG,width=21,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_4W = c(rep(NA,27),rollapply(fulldata$TAVG,width=28,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_1M = c(rep(NA,29),rollapply(fulldata$TAVG,width=30,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_2M = c(rep(NA,59),rollapply(fulldata$TAVG,width=60,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_3M = c(rep(NA,89),rollapply(fulldata$TAVG,width=90,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_4M = c(rep(NA,119),rollapply(fulldata$TAVG,width=120,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_5M = c(rep(NA,149),rollapply(fulldata$TAVG,width=150,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_6M = c(rep(NA,179),rollapply(fulldata$TAVG,width=180,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_12M = c(rep(NA,364),rollapply(fulldata$TAVG,width=365,FUN=mean,na.rm=TRUE))
  fulldata$TAVG_24M = c(rep(NA,729),rollapply(fulldata$TAVG,width=730,FUN=mean,na.rm=TRUE))
  
  fulldata$PRCPAVG_3D = c(rep(NA,2),rollapply(fulldata$PRCP,width=3,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_1W = c(rep(NA,6),rollapply(fulldata$PRCP,width=7,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_2W = c(rep(NA,13),rollapply(fulldata$PRCP,width=14,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_3W = c(rep(NA,20),rollapply(fulldata$PRCP,width=21,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_4W = c(rep(NA,27),rollapply(fulldata$PRCP,width=28,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_1M = c(rep(NA,29),rollapply(fulldata$PRCP,width=30,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_2M = c(rep(NA,59),rollapply(fulldata$PRCP,width=60,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_3M = c(rep(NA,89),rollapply(fulldata$PRCP,width=90,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_4M = c(rep(NA,119),rollapply(fulldata$PRCP,width=120,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_5M = c(rep(NA,149),rollapply(fulldata$PRCP,width=150,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_6M = c(rep(NA,179),rollapply(fulldata$PRCP,width=180,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_12M = c(rep(NA,364),rollapply(fulldata$PRCP,width=365,FUN=mean,na.rm=TRUE))
  fulldata$PRCPAVG_24M = c(rep(NA,729),rollapply(fulldata$PRCP,width=730,FUN=mean,na.rm=TRUE))
  
  fulldata$PRCP_3D = c(rep(NA,2),rollapply(fulldata$PRCP,width=3,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_1W = c(rep(NA,6),rollapply(fulldata$PRCP,width=7,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_2W = c(rep(NA,13),rollapply(fulldata$PRCP,width=14,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_3W = c(rep(NA,20),rollapply(fulldata$PRCP,width=21,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_4W = c(rep(NA,27),rollapply(fulldata$PRCP,width=28,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_1M = c(rep(NA,29),rollapply(fulldata$PRCP,width=30,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_2M = c(rep(NA,59),rollapply(fulldata$PRCP,width=60,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_3M = c(rep(NA,89),rollapply(fulldata$PRCP,width=90,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_4M = c(rep(NA,119),rollapply(fulldata$PRCP,width=120,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_5M = c(rep(NA,149),rollapply(fulldata$PRCP,width=150,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_6M = c(rep(NA,179),rollapply(fulldata$PRCP,width=180,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_12M = c(rep(NA,364),rollapply(fulldata$PRCP,width=365,FUN=sum,na.rm=TRUE))
  fulldata$PRCP_24M = c(rep(NA,729),rollapply(fulldata$PRCP,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$MAXPRCP1W_2W = c(rep(NA,13),rollapply(fulldata$PRCP_1W,width=14,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_3W = c(rep(NA,20),rollapply(fulldata$PRCP_1W,width=21,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_4W = c(rep(NA,27),rollapply(fulldata$PRCP_1W,width=28,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_1M = c(rep(NA,29),rollapply(fulldata$PRCP_1W,width=30,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_2M = c(rep(NA,59),rollapply(fulldata$PRCP_1W,width=60,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_3M = c(rep(NA,89),rollapply(fulldata$PRCP_1W,width=90,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_4M = c(rep(NA,119),rollapply(fulldata$PRCP_1W,width=120,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_5M = c(rep(NA,149),rollapply(fulldata$PRCP_1W,width=150,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_6M = c(rep(NA,179),rollapply(fulldata$PRCP_1W,width=180,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_12M = c(rep(NA,364),rollapply(fulldata$PRCP_1W,width=365,FUN=max,na.rm=TRUE))
  fulldata$MAXPRCP1W_24M = c(rep(NA,729),rollapply(fulldata$PRCP_1W,width=730,FUN=max,na.rm=TRUE))
  
  fulldata$MDRN = ifelse(fulldata$PRCP<0.01,0,1)
  fulldata$MDRN_1W = c(rep(NA,6),rollapply(fulldata$MDRN,width=7,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_2W = c(rep(NA,13),rollapply(fulldata$MDRN,width=14,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_3W = c(rep(NA,20),rollapply(fulldata$MDRN,width=21,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_4W = c(rep(NA,27),rollapply(fulldata$MDRN,width=28,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_1M = c(rep(NA,29),rollapply(fulldata$MDRN,width=30,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_2M = c(rep(NA,59),rollapply(fulldata$MDRN,width=60,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_3M = c(rep(NA,89),rollapply(fulldata$MDRN,width=90,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_4M = c(rep(NA,119),rollapply(fulldata$MDRN,width=120,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_5M = c(rep(NA,149),rollapply(fulldata$MDRN,width=150,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_6M = c(rep(NA,179),rollapply(fulldata$MDRN,width=180,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_12M = c(rep(NA,364),rollapply(fulldata$MDRN,width=365,FUN=sum,na.rm=TRUE))
  fulldata$MDRN_24M = c(rep(NA,729),rollapply(fulldata$MDRN,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$TMAXAVG_1W = c(rep(NA,6),rollapply(fulldata$TMAX,width=7,FUN=mean,na.rm=TRUE))
  
  fulldata$MAXTMAX1W_2W = c(rep(NA,13),rollapply(fulldata$TMAXAVG_1W,width=14,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_3W = c(rep(NA,20),rollapply(fulldata$TMAXAVG_1W,width=21,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_4W = c(rep(NA,27),rollapply(fulldata$TMAXAVG_1W,width=28,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_1M = c(rep(NA,29),rollapply(fulldata$TMAXAVG_1W,width=30,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_2M = c(rep(NA,59),rollapply(fulldata$TMAXAVG_1W,width=60,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_3M = c(rep(NA,89),rollapply(fulldata$TMAXAVG_1W,width=90,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_4M = c(rep(NA,119),rollapply(fulldata$TMAXAVG_1W,width=120,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_5M = c(rep(NA,149),rollapply(fulldata$TMAXAVG_1W,width=150,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_6M = c(rep(NA,179),rollapply(fulldata$TMAXAVG_1W,width=180,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_12M = c(rep(NA,364),rollapply(fulldata$TMAXAVG_1W,width=365,FUN=max,na.rm=TRUE))
  fulldata$MAXTMAX1W_24M = c(rep(NA,729),rollapply(fulldata$TMAXAVG_1W,width=730,FUN=max,na.rm=TRUE))
  
  fulldata$TMINAVG_1W = c(rep(NA,6),rollapply(fulldata$TMAX,width=7,FUN=mean,na.rm=TRUE))
  
  fulldata$MINTMIN1W_2W = c(rep(NA,13),rollapply(fulldata$TMINAVG_1W,width=14,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_3W = c(rep(NA,20),rollapply(fulldata$TMINAVG_1W,width=21,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_4W = c(rep(NA,27),rollapply(fulldata$TMINAVG_1W,width=28,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_1M = c(rep(NA,29),rollapply(fulldata$TMINAVG_1W,width=30,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_2M = c(rep(NA,59),rollapply(fulldata$TMINAVG_1W,width=60,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_3M = c(rep(NA,89),rollapply(fulldata$TMINAVG_1W,width=90,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_4M = c(rep(NA,119),rollapply(fulldata$TMINAVG_1W,width=120,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_5M = c(rep(NA,149),rollapply(fulldata$TMINAVG_1W,width=150,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_6M = c(rep(NA,179),rollapply(fulldata$TMINAVG_1W,width=180,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_12M = c(rep(NA,364),rollapply(fulldata$TMINAVG_1W,width=365,FUN=min,na.rm=TRUE))
  fulldata$MINTMIN1W_24M = c(rep(NA,729),rollapply(fulldata$TMINAVG_1W,width=730,FUN=min,na.rm=TRUE))
  
  fulldata$HOT1S_1W = c(rep(NA,6),rollapply(fulldata$HOT1S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_2W = c(rep(NA,13),rollapply(fulldata$HOT1S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_3W = c(rep(NA,20),rollapply(fulldata$HOT1S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_4W = c(rep(NA,27),rollapply(fulldata$HOT1S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_1M = c(rep(NA,29),rollapply(fulldata$HOT1S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_2M = c(rep(NA,59),rollapply(fulldata$HOT1S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_3M = c(rep(NA,89),rollapply(fulldata$HOT1S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_4M = c(rep(NA,119),rollapply(fulldata$HOT1S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_5M = c(rep(NA,149),rollapply(fulldata$HOT1S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_6M = c(rep(NA,179),rollapply(fulldata$HOT1S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_12M = c(rep(NA,364),rollapply(fulldata$HOT1S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT1S_24M = c(rep(NA,729),rollapply(fulldata$HOT1S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$HOT2S_1W = c(rep(NA,6),rollapply(fulldata$HOT2S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_2W = c(rep(NA,13),rollapply(fulldata$HOT2S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_3W = c(rep(NA,20),rollapply(fulldata$HOT2S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_4W = c(rep(NA,27),rollapply(fulldata$HOT2S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_1M = c(rep(NA,29),rollapply(fulldata$HOT2S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_2M = c(rep(NA,59),rollapply(fulldata$HOT2S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_3M = c(rep(NA,89),rollapply(fulldata$HOT2S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_4M = c(rep(NA,119),rollapply(fulldata$HOT2S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_5M = c(rep(NA,149),rollapply(fulldata$HOT2S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_6M = c(rep(NA,179),rollapply(fulldata$HOT2S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_12M = c(rep(NA,364),rollapply(fulldata$HOT2S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT2S_24M = c(rep(NA,729),rollapply(fulldata$HOT2S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$HOT3S_1W = c(rep(NA,6),rollapply(fulldata$HOT3S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_2W = c(rep(NA,13),rollapply(fulldata$HOT3S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_3W = c(rep(NA,20),rollapply(fulldata$HOT3S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_4W = c(rep(NA,27),rollapply(fulldata$HOT3S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_1M = c(rep(NA,29),rollapply(fulldata$HOT3S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_2M = c(rep(NA,59),rollapply(fulldata$HOT3S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_3M = c(rep(NA,89),rollapply(fulldata$HOT3S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_4M = c(rep(NA,119),rollapply(fulldata$HOT3S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_5M = c(rep(NA,149),rollapply(fulldata$HOT3S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_6M = c(rep(NA,179),rollapply(fulldata$HOT3S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_12M = c(rep(NA,364),rollapply(fulldata$HOT3S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT3S_24M = c(rep(NA,729),rollapply(fulldata$HOT3S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$HOT90_1W = c(rep(NA,6),rollapply(fulldata$HOT90,width=7,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_2W = c(rep(NA,13),rollapply(fulldata$HOT90,width=14,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_3W = c(rep(NA,20),rollapply(fulldata$HOT90,width=21,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_4W = c(rep(NA,27),rollapply(fulldata$HOT90,width=28,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_1M = c(rep(NA,29),rollapply(fulldata$HOT90,width=30,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_2M = c(rep(NA,59),rollapply(fulldata$HOT90,width=60,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_3M = c(rep(NA,89),rollapply(fulldata$HOT90,width=90,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_4M = c(rep(NA,119),rollapply(fulldata$HOT90,width=120,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_5M = c(rep(NA,149),rollapply(fulldata$HOT90,width=150,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_6M = c(rep(NA,179),rollapply(fulldata$HOT90,width=180,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_12M = c(rep(NA,364),rollapply(fulldata$HOT90,width=365,FUN=sum,na.rm=TRUE))
  fulldata$HOT90_24M = c(rep(NA,729),rollapply(fulldata$HOT90,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD1S_1W = c(rep(NA,6),rollapply(fulldata$COLD1S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_2W = c(rep(NA,13),rollapply(fulldata$COLD1S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_3W = c(rep(NA,20),rollapply(fulldata$COLD1S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_4W = c(rep(NA,27),rollapply(fulldata$COLD1S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_1M = c(rep(NA,29),rollapply(fulldata$COLD1S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_2M = c(rep(NA,59),rollapply(fulldata$COLD1S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_3M = c(rep(NA,89),rollapply(fulldata$COLD1S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_4M = c(rep(NA,119),rollapply(fulldata$COLD1S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_5M = c(rep(NA,149),rollapply(fulldata$COLD1S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_6M = c(rep(NA,179),rollapply(fulldata$COLD1S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_12M = c(rep(NA,364),rollapply(fulldata$COLD1S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD1S_24M = c(rep(NA,729),rollapply(fulldata$COLD1S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD2S_1W = c(rep(NA,6),rollapply(fulldata$COLD2S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_2W = c(rep(NA,13),rollapply(fulldata$COLD2S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_3W = c(rep(NA,20),rollapply(fulldata$COLD2S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_4W = c(rep(NA,27),rollapply(fulldata$COLD2S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_1M = c(rep(NA,29),rollapply(fulldata$COLD2S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_2M = c(rep(NA,59),rollapply(fulldata$COLD2S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_3M = c(rep(NA,89),rollapply(fulldata$COLD2S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_4M = c(rep(NA,119),rollapply(fulldata$COLD2S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_5M = c(rep(NA,149),rollapply(fulldata$COLD2S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_6M = c(rep(NA,179),rollapply(fulldata$COLD2S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_12M = c(rep(NA,364),rollapply(fulldata$COLD2S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD2S_24M = c(rep(NA,729),rollapply(fulldata$COLD2S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD3S_1W = c(rep(NA,6),rollapply(fulldata$COLD3S,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_2W = c(rep(NA,13),rollapply(fulldata$COLD3S,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_3W = c(rep(NA,20),rollapply(fulldata$COLD3S,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_4W = c(rep(NA,27),rollapply(fulldata$COLD3S,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_1M = c(rep(NA,29),rollapply(fulldata$COLD3S,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_2M = c(rep(NA,59),rollapply(fulldata$COLD3S,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_3M = c(rep(NA,89),rollapply(fulldata$COLD3S,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_4M = c(rep(NA,119),rollapply(fulldata$COLD3S,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_5M = c(rep(NA,149),rollapply(fulldata$COLD3S,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_6M = c(rep(NA,179),rollapply(fulldata$COLD3S,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_12M = c(rep(NA,364),rollapply(fulldata$COLD3S,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD3S_24M = c(rep(NA,729),rollapply(fulldata$COLD3S,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$COLD32_1W = c(rep(NA,6),rollapply(fulldata$COLD32,width=7,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_2W = c(rep(NA,13),rollapply(fulldata$COLD32,width=14,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_3W = c(rep(NA,20),rollapply(fulldata$COLD32,width=21,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_4W = c(rep(NA,27),rollapply(fulldata$COLD32,width=28,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_1M = c(rep(NA,29),rollapply(fulldata$COLD32,width=30,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_2M = c(rep(NA,59),rollapply(fulldata$COLD32,width=60,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_3M = c(rep(NA,89),rollapply(fulldata$COLD32,width=90,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_4M = c(rep(NA,119),rollapply(fulldata$COLD32,width=120,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_5M = c(rep(NA,149),rollapply(fulldata$COLD32,width=150,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_6M = c(rep(NA,179),rollapply(fulldata$COLD32,width=180,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_12M = c(rep(NA,364),rollapply(fulldata$COLD32,width=365,FUN=sum,na.rm=TRUE))
  fulldata$COLD32_24M = c(rep(NA,729),rollapply(fulldata$COLD32,width=730,FUN=sum,na.rm=TRUE))
  
  fulldata$YEAR = as.numeric(substr(fulldata$DATE,1,4))
  years = unique(fulldata$YEAR)
  
  fulldata2 = subset(fulldata,YEAR>=futureperiod[1] & YEAR<=futureperiod[2])
  
  #####
  
  load(paste("/home/woot0002/streamflowtests_3daymovingmean_",gaugenum,".Rdata",sep=""))
  #load(paste("/home/woot0002/streamflowtests_3daymovingmean_08146000.Rdata",sep=""))
  corsmaxidx = c()
  corsmaxval = c()
  for(I in 1:length(corsusedall)){
    
    tmp = do.call("c",corsusedall[[I]])
    corsmaxidx[I] = which(tmp==max(tmp))
    corsmaxval[I] = max(tmp)
    
  }
  
  maxidx = which(corsmaxval==max(corsmaxval))
  varnamesin = varsusedall[[maxidx]][[corsmaxidx[maxidx]]]
  
  allnames = names(fulldata2)
  usevaridx = which(allnames %in% varnamesin)
  
  #####
  # All vars
  
  xmat = as.matrix(fulldata2[,usevaridx])
  
  load(paste("/home/woot0002/",projfilebreakdown$obs[i],"streamflowfit_",gaugenum,".Rdata",sep=""))
  #load(paste("/home/woot0002/",projfilebreakdown$obs[i],"streamflowfit_08146000.Rdata",sep=""))
  best_step <- fit$df[which.min(fit$RSS)]
  predictions <- predict(fit, xmat, s=best_step, type="fit")$fit
  predstreamflow = ifelse(predictions<0,0,predictions)
  #plot(density(predstreamflow))
  
  projresults[[i]] = list(predstreamflow)
  message("Finished gathering results for file ",i," / ",length(prfile_proj))
}

################
# Save results

save(list=c("histfilebreakdown","projfilebreakdown","histresults","projresults"),file=paste("/home/woot0002/CPREP_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))
} else {
  load(file=paste("/home/woot0002/CPREP_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))
}

################
# projected changes and comparison with obs

streamflow = read.table(paste("/home/woot0002/streamflow_",gaugenum,sep=""),header=TRUE,sep="\t",fill=TRUE)

if(gaugenum!="08151500" & gaugenum!="08144500" & gaugenum!="08148500"){
  names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")
  streamidx = 4
} else {
  names(streamflow) = c("agency","site","DATE","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC","streamflow_mean","streamflow_mean_QC")
  streamidx = 8
}


streamflow$DATE = as.character(streamflow$DATE)
streamflow$streamflow_mean=as.numeric(streamflow$streamflow_mean)
tmp = c(NA,rollapply(streamflow$streamflow_mean,3,mean,na.rm=TRUE),NA)
streamflow$streamflow_mean=tmp

streamflow$YEAR = as.numeric(substr(streamflow$DATE,1,4))
streamflow = subset(streamflow,YEAR>=1981 & YEAR<=2005)

#plot(density(log(streamflow$streamflow_mean[-which(is.na(streamflow$streamflow_mean)==TRUE)])),lty=2,main="Obs vs. Historical DS streamflow 1981-2005",xlab="ln(streamflow)",ylim=c(0,1),xlim=c(0,10))
plot(density(log(streamflow$streamflow_mean)),lty=2,main="Obs vs. Historical DS streamflow 1981-2005",xlab="ln(streamflow)",ylim=c(0,1),xlim=c(0,10))

for(I in 1:length(prfile_hist)){
  tmpdat = do.call("c",histresults[[I]])
  idxout = which(is.na(tmpdat)==TRUE)
  if(length(idxout)>0){
    lines(density(log(tmpdat[-idxout])),col="gray",lwd=0.5)
  } else {
    lines(density(log(tmpdat)),col="gray",lwd=0.5)
  }
}

legend("topleft",legend=c("Observed","Modeled (N=27)"),col=c("black","gray"),lwd==c(1,0.5),lty=c(2,1))


#plot(density(log(streamflow$streamflow_mean[-which(is.na(streamflow$streamflow_mean)==TRUE)])),lty=2,main="Obs vs. Downscaled Streamflow",xlab="ln(streamflow)",xlim=c(-5,10),ylim=c(0,1))
plot(density(log(streamflow$streamflow_mean)),lty=2,main="Obs vs. Downscaled Streamflow",xlab="ln(streamflow)",xlim=c(0,10),ylim=c(0,1))

for(I in 1:length(prfile_hist)){
  tmpdat = do.call("c",histresults[[I]])
  idxout = which(is.na(tmpdat)==TRUE)
  if(length(idxout)>0){
    lines(density(log(tmpdat[-idxout])),col="gray",lwd=0.5)
  } else {
    lines(density(log(tmpdat)),col="gray",lwd=0.5)
  }
}

for(I in 1:length(prfile_proj)){
  tmpdat = do.call("c",projresults[[I]])
  idxout = which(is.na(tmpdat)==TRUE)
  if(projfilebreakdown$scen[I]=="rcp26"){
    colval = "green"
  }
  if(projfilebreakdown$scen[I]=="rcp45"){
    colval = "blue"
  }
  if(projfilebreakdown$scen[I]=="rcp85"){
    colval = "red"
  }
  if(length(idxout)>0){
    lines(density(log(tmpdat[-idxout])),col=colval,lwd=0.5)
  } else {
    lines(density(log(tmpdat)),col=colval,lwd=0.5)
  }
}


legend("topleft",legend=c("Observed (1981-2005)","Historical DS (1981-2005)",paste("RCP 2.6 DS (",futureperiod[1],"-",futureperiod[2],")",sep=""),paste("RCP 4.5 DS (",futureperiod[1],"-",futureperiod[2],")",sep=""),paste("RCP 8.5 DS (",futureperiod[1],"-",futureperiod[2],")",sep="")),col=c("black","gray","green","blue","red"),lwd==c(1,0.5,0.5,0.5,0.5),lty=c(2,1,1,1,1))

###########
# boxplots

# combined historical daily obs and modeled

namesidx = which(names(streamflow)=="DATE" | names(streamflow)=="streamflow_mean")
streamflow_hist = streamflow[,namesidx]

for(i in 1:nrow(histfilebreakdown)){
  tmp = as.vector(histresults[[i]])[[1]]
  if(length(tmp)<nrow(streamflow_hist)){
    tmp2 = rep(NA,nrow(streamflow_hist))
    idxin = which(substr(streamflow_hist$DATE,6,10)!="02-29")
    tmp2[idxin]=tmp
    streamflow_hist = cbind(streamflow_hist,tmp2)
  } else {
    streamflow_hist = cbind(streamflow_hist,tmp)
  }
  names(streamflow_hist)[length(streamflow_hist)]=paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")  
}


projdates = seq(as.Date(paste(futureperiod[1],"-01-01",sep="")),as.Date(paste(futureperiod[2],"-12-31",sep="")),by="day")

for(i in 1:nrow(projfilebreakdown)){
  tmp = as.vector(projresults[[i]])[[1]]
  
  if(length(tmp)<length(projdates)){
    tmp2 = rep(NA,length(projdates))
    idxin = which(substr(projdates,6,10)!="02-29")
    tmp2[idxin]=tmp
    
    if(i==1){
      streamflow_proj = data.frame(projdates,tmp2)
      names(streamflow_proj) = c("DATE",paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_") )
    } else {
      streamflow_proj = cbind(streamflow_proj,tmp2)
      names(streamflow_proj)[length(streamflow_proj)]=paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_")
    }
    
  } else {
    
    if(i==1){
      streamflow_proj = data.frame(projdates,tmp)
      names(streamflow_proj) = c("DATE",paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_") )
    } else {
      streamflow_proj = cbind(streamflow_proj,tmp)
      names(streamflow_proj)[length(streamflow_proj)]=paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_")
    }
  }
  
}

# anomaly calcs daily

streamflow_anom = streamflow_proj
for(i in 2:ncol(streamflow_proj)){
  nametmp = names(streamflow_proj)[i]
  histidx = which(names(streamflow_hist)==substr(nametmp,7,nchar(nametmp)))
  
  streamflow_anom[,i] = streamflow_proj[,i]-mean(streamflow_hist[,histidx],na.rm=TRUE)
}

###
# seasonal calcs

streamflow_hist$season = NA
DJFidx = which(as.numeric(substr(streamflow_hist$DATE,6,7))==12 | as.numeric(substr(streamflow_hist$DATE,6,7))<=2)
MAMidx = which(as.numeric(substr(streamflow_hist$DATE,6,7))>=3 & as.numeric(substr(streamflow_hist$DATE,6,7))<=5)
JJAidx = which(as.numeric(substr(streamflow_hist$DATE,6,7))>=6 & as.numeric(substr(streamflow_hist$DATE,6,7))<=8)
SONidx = which(as.numeric(substr(streamflow_hist$DATE,6,7))>=9 & as.numeric(substr(streamflow_hist$DATE,6,7))<=11)
streamflow_hist$season[DJFidx] = 1
streamflow_hist$season[MAMidx] = 2
streamflow_hist$season[JJAidx] = 3
streamflow_hist$season[SONidx] = 4
streamflow_hist$year = as.numeric(substr(streamflow_hist$DATE,1,4))
streamflow_hist$seayear = paste(streamflow_hist$year,streamflow_hist$season,sep="-")
streamflow_hist_season =  aggregate(streamflow_hist[,2:29],by=list(seayear=streamflow_hist$seayear),mean,na.rm=TRUE)

streamflow_proj$season = NA
DJFidx = which(as.numeric(substr(streamflow_proj$DATE,6,7))==12 | as.numeric(substr(streamflow_proj$DATE,6,7))<=2)
MAMidx = which(as.numeric(substr(streamflow_proj$DATE,6,7))>=3 & as.numeric(substr(streamflow_proj$DATE,6,7))<=5)
JJAidx = which(as.numeric(substr(streamflow_proj$DATE,6,7))>=6 & as.numeric(substr(streamflow_proj$DATE,6,7))<=8)
SONidx = which(as.numeric(substr(streamflow_proj$DATE,6,7))>=9 & as.numeric(substr(streamflow_proj$DATE,6,7))<=11)
streamflow_proj$season[DJFidx] = 1
streamflow_proj$season[MAMidx] = 2
streamflow_proj$season[JJAidx] = 3
streamflow_proj$season[SONidx] = 4
streamflow_proj$year = as.numeric(substr(streamflow_proj$DATE,1,4))
streamflow_proj$seayear = paste(streamflow_proj$year,streamflow_proj$season,sep="-")
streamflow_proj_season =  aggregate(streamflow_proj[,2:82],by=list(seayear=streamflow_proj$seayear),mean,na.rm=TRUE)

# anomaly calcs seasonal
streamflow_anom_season = streamflow_proj_season
for(i in 2:ncol(streamflow_proj_season)){
  nametmp = names(streamflow_proj_season)[i]
  histidx = which(names(streamflow_hist_season)==substr(nametmp,7,nchar(nametmp)))
  
  for(s in 1:4){
    
    seasonidx_proj = which(as.numeric(substr(streamflow_proj_season$seayear,6,6))==s)
    seasonidx_hist = which(as.numeric(substr(streamflow_hist_season$seayear,6,6))==s)
    
    streamflow_anom_season[seasonidx_proj,i] = streamflow_proj_season[seasonidx_proj,i]-mean(streamflow_hist_season[seasonidx_hist,histidx],na.rm=TRUE)
  }
}

###
# yearly calcs

streamflow_hist_year =  aggregate(streamflow_hist[,2:29],by=list(year=streamflow_hist$year),mean,na.rm=TRUE)
streamflow_proj_year =  aggregate(streamflow_proj[,2:82],by=list(year=streamflow_proj$year),mean,na.rm=TRUE)

streamflow_anom_year = streamflow_proj_year
for(i in 2:ncol(streamflow_proj_year)){
  nametmp = names(streamflow_proj_year)[i]
  histidx = which(names(streamflow_hist_year)==substr(nametmp,7,nchar(nametmp)))
  
  streamflow_anom_year[,i] = streamflow_proj_year[,i]-mean(streamflow_hist_year[,histidx],na.rm=TRUE)
}

###
# yearly climo calcs

streamflow_hist_yearclimo =  apply(streamflow_hist_year[,2:29],2,mean,na.rm=TRUE)
streamflow_proj_yearclimo =  apply(streamflow_proj_year[,2:82],2,mean,na.rm=TRUE)

streamflow_anom_yearclimo = streamflow_proj_yearclimo
for(i in 1:length(streamflow_proj_yearclimo)){
  nametmp = names(streamflow_proj_yearclimo)[i]
  histidx = which(names(streamflow_hist_yearclimo)==substr(nametmp,7,nchar(nametmp)))
  
  streamflow_anom_yearclimo[i] = streamflow_proj_yearclimo[i]-streamflow_hist_yearclimo[histidx]
}


####
# rearrange tables

# yearly climo anomaly

sanom_climo = NULL

for(i in 1:81){
  
  vnamesplit = strsplit(names(streamflow_anom_yearclimo)[i],"_",fixed=TRUE)[[1]]
  scen = vnamesplit[1]
  GCM = vnamesplit[2]
  DS = vnamesplit[3]
  obs = vnamesplit[4]
  
  sanom = streamflow_anom_yearclimo[i]
  names(sanom)="sanom"
  
  tmp = data.frame(scen,GCM,DS,obs,sanom)
  sanom_climo = rbind(sanom_climo,tmp)
  
}

####
#  historical

shist_daily = shist_seasonal = shist_yearly =NULL

for(i in 2:29){
  
  #daily
  tmp = streamflow_hist[,c(1,i)]
  if(i==2){
    tmp$sim = "obs"
    tmp$group = "obs"
  } else {
    tmp$sim = names(streamflow_hist)[i]
    tmp$group = "model"
    names(tmp) = names(shist_daily)
  }
  shist_daily=rbind(shist_daily,tmp)
  
  #seasonal
  tmp = streamflow_hist_season[,c(1,i)]
  if(i==2){
    tmp$sim = "obs"
    tmp$group = "obs"
  } else {
    tmp$sim = names(streamflow_hist_season)[i]
    tmp$group = "model"
    names(tmp) = names(shist_seasonal)
  }
  shist_seasonal=rbind(shist_seasonal,tmp)
  
  #yearly
  tmp = streamflow_hist_year[,c(1,i)]
  if(i==2){
    tmp$sim = "obs"
    tmp$group = "obs"
  } else {
    tmp$sim = names(streamflow_hist_year)[i]
    tmp$group = "model"
    names(tmp) = names(shist_yearly)
  }
  shist_yearly=rbind(shist_yearly,tmp)
  
}

### projected and anoms

sproj_daily = sproj_seasonal = sproj_yearly = sanom_daily = sanom_seasonal = sanom_yearly =NULL

for(i in 2:82){
  
  #daily
  tmp = streamflow_proj[,c(1,i)]
  tmp$sim = names(streamflow_proj)[i]
  tmp$group = "model"
  names(tmp) = names(shist_daily)
  sproj_daily=rbind(sproj_daily,tmp)
  
  tmp = streamflow_anom[,c(1,i)]
  tmp$sim = names(streamflow_anom)[i]
  tmp$group = "model"
  names(tmp) = names(shist_daily)
  sanom_daily=rbind(sanom_daily,tmp)
  
  #seasonal
  tmp = streamflow_proj_season[,c(1,i)]
  tmp$sim = names(streamflow_proj_season)[i]
  tmp$group = "model"
  names(tmp) = names(shist_seasonal)
  sproj_seasonal=rbind(sproj_seasonal,tmp)
  
  tmp = streamflow_anom_season[,c(1,i)]
  tmp$sim = names(streamflow_anom_season)[i]
  tmp$group = "model"
  names(tmp) = names(shist_seasonal)
  sanom_seasonal=rbind(sanom_seasonal,tmp)
  
  
  #yearly
  tmp = streamflow_proj_year[,c(1,i)]
  tmp$sim = names(streamflow_proj_year)[i]
  tmp$group = "model"
  names(tmp) = names(shist_yearly)
  sproj_yearly=rbind(sproj_yearly,tmp)
  
  tmp = streamflow_anom_year[,c(1,i)]
  tmp$sim = names(streamflow_anom_year)[i]
  tmp$group = "model"
  names(tmp) = names(shist_yearly)
  sanom_yearly=rbind(sanom_yearly,tmp)
  
}

####
# cleanup

shist_seasonal$season = as.numeric(substr(shist_seasonal$seayear,6,6))
sproj_seasonal$season = as.numeric(substr(sproj_seasonal$seayear,6,6))
sanom_seasonal$season = as.numeric(substr(sanom_seasonal$seayear,6,6))

shist_seasonal$season[which(shist_seasonal$season==1)]="DJF"
shist_seasonal$season[which(shist_seasonal$season==2)]="MAM"
shist_seasonal$season[which(shist_seasonal$season==3)]="JJA"
shist_seasonal$season[which(shist_seasonal$season==4)]="SON"
shist_seasonal$season <- factor(shist_seasonal$season,levels = c('DJF','MAM','JJA','SON'),ordered = TRUE)

sproj_seasonal$season[which(sproj_seasonal$season==1)]="DJF"
sproj_seasonal$season[which(sproj_seasonal$season==2)]="MAM"
sproj_seasonal$season[which(sproj_seasonal$season==3)]="JJA"
sproj_seasonal$season[which(sproj_seasonal$season==4)]="SON"
sproj_seasonal$season <- factor(sproj_seasonal$season,levels = c('DJF','MAM','JJA','SON'),ordered = TRUE)

sanom_seasonal$season[which(sanom_seasonal$season==1)]="DJF"
sanom_seasonal$season[which(sanom_seasonal$season==2)]="MAM"
sanom_seasonal$season[which(sanom_seasonal$season==3)]="JJA"
sanom_seasonal$season[which(sanom_seasonal$season==4)]="SON"
sanom_seasonal$season <- factor(sanom_seasonal$season,levels = c('DJF','MAM','JJA','SON'),ordered = TRUE)



sproj_yearly$scen = substr(sproj_yearly$sim,1,5)
sanom_yearly$scen = substr(sanom_yearly$sim,1,5)
sproj_seasonal$scen = substr(sproj_seasonal$sim,1,5)
sanom_seasonal$scen = substr(sanom_seasonal$sim,1,5)
sproj_daily$scen = substr(sproj_daily$sim,1,5)
sanom_daily$scen = substr(sanom_daily$sim,1,5)




aggregate(sanom_climo$sanom,by=list(scen=sanom_climo$scen),mean,na.rm=TRUE)
aggregate(sanom_climo$sanom,by=list(scen=sanom_climo$scen),sd,na.rm=TRUE)

aggregate(sanom_yearly$streamflow_mean,by=list(scen=sanom_yearly$scen),mean,na.rm=TRUE)
aggregate(sanom_yearly$streamflow_mean,by=list(scen=sanom_yearly$scen),sd,na.rm=TRUE)


####
# boxplots

library(ggplot2)

ggplot(shist_yearly, aes(x=group, y=streamflow_mean,fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Yearly Average Streamflow (1981-2005)")+xlab("")+ylab("Average Streamflow (ft^3/s)") 

ggplot(shist_seasonal, aes(x=group, y=streamflow_mean,fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Seasonal Average Streamflow (1981-2005)")+xlab("")+ylab("Average Streamflow (ft^3/s)")+facet_wrap(vars(season))


ggplot(shist_daily, aes(x=group, y=streamflow_mean,fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Daily Average Streamflow (1981-2005)")+xlab("")+ylab("Average Streamflow (ft^3/s)") 



ggplot(sanom_yearly, aes(x=scen, y=streamflow_mean,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed")+ scale_fill_manual(values=c( "green", "blue","red")) +
  ggtitle(paste("Projected Change - yearly average \nstreamflow (",futureperiod[1],"-",futureperiod[2],")",sep=""))+xlab("")+ylab("Change in Mean Streamflow (ft^3/s)") 
ggsave(paste("/home/woot0002/streamflow_projchangeyearly_",gaugenum,".pdf",sep=""),device = "pdf",width=5,height=5)

ggplot(sanom_seasonal, aes(x=scen, y=streamflow_mean,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed")+ scale_fill_manual(values=c( "green", "blue","red")) +
  ggtitle(paste("Projected Change - seasonal average \nstreamflow (",futureperiod[1],"-",futureperiod[2],")",sep=""))+xlab("")+ylab("Change in Mean Streamflow (ft^3/s)") +facet_wrap(vars(season))

ggplot(sanom_daily, aes(x=scen, y=streamflow_mean,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed")+ scale_fill_manual(values=c( "green", "blue","red")) +
  ggtitle("Projected Change - daily average streamflow (2070-2099)")+xlab("")+ylab("Change in Mean Streamflow (ft^3/s)")

####
ggplot(sanom_climo, aes(x=scen, y=sanom,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed")+ scale_fill_manual(values=c( "green", "blue","red")) +
  ggtitle(paste("Projected Change - yearly average \nstreamflow (",futureperiod[1],"-",futureperiod[2],")",sep=""))+xlab("")+ylab("Change in Mean Streamflow (ft^3/s)") 
ggsave(paste("/home/woot0002/streamflow_projchangeyearlyclimo_",gaugenum,".pdf",sep=""),device = "pdf",width=5,height=5)



