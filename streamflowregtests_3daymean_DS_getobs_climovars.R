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

gaugenum = "08144500"
loclon = -99.78561401 # gauging station location 
loclat = 30.91913605
futureperiod = c(2006,2099)
runcalcs=TRUE
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
  
  histresults[[i]] = fulldata2
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
  
  projresults[[i]] = fulldata2
  message("Finished gathering results for file ",i," / ",length(prfile_proj))
}

################
# Save results

save(list=c("histfilebreakdown","projfilebreakdown","histresults","projresults"),file=paste("/home/woot0002/CPREP_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed_climovars.Rdata",sep=""))
} 