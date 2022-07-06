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
library(mailR)

prfile_hist = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr*historical*.nc",intern=TRUE))

tasmaxfile_hist = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/tasmax/EDQM/tasmax*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/tasmax/PARM/tasmax*historical*.nc",intern=TRUE))

tasminfile_hist = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin*historical*.nc",intern=TRUE))

prfile = "/data2/3to5/I35/pr/PRISM/pr_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
tasmaxfile = "/data2/3to5/I35/tasmax/PRISM/tasmax_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
tasminfile = "/data2/3to5/I35/tasmin/PRISM/tasmin_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

streamflowID = "08144600"
loclon = -99.26885986 # gauging station location 
loclat = 30.91913605

streamflow = read.table(paste("/home/woot0002/streamflow_",streamflowID,sep=""),header=TRUE,sep="\t",fill=TRUE) # watch it with reordered tables based on gauge ID

if(streamflowID!="08144500"){
  names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")
  streamidx = 4
} else {
  names(streamflow) = c("agency","site","DATE","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC","streamflow_mean","streamflow_mean_QC")
  streamidx = 8
}


####
# gridded obs data - find data point

nctest = nc_open(prfile)
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

###
# Get PRCP, TMAX, and TMIN data

nctest = nc_open(prfile)
PRCP = ncvar_get(nctest,"pr",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))*86400
nc_close(nctest)

nctest = nc_open(tasmaxfile)
TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
nc_close(nctest)

nctest = nc_open(tasminfile)
TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
nc_close(nctest)

#####
# unit conversion

PRCP = PRCP/25.4
TMAX = TMAX-273.15
TMAX = (TMAX*9/5)+32
TMIN = TMIN-273.15
TMIN = (TMIN*9/5)+32
DATE = as.character(seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"))
climate = data.frame(DATE,PRCP,TMAX,TMIN)

#####
# combine data

streamflow$DATE = as.character(streamflow$DATE)
fulldata = merge(streamflow,climate,by="DATE")

fulldata$streamflow_mean=as.numeric(fulldata$streamflow_mean)
tmp = c(NA,rollapply(fulldata$streamflow_mean,3,mean,na.rm=TRUE),NA)
fulldata$streamflow_mean=tmp
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

percentile_bottom = quantile(fulldata$PRCP_1M,probs=0.25,na.rm=TRUE)
percentile_top = quantile(fulldata$PRCP_1M,probs=0.75,na.rm=TRUE)

precip_hardcode = fulldata$PRCP_1M
precip_hardcode=ifelse(fulldata$PRCP_1M<=percentile_bottom | fulldata$PRCP_1M>=percentile_top,ifelse(fulldata$PRCP_1M<=percentile_bottom,0,2),1)

tmp = precip_hardcode[1:(length(fulldata$PRCP_1M)-10)]-precip_hardcode[11:length(fulldata$PRCP_1M)]
tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,10))
tmppd = c(ifelse(tmp==2,1,0),rep(NA,10))
PDcount_PRCP = sum(tmppd,na.rm=TRUE)
DPcount_PRCP = sum(tmpdp,na.rm=TRUE)

fulldata$PRCP_DP = tmpdp
fulldata$PRCP_PD = tmppd
fulldata$PRCP_WP = c(tmp,rep(NA,10))


fulldata$SFM_1M = c(rep(NA,29),rollapply(fulldata$streamflow_mean,width=30,FUN=mean,na.rm=TRUE))
#fulldata$SFM_1M = c(rep(NA,29),rollapply(fulldata$streamflow_mean,width=30,FUN=sum,na.rm=TRUE))

percentile_bottom = quantile(fulldata$SFM_1M,probs=0.25,na.rm=TRUE)
percentile_top = quantile(fulldata$SFM_1M,probs=0.75,na.rm=TRUE)

precip_hardcode = fulldata$SFM_1M
precip_hardcode=ifelse(fulldata$SFM_1M<=percentile_bottom | fulldata$SFM_1M>=percentile_top,ifelse(fulldata$SFM_1M<=percentile_bottom,0,2),1)

tmp = precip_hardcode[1:(length(fulldata$SFM_1M)-10)]-precip_hardcode[11:length(fulldata$SFM_1M)]
tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,10))
tmppd = c(ifelse(tmp==2,1,0),rep(NA,10))
PDcount_SFM = sum(tmppd,na.rm=TRUE)
DPcount_SFM = sum(tmpdp,na.rm=TRUE)

fulldata$SFM_DP = tmpdp
fulldata$SFM_PD = tmppd
fulldata$SFM_WP = c(tmp,rep(NA,10))

#####
# correlation table

fulldata$YEAR = as.numeric(substr(fulldata$DATE,1,4))
years = unique(fulldata$YEAR)

fulldata2 = fulldata[which(fulldata$YEAR>=1950 & fulldata$YEAR<=2015),]
years2 = unique(fulldata2$YEAR)

startidx = which(names(fulldata2)=="PRCP")

cortable = NULL

for(j in startidx:ncol(fulldata2)){
  J=j
  varname =  names(fulldata2)[j]
  corval = cor(as.numeric(fulldata2$streamflow_mean),as.numeric(fulldata2[,j]),method="spearman",use="pairwise.complete.obs")
  tmp=data.frame(J,varname,corval)
  cortable=rbind(cortable,tmp)
}
cortable$abscor = abs(cortable$corval)
sortidx = order(cortable$abscor,decreasing=TRUE)
cortable= cortable[sortidx,]

TMAXGroup1 = c("TAVG","TMAX","TMIN","TOBS","TAVG_3D","HOT1S","HOT2S","HOT3S","HOT90","COLD1S","COLD2S","COLD3S","COLD32")
PRCPGroup1 = c("SNOW","SNWD","PRCP","PRCP_3D","MDRN","PRCPAVG_3D")

TMAXGroup2 = c("TAVG_1W","TAVG_2W","TAVG_3W","TAVG_4W","TMAXAVG_1W","MAXTMAX1W_2W","MAXTMAX1W_3W","MAXTMAX1W_4W","TMINAVG_1W","MINTMIN1W_2W","MINTMIN1W_3W","MINTMIN1W_4W"
               ,"HOT1S_1W","HOT1S_2W","HOT1S_3W","HOT1S_4W","HOT2S_1W","HOT2S_2W","HOT2S_3W","HOT2S_4W","HOT3S_1W","HOT3S_2W","HOT3S_3W","HOT3S_4W"
               ,"HOT90_1W","HOT90_2W","HOT90_3W","HOT90_4W","COLD1S_1W","COLD1S_2W","COLD1S_3W","COLD1S_4W","COLD2S_1W","COLD2S_2W","COLD2S_3W","COLD2S_4W"
               ,"COLD3S_1W","COLD3S_2W","COLD3S_3W","COLD3S_4W","COLD32_1W","COLD32_2W","COLD32_3W","COLD32_4W")
PRCPGroup2 = c("PRCPAVG_1W","PRCPAVG_2W","PRCPAVG_3W","PRCPAVG_4W","PRCP_1W","PRCP_2W","PRCP_3W","PRCP_4W"
               ,"MDRN_1W","MDRN_2W","MDRN_3W","MDRN_4W","MAXPRCP1W_2W","MAXPRCP1W_3W","MAXPRCP1W_4W")

TMAXGroup3 = c("TAVG_1M","TAVG_2M","TAVG_3M","TAVG_4M","TAVG_5M","TAVG_6M","MAXTMAX1W_1M","MAXTMAX1W_2M","MAXTMAX1W_3M"
               ,"MAXTMAX1W_4M","MAXTMAX1W_5M","MAXTMAX1W_6M","MINTMIN1W_1M","MINTMIN1W_2M","MINTMIN1W_3M","MINTMIN1W_4M"
               ,"MINTMIN1W_5M","MINTMIN1W_6M","HOT1S_1M","HOT1S_2M","HOT1S_3M","HOT1S_4M","HOT1S_5M","HOT1S_6M"
               ,"HOT2S_1M","HOT2S_2M","HOT2S_3M","HOT2S_4M","HOT2S_5M","HOT2S_6M"
               ,"HOT3S_1M","HOT3S_2M","HOT3S_3M","HOT3S_4M","HOT3S_5M","HOT3S_6M"
               ,"HOT90_1M","HOT90_2M","HOT90_3M","HOT90_4M","HOT90_5M","HOT90_6M"
               ,"COLD1S_1M","COLD1S_2M","COLD1S_3M","COLD1S_4M","COLD1S_5M","COLD1S_6M"
               ,"COLD2S_1M","COLD2S_2M","COLD2S_3M","COLD2S_4M","COLD2S_5M","COLD2S_6M"
               ,"COLD3S_1M","COLD3S_2M","COLD3S_3M","COLD3S_4M","COLD3S_5M","COLD3S_6M"
               ,"COLD32_1M","COLD32_2M","COLD32_3M","COLD32_4M","COLD32_5M","COLD32_6M")
PRCPGroup3 = c("PRCPAVG_1M","PRCPAVG_2M","PRCPAVG_3M","PRCPAVG_4M","PRCPAVG_5M","PRCPAVG_6M"
               ,"PRCP_1M","PRCP_2M","PRCP_3M","PRCP_4M","PRCP_5M","PRCP_6M"
               ,"MAXPRCP1W_1M","MAXPRCP1W_2M","MAXPRCP1W_3M","MAXPRCP1W_4M","MAXPRCP1W_5M","MAXPRCP1W_6M"
               ,"MDRN_1M","MDRN_2M","MDRN_3M","MDRN_4M","MDRN_5M","MDRN_6M")

TMAXGroup4 = c("TAVG_12M","TAVG_24M","MAXTMAX1W_12M","MAXTMAX1W_24M","MINTMIN1W_12M","MINTMIN1W_24M"
               ,"HOT1S_12M","HOT1S_24M","HOT2S_12M","HOT2S_24M","HOT3S_12M","HOT3S_24M","HOT90_12M","HOT90_24M"
               ,"COLD1S_12M","COLD1S_24M","COLD2S_12M","COLD2S_24M","COLD2S_12M","COLD2S_24M","COLD32_12M","COLD32_24M")
PRCPGroup4 = c("PRCPAVG_12M","PRCPAVG_24M","PRCP_12M","PRCP_24M","MAXPRCP1W_12M","MAXPRCP1W_24M"
               ,"MDRN_12M","MDRN_24M")

corTMAXgroup1 = cortable[which(cortable$varname %in% TMAXGroup1),]
corPRCPgroup1 = cortable[which(cortable$varname %in% PRCPGroup1),]

corTMAXgroup2 = cortable[which(cortable$varname %in% TMAXGroup2),]
corPRCPgroup2 = cortable[which(cortable$varname %in% PRCPGroup2),]

corTMAXgroup3 = cortable[which(cortable$varname %in% TMAXGroup3),]
corPRCPgroup3 = cortable[which(cortable$varname %in% PRCPGroup3),]

corTMAXgroup4 = cortable[which(cortable$varname %in% TMAXGroup4),]
corPRCPgroup4 = cortable[which(cortable$varname %in% PRCPGroup4),]

namesused=c("PRCPAVG_3D","MDRN_4W","MAXPRCP1W_2W","HOT1S_1W","TMAXAVG_1W","HOT1S_3M","TAVG_3M","MDRN_3M","PRCPAVG_6M","PRCPAVG_12M","HOT90_12M","TAVG_12M")

usevaridx = c(11,35,39,61,88,85,91,29,45,69,46,71,119,34)
#names(traingroup)[usevaridx]
#####
# Make regression models to test

is.odd <- function(x) x %% 2 != 0

oddix = which(is.odd(years2)==TRUE)
evenix = which(is.odd(years2)==FALSE)

traingroup = NULL
predgroup = NULL
for(i in 1:length(oddix)){
  traingroup = rbind(traingroup,fulldata2[which(fulldata2$YEAR==years2[oddix[i]]),])
}

for(i in 1:length(evenix)){
  predgroup = rbind(predgroup,fulldata2[which(fulldata2$YEAR==years2[evenix[i]]),])
}

#####

#load(paste("/home/woot0002/streamflowtests_3daymovingmean_",streamflowID,".Rdata",sep=""))
load("/home/woot0002/streamflowtests_3daymovingmean_08146000.Rdata")
corsmaxidx = c()
corsmaxval = c()
for(i in 1:14){
  
  tmp = do.call("c",corsusedall[[i]])
  corsmaxidx[i] = which(tmp==max(tmp))
  corsmaxval[i] = max(tmp)
  
}

maxidx = which(corsmaxval==max(corsmaxval))
varnamesin = varsusedall[[maxidx]][[corsmaxidx[maxidx]]]

allnames = names(traingroup)
usevaridx = which(allnames %in% varnamesin)

#####
# All vars

xmat = as.matrix(traingroup[,usevaridx])
#xmat = xmat[,which(allnames %in% varnamesin)]
ymat = as.matrix(traingroup[,streamidx])

xmat2 = as.matrix(predgroup[,usevaridx])
#xmat2 = xmat2[,which(allnames %in% varnamesin)]
ymat2 = as.matrix(predgroup[,streamidx])

NAidx = which(is.na(ymat)==TRUE)
xmat = xmat[-NAidx,]
ymat = ymat[-NAidx,]

NAidx = which(is.na(xmat)==TRUE,arr.ind = TRUE)
NAidxu = unique(NAidx[,1])

xmat = xmat[-NAidxu,]
ymat = ymat[-NAidxu]

NAidx = which(is.na(ymat2)==TRUE)
xmat2 = xmat2[-NAidx,]
ymat2 = ymat2[-NAidx,]

NAidx = which(is.na(xmat2)==TRUE,arr.ind = TRUE)
NAidxu = unique(NAidx[,1])

if(length(NAidxu)>0){
  xmat2 = xmat2[-NAidxu,]
  ymat2 = ymat2[-NAidxu]
}

#fit <- lars(xmat,ymat , type="lasso")
#save("fit",file=paste("/home/woot0002/PRISMstreamflowfit_",streamflowID,".Rdata",sep=""))
#plot(fit)
load("/home/woot0002/PRISMstreamflowfit_08146000.Rdata")
best_step <- fit$df[which.min(fit$RSS)]
predictions <- predict(fit, xmat, s=best_step, type="fit")$fit
predictions = ifelse(predictions<0,0,predictions)

predictions2 <- predict(fit, xmat2, s=best_step, type="fit")$fit
predictions2 = ifelse(predictions2<0,0,predictions2)

RMSE_train = sqrt(mean((predictions-ymat)^2))
cor_train = cor(ymat,predictions)

RMSE_pred = sqrt(mean((predictions2-ymat2)^2))
cor_pred = cor(ymat2,predictions2)

logpred = log(predictions)
NAidx2=which(is.na(logpred)==TRUE | logpred == -Inf)

plot(density(log(ymat)),xlab="ln(streamflow)",main="Training period streamflow pdfs",ylim=c(0,1))
lines(density(logpred),col="blue")
legend("topright",legend=c("Observed","Modeled"),col=c("black","blue"),lwd=1)

logpred2 = log(predictions2)
NAidx2=which(is.na(logpred2)==TRUE | logpred2 == -Inf)

plot(density(log(ymat2)),xlab="ln(streamflow)",main="Prediction period streamflow pdfs",ylim=c(0,1))
lines(density(logpred2[-NAidx2]),col="blue")
legend("topright",legend=c("Observed","Modeled"),col=c("black","blue"),lwd=1)


plot(ymat,type="l")
lines(predictions,col="blue")
legend("topright",legend=c("Observed","Modeled"),col=c("black","blue"),lwd=1)

plot(ymat[400:600],type="l")
lines(predictions[400:600],col="blue")
legend("topright",legend=c("Observed","Modeled"),col=c("black","blue"),lwd=1)


plot(ymat2,type="l")
lines(predictions2,col="blue")
legend("topright",legend=c("Observed","Modeled"),col=c("black","blue"),lwd=1)

plot(ymat2[400:600],type="l")
lines(predictions2[400:600],col="blue")
legend("topright",legend=c("Observed","Modeled"),col=c("black","blue"),lwd=1)

######
######




