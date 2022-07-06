#####
#
# streamflow testing - regression

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

####

streamflowID = "08151500"
climateID = "2562869"

streamflow = read.table(paste("/home/woot0002/streamflow_",streamflowID,sep=""),header=TRUE,sep="\t",fill=TRUE) # watch it with reordered tables based on gauge ID
climate = read.table(paste("/home/woot0002/",climateID,".csv",sep=""),header=TRUE,sep=",")

if(streamflowID!="08151500" & streamflowID!="08144500" & streamflowID!="08148500"){
  names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")
  streamidx = 4
} else {
  names(streamflow) = c("agency","site","DATE","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC","streamflow_mean","streamflow_mean_QC")
  streamidx = 8
}

fulldata = merge(streamflow,climate,by="DATE")
fulldata$streamflow_mean=as.numeric(fulldata$streamflow_mean)

tmp = c(NA,rollapply(fulldata$streamflow_mean,3,mean,na.rm=TRUE),NA)

fulldata$streamflow_mean=tmp
fulldata$PRCP=as.numeric(fulldata$PRCP)
fulldata$TMAX=as.numeric(fulldata$TMAX)
fulldata$TMIN=as.numeric(fulldata$TMIN)
fulldata$TAVG=(as.numeric(fulldata$TMAX)+as.numeric(fulldata$TMIN))/2
fulldata$SNOW=as.numeric(fulldata$SNOW)
fulldata$DAPR=as.numeric(fulldata$DAPR)
fulldata$DASF=as.numeric(fulldata$DASF)
fulldata$MDPR=as.numeric(fulldata$MDPR)
fulldata$MDSF=as.numeric(fulldata$MDSF)
fulldata$SNWD=as.numeric(fulldata$SNWD)

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

fulldata2 = fulldata[which(fulldata$YEAR>=1981 & fulldata$YEAR<=2005),] # normal range 1950-2015
years2 = unique(fulldata2$YEAR)

#####
# regression information


load(paste("/home/woot0002/streamflowtests_3daymovingmean_",streamflowID,".Rdata",sep=""))
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

namesused = names(fulldata2)[usevaridx]

xmat = as.matrix(fulldata2[,usevaridx])
ymat = as.matrix(fulldata2[,streamidx])

####

NAidx1 = which(is.na(ymat)==TRUE)
if(length(NAidx1)>0){
  xmat = xmat[-NAidx1,]
  ymat = ymat[-NAidx1,]
}

NAidx2 = which(is.na(xmat)==TRUE,arr.ind = TRUE)
NAidxu = unique(NAidx2[,1])

if(length(NAidxu)>0){
  xmat = xmat[-NAidxu,]
  ymat = ymat[-NAidxu]
}

####

fit <- lars(xmat,ymat , type="lasso")
best_step <- fit$df[which.min(fit$RSS)]
predictions <- predict(fit, xmat, s=best_step, type="fit")$fit

fulldata2$streamflowmod = NA

idx =1:nrow(fulldata2)
idxuse = idx[-c(NAidx1,unique(NAidx2[,1]))]

fulldata2$streamflowmod[idxuse] = predictions

save(list=c("fulldata2"),file=paste("/home/woot0002/streamflowreg_3daymovingmean_",streamflowID,".Rdata",sep=""))




