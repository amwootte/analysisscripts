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
library(mailR)

# email settings
emadd = "amwootte@ou.edu"
pswd = "xKFkcBZ6oN"


streamflow = read.table("/home/woot0002/streamflow_08146000",header=TRUE,sep="\t",fill=TRUE)
climate = read.table("/home/woot0002/2444385.csv",header=TRUE,sep=",")

names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")

fulldata = merge(streamflow,climate,by="DATE")
fulldata$streamflow_mean=as.numeric(fulldata$streamflow_mean)
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

fulldata2 = fulldata[which(fulldata$YEAR>=1950 & fulldata$YEAR<=2015),]
years2 = unique(fulldata2$YEAR)

cortable = NULL

for(j in 16:ncol(fulldata2)){
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

usevaridx = c(19,44,85,70,118,94,124,38,88,75,55,79,164,42)

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

######
# All vars and combos

xmato = as.matrix(traingroup[,usevaridx])
xmats = xmato^2
xmatr = sqrt(xmato)
xmata = xmato
xmean = apply(xmato,2,mean,na.rm=TRUE)
for(i in 1:ncol(xmata)) xmata[,i]=xmata[,i]-xmean[i]

xmatall = cbind(xmato,xmats)
xmatall = cbind(xmatall,xmatr)
xmatall = cbind(xmatall,xmata)

allnames = c(names(traingroup[,usevaridx]),
             paste(names(traingroup[,usevaridx]),"^2",sep=""),
             paste(names(traingroup[,usevaridx]),"_sqrt",sep=""),
             paste(names(traingroup[,usevaridx]),"_anom",sep=""))

#names(xmatall) = allnames


#options(show.error.messages = TRUE)
for(c in 0:(length(allnames)-1)){
  combos = combn(1:length(allnames),c)
  varswithin = list()
  corswithin = list()
  RMSEswithin = list()
  #pdf(paste("/home/woot0002/complexstreamflowregressiontests_",c,"varsout.pdf",sep=""),onefile=TRUE,width=5,height=5)
  
  if(c==0){
    endpoint=1
  } else {
    endpoint=ncol(combos)
  }
  
  for(i in 1:ncol(combos)){
    if(c!=0){
      xmat = xmatall[,-c(combos[,i])]
      varswithin[[i]] = allnames[-c(combos[,i])]
    } else {
      xmat=xmatall
      varswithin[[i]] = allnames
    }
    ymat = as.matrix(traingroup[,4])
    
    NAidx = which(is.na(ymat)==TRUE)
    if(length(NAidx)>0){
      xmat = xmat[-NAidx,]
      ymat = ymat[-NAidx,]
    }
    NAidx = which(is.na(xmat)==TRUE,arr.ind = TRUE)
    if(c<(length(allnames)-1)){
      NAidxu = unique(NAidx[,1])
    } else {
      NAidxu = NAidx
    }
    if(c==(length(allnames)-1)){
      xmat = as.matrix(xmat)
      ymat = as.matrix(ymat)
    }
    if(length(NAidxu)>0){
      xmat = xmat[-NAidxu,]
      ymat = ymat[-NAidxu]
    }
    if(c==(length(allnames)-1)){
      xmat = as.matrix(xmat)
      
    }
    ymat = as.matrix(ymat)
    try(fit <- lars(xmat,ymat , type="lasso"))
    #plot(fit)
    best_step <- fit$df[which.min(fit$RSS)]
    predictions <- predict(fit, xmat, s=best_step, type="fit")$fit
    predictions = ifelse(predictions<0,0,predictions)
    RMSEswithin[[i]] = sqrt(mean((predictions-ymat)^2))
    corswithin[[i]] = cor(ymat,predictions)
    #logpred = log(predictions)
    #plot(density(log(ymat)),main=paste("Observed vs. Modeled Streamflow \n",c," vars left out, version ",i," / ",ncol(combos)),xlab="ln(streamflow)")
    #try(lines(density(logpred),col="red"))
    #text(-2,0.5,labels=paste("cor = ",round(corswithin[[i]],4),sep=""),cex=1,pos = 4)
    #text(-2,0.45,labels=paste("RMSE = ",round(RMSEswithin[[i]],4),sep=""),cex=1,pos = 4)
    #legend("topright",legend=c("Observed","Modeled"),col=c("black","red"),lty=1)
  }
  save(list=c("varswithin","corswithin","RMSEswithin"),file=paste("/home/woot0002/streamflowtests_complexresults_",c,"varsleftout.Rdata",sep=""))
  
  send.mail(from = emadd,
            to = emadd,
            subject = "message from R on climatedata",
            body = paste("streamflowregtests_all.R has finished running for ",c," vars out",sep=""), 
            authenticate = TRUE,
            smtp = list(host.name = "smtp.office365.com", port = 587,
                        user.name = emadd, passwd = pswd, tls = TRUE))
  
  #varsusedall[[(c+1)]] = varswithin
  #corsusedall[[(c+1)]] = corswithin
  #RMSEsusedall[[(c+1)]] = RMSEswithin
  #dev.off()
}

#save(list=c("varsusedall","corsusedall","RMSEsusedall"),file="/home/woot0002/streamflowtests_complexresults.Rdata")
