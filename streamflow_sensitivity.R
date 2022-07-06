
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

gaugenum = "08146000"
loclon = -98.71947479 # gauging station location 
loclat = 31.2131691
futureperiod = c(2070,2099)

load(paste("/home/woot0002/CPREP_climo_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))
load(paste("/home/woot0002/CPREP_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))

######
# precipitation climo calcs - historical

climohist = histfilebreakdown
climohist$PRCP = NA
climohist$PRCPSD = NA
climohist$TMAX = NA
climohist$TMAXSD = NA
climohist$TMIN = NA
climohist$TMINSD = NA

climohist$HOT1S_1W = NA
climohist$HOT1S_1WSD = NA
climohist$HOT1S_3M = NA
climohist$HOT1S_3MSD = NA

climohist$HOT90_12M = NA
climohist$HOT90_12MSD = NA
climohist$HOT90_24M = NA
climohist$HOT90_24MSD = NA
climohist$HOT90_3M = NA
climohist$HOT90_3MSD = NA

climohist$MAXPRCP1W_12M = NA
climohist$MAXPRCP1W_12MSD = NA
climohist$MAXPRCP1W_6M = NA
climohist$MAXPRCP1W_6MSD = NA
climohist$MAXPRCP1W_2W = NA
climohist$MAXPRCP1W_2WSD = NA

climohist$MAXTMAX1W_2W = NA
climohist$MAXTMAX1W_2WSD = NA
climohist$MAXTMAX1W_3W = NA
climohist$MAXTMAX1W_3WSD = NA

climohist$MDRN_3M = NA
climohist$MDRN_3MSD = NA
climohist$MDRN_4W = NA
climohist$MDRN_4WSD = NA

climohist$PRCPAVG_3D = NA
climohist$PRCPAVG_3DSD = NA
climohist$PRCPAVG_6M = NA
climohist$PRCPAVG_6MSD = NA
climohist$PRCPAVG_12M = NA
climohist$PRCPAVG_12MSD = NA

climohist$TAVG_3M = NA
climohist$TAVG_3MSD = NA
climohist$TAVG_12M = NA
climohist$TAVG_12MSD = NA
climohist$TMAXAVG_1W = NA
climohist$TMAXAVG_1WSD = NA

for(i in 1:nrow(histfilebreakdown)){
  
  tmp = aggregate(histresults_climo[[i]]$TAVG_3M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$TAVG_3M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$TAVG_3M[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$TAVG_12M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$TAVG_12M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$TAVG_12MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$TMAXAVG_1W,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$TMAXAVG_1W[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$TMAXAVG_1W[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$PRCPAVG_3D,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$PRCPAVG_3D[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$PRCPAVG_3DSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$PRCPAVG_6M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$PRCPAVG_6M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$PRCPAVG_6MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$PRCPAVG_12M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$PRCPAVG_12M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$PRCPAVG_12MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$MDRN_3M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$MDRN_3M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$MDRN_3MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$MDRN_4W,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$MDRN_4W[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$MDRN_4WSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$MAXTMAX1W_2W,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$MAXTMAX1W_2W[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$MAXTMAX1W_2WSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$MAXTMAX1W_3W,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$MAXTMAX1W_3W[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$MAXTMAX1W_3WSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$MAXPRCP1W_12M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$MAXPRCP1W_12M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$MAXPRCP1W_12MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$MAXPRCP1W_6M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$MAXPRCP1W_6M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$MAXPRCP1W_6MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$MAXPRCP1W_2W,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$MAXPRCP1W_2W[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$MAXPRCP1W_2WSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$HOT90_12M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$HOT90_12M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$HOT90_12MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$HOT90_24M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$HOT90_24M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$HOT90_24MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$HOT90_3M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$HOT90_3M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$HOT90_3MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$HOT1S_1W,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$HOT1S_1W[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$HOT1S_1WSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$HOT1S_3M,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$HOT1S_3M[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$HOT1S_3MSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$PRCP,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),sum,na.rm=TRUE)
  climohist$PRCP[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$PRCPSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$TMAX,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$TMAX[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$TMAXSD[i] = sd(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(histresults_climo[[i]]$TMIN,by=list(year=as.numeric(substr(histresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climohist$TMIN[i] = mean(tmp[,2],na.rm=TRUE)
  climohist$TMINSD[i] = sd(tmp[,2],na.rm=TRUE)
}

######
# precipitation climo calcs - projected

climoproj = projfilebreakdown
climoproj$PRCP = NA
climoproj$R1MM = NA
climoproj$TMAX = NA
climoproj$TMIN = NA

climoproj$PRCP_anom = NA
climoproj$PRCP_Nanom = NA
climoproj$TMAX_anom = NA
climoproj$TMAX_Nanom = NA
climoproj$TMIN_anom = NA
climoproj$TMIN_Nanom = NA

climoproj$HOT90_12M_anom = NA
climoproj$HOT90_12M_Nanom = NA
climoproj$HOT90_24M_anom = NA
climoproj$HOT90_24M_Nanom = NA
climoproj$HOT90_3M_anom = NA
climoproj$HOT90_3M_Nanom = NA

climoproj$HOT1S_1W_anom = NA
climoproj$HOT1S_1W_Nanom = NA
climoproj$HOT1S_3M_anom = NA
climoproj$HOT1S_3M_Nanom = NA

climoproj$MAXPRCP1W_12M_anom = NA
climoproj$MAXPRCP1W_12M_Nanom = NA
climoproj$MAXPRCP1W_6M_anom = NA
climoproj$MAXPRCP1W_6M_Nanom = NA
climoproj$MAXPRCP1W_2W_anom = NA
climoproj$MAXPRCP1W_2W_Nanom = NA

climoproj$MAXTMAX1W_2W_anom = NA
climoproj$MAXTMAX1W_2W_Nanom = NA
climoproj$MAXTMAX1W_3W_anom = NA
climoproj$MAXTMAX1W_3W_Nanom = NA

climoproj$MDRN_3M_anom = NA
climoproj$MDRN_3M_Nanom = NA
climoproj$MDRN_4W_anom = NA
climoproj$MDRN_4W_Nanom = NA

climoproj$PRCPAVG_12M_anom = NA
climoproj$PRCPAVG_12M_Nanom = NA
climoproj$PRCPAVG_6M_anom = NA
climoproj$PRCPAVG_6M_Nanom = NA
climoproj$PRCPAVG_3D_anom = NA
climoproj$PRCPAVG_3D_Nanom = NA

climoproj$TAVG_12M_anom = NA
climoproj$TAVG_12M_Nanom = NA
climoproj$TAVG_3M_anom = NA
climoproj$TAVG_3M_Nanom = NA
climoproj$TMAXAVG_1W_anom = NA
climoproj$TMAXAVG_1W_Nanom = NA

for(i in 1:nrow(projfilebreakdown)){
  
  tmp = aggregate(projresults_climo[[i]]$TAVG_12M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$TAVG_12M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$TAVG_3M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$TAVG_3M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$TMAXAVG_1W,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$TMAXAVG_1W[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$PRCPAVG_12M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$PRCPAVG_12M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$PRCPAVG_6M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$PRCPAVG_6M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$PRCPAVG_3D,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$PRCPAVG_3D[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$MDRN_3M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$MDRN_3M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$MDRN_4W,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$MDRN_4W[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$MAXTMAX1W_2W,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$MAXTMAX1W_2W[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$MAXTMAX1W_3W,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$MAXTMAX1W_3W[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$MAXPRCP1W_12M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$MAXPRCP1W_12M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$MAXPRCP1W_6M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$MAXPRCP1W_6M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$MAXPRCP1W_2W,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$MAXPRCP1W_2W[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$HOT1S_1W,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$HOT1S_1W[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$HOT1S_3M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$HOT1S_3M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$HOT90_12M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$HOT90_12M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$HOT90_24M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$HOT90_24M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$HOT90_3M,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$HOT90_3M[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$PRCP,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),sum,na.rm=TRUE)
  climoproj$PRCP[i] = mean(tmp[,2],na.rm=TRUE)
 
  tmp = aggregate(projresults_climo[[i]]$TMAX,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$TMAX[i] = mean(tmp[,2],na.rm=TRUE)
  
  tmp = aggregate(projresults_climo[[i]]$TMIN,by=list(year=as.numeric(substr(projresults_climo[[i]]$DATE,1,4))),mean,na.rm=TRUE)
  climoproj$TMIN[i] = mean(tmp[,2],na.rm=TRUE)
  
  histidx = which(histfilebreakdown$GCM==projfilebreakdown$GCM[i] & histfilebreakdown$DS==projfilebreakdown$DS[i] & histfilebreakdown$obs==projfilebreakdown$obs[i])
  
  climoproj$PRCP_anom[i] = climoproj$PRCP[i]-climohist$PRCP[histidx]
  climoproj$PRCP_Nanom[i] = climoproj$PRCP_anom[i]/climohist$PRCPSD[histidx]
  
  climoproj$TMAX_anom[i] = climoproj$TMAX[i]-climohist$TMAX[histidx]
  climoproj$TMAX_Nanom[i] = climoproj$TMAX_anom[i]/climohist$TMAXSD[histidx]
  
  climoproj$TMIN_anom[i] = climoproj$TMIN[i]-climohist$TMIN[histidx]
  climoproj$TMIN_Nanom[i] = climoproj$TMIN_anom[i]/climohist$TMINSD[histidx]
  
  climoproj$HOT1S_1W_anom[i] = climoproj$HOT1S_1W[i]-climohist$HOT1S_1W[histidx]
  climoproj$HOT1S_1W_Nanom[i] = climoproj$HOT1S_1W_anom[i]/climohist$HOT1S_1WSD[histidx]
  
  climoproj$HOT1S_3M_anom[i] = climoproj$HOT1S_3M[i]-climohist$HOT1S_3M[histidx]
  climoproj$HOT1S_3M_Nanom[i] = climoproj$HOT1S_3M_anom[i]/climohist$HOT1S_3MSD[histidx]
  
  climoproj$HOT90_24M_anom[i] = climoproj$HOT90_24M[i]-climohist$HOT90_24M[histidx]
  climoproj$HOT90_24M_Nanom[i] = climoproj$HOT90_24M_anom[i]/climohist$HOT90_24MSD[histidx]
  
  climoproj$HOT90_12M_anom[i] = climoproj$HOT90_12M[i]-climohist$HOT90_12M[histidx]
  climoproj$HOT90_12M_Nanom[i] = climoproj$HOT90_12M_anom[i]/climohist$HOT90_12MSD[histidx]
  
  climoproj$HOT90_3M_anom[i] = climoproj$HOT90_3M[i]-climohist$HOT90_3M[histidx]
  climoproj$HOT90_3M_Nanom[i] = climoproj$HOT90_3M_anom[i]/climohist$HOT90_3MSD[histidx]
  
  climoproj$MAXPRCP1W_12M_anom[i] = climoproj$MAXPRCP1W_12M[i]-climohist$MAXPRCP1W_12M[histidx]
  climoproj$MAXPRCP1W_12M_Nanom[i] = climoproj$MAXPRCP1W_12M_anom[i]/climohist$MAXPRCP1W_12MSD[histidx]
  
  climoproj$MAXPRCP1W_2W_anom[i] = climoproj$MAXPRCP1W_2W[i]-climohist$MAXPRCP1W_2W[histidx]
  climoproj$MAXPRCP1W_2W_Nanom[i] = climoproj$MAXPRCP1W_2W_anom[i]/climohist$MAXPRCP1W_2WSD[histidx]
  
  climoproj$MAXPRCP1W_6M_anom[i] = climoproj$MAXPRCP1W_6M[i]-climohist$MAXPRCP1W_6M[histidx]
  climoproj$MAXPRCP1W_6M_Nanom[i] = climoproj$MAXPRCP1W_6M_anom[i]/climohist$MAXPRCP1W_6MSD[histidx]
  
  climoproj$MAXTMAX1W_2W_anom[i] = climoproj$MAXTMAX1W_2W[i]-climohist$MAXTMAX1W_2W[histidx]
  climoproj$MAXTMAX1W_2W_Nanom[i] = climoproj$MAXTMAX1W_2W_anom[i]/climohist$MAXTMAX1W_2WSD[histidx]
  
  climoproj$MAXTMAX1W_3W_anom[i] = climoproj$MAXTMAX1W_3W[i]-climohist$MAXTMAX1W_3W[histidx]
  climoproj$MAXTMAX1W_3W_Nanom[i] = climoproj$MAXTMAX1W_3W_anom[i]/climohist$MAXTMAX1W_3WSD[histidx]
  
  climoproj$MDRN_3M_anom[i] = climoproj$MDRN_3M[i]-climohist$MDRN_3M[histidx]
  climoproj$MDRN_3M_Nanom[i] = climoproj$MDRN_3M_anom[i]/climohist$MDRN_3MSD[histidx]
  
  climoproj$MDRN_4W_anom[i] = climoproj$MDRN_4W[i]-climohist$MDRN_4W[histidx]
  climoproj$MDRN_4W_Nanom[i] = climoproj$MDRN_4W_anom[i]/climohist$MDRN_4WSD[histidx]
  
  climoproj$PRCPAVG_12M_anom[i] = climoproj$PRCPAVG_12M[i]-climohist$PRCPAVG_12M[histidx]
  climoproj$PRCPAVG_12M_Nanom[i] = climoproj$PRCPAVG_12M_anom[i]/climohist$PRCPAVG_12MSD[histidx]
  
  climoproj$PRCPAVG_6M_anom[i] = climoproj$PRCPAVG_6M[i]-climohist$PRCPAVG_6M[histidx]
  climoproj$PRCPAVG_6M_Nanom[i] = climoproj$PRCPAVG_6M_anom[i]/climohist$PRCPAVG_6MSD[histidx]
  
  climoproj$PRCPAVG_3D_anom[i] = climoproj$PRCPAVG_3D[i]-climohist$PRCPAVG_3D[histidx]
  climoproj$PRCPAVG_3D_Nanom[i] = climoproj$PRCPAVG_3D_anom[i]/climohist$PRCPAVG_3DSD[histidx]
  
  
  climoproj$TAVG_12M_anom[i] = climoproj$TAVG_12M[i]-climohist$TAVG_12M[histidx]
  climoproj$TAVG_12M_Nanom[i] = climoproj$TAVG_12M_anom[i]/climohist$TAVG_12MSD[histidx]
  
  climoproj$TAVG_3M_anom[i] = climoproj$TAVG_3M[i]-climohist$TAVG_3M[histidx]
  climoproj$TAVG_3M_Nanom[i] = climoproj$TAVG_3M_anom[i]/climohist$TAVG_3MSD[histidx]
  
  climoproj$TMAXAVG_1W_anom[i] = climoproj$TMAXAVG_1W[i]-climohist$TMAXAVG_1W[histidx]
  climoproj$TMAXAVG_1W_Nanom[i] = climoproj$TMAXAVG_1W_anom[i]/climohist$TMAXAVG_1WSD[histidx]
  
}


######
# streamflow calcs - yearly

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

###
# yearly calcs

streamflow_hist_year =  aggregate(streamflow_hist[,2:29],by=list(year=as.numeric(substr(streamflow_hist$DATE,1,4))),mean,na.rm=TRUE)

streamflow_proj_year =  aggregate(streamflow_proj[,2:82],by=list(year=as.numeric(substr(streamflow_proj$DATE,1,4))),mean,na.rm=TRUE)

climoproj$SF = climoproj$SF_anom = climoproj$SF_Nanom = NA
climohist$SF = climohist$SFSD = NA
for(i in 2:ncol(streamflow_proj_year)){
  nametmp = names(streamflow_proj_year)[i]
  histidx = which(names(streamflow_hist_year)==substr(nametmp,7,nchar(nametmp)))
  
  climoproj$SF[i-1] = mean(streamflow_proj_year[,i],na.rm=TRUE)
  climohist$SF[histidx-2] = mean(streamflow_hist_year[,histidx],na.rm=TRUE)
  climohist$SFSD[histidx-2] = sd(streamflow_hist_year[,histidx],na.rm=TRUE)
  
  climoproj$SF_anom[i-1] = climoproj$SF[i-1]-climohist$SF[histidx-2]
  climoproj$SF_Nanom[i-1] = climoproj$SF_anom[i-1]/climohist$SFSD[histidx-2]
}

#####

climoprojin = subset(climoproj,scen=="rcp85")

#library(ggplot2)

#ggplot(climoprojin, aes(x=PRCP_Nanom, y=SF_Nanom))+geom_point(aes(fill=factor(GCM), shape=factor(DS),colour=factor(obs)),size=4,stroke=1.5)+
#  geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed") + 
#  scale_shape_manual(values=c(21, 22, 25)) +
#  scale_color_manual(aesthetics="colour",values=c('seashell4','yellow', 'purple')) +
#  scale_color_manual(aesthetics="fill",values=c('red','blue', 'green')) + ylab("Standardized Streamflow Anomaly") + 
#  xlab("Standardized Precipitation Anomaly") + ggtitle(paste("Standardized Anomalies: ",gaugenum,sep=""))+ylim(-3,0.5)+xlim(-3,0.5)
#ggsave(paste("/home/woot0002/streamflow_stdanomalies_",gaugenum,".pdf",sep=""),device = "pdf",width=8,height=7)

#ggplot(climoprojin, aes(x=PRCP_Nanom, y=SF_Nanom))+geom_point(aes(shape=factor(DS),fill=factor(obs)),size=3)+
#  geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed") + 
#  scale_shape_manual(values=c(21, 22, 25)) +
#  scale_color_manual(aesthetics="fill",values=c('seashell4','yellow', 'purple')) +
#  ylab("Standardized Streamflow Anomaly") + 
#  xlab("Standardized Precipitation Anomaly") + ggtitle(paste("Standardized Anomalies: ",gaugenum,sep=""))+facet_wrap(vars(GCM))+ylim(-3,0.5)+xlim(-3,0.5)


write.table(climoproj,file=paste("/home/woot0002/climotable_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],".csv",sep=""),row.names=FALSE,sep=",")
write.table(climohist,file=paste("/home/woot0002/climotable_",gaugenum,"_1981-2005.csv",sep=""),row.names=FALSE,sep=",")

