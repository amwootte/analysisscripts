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

streamflowID = "08144500"
climateID = "2562869"
loclon = -99.78561401
loclat = 30.91913605

prfile_hist = c(system("ls /data2/3to5/I35/pr/Daymet/pr*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/Livneh/pr*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PRISM/pr*historical*.nc",intern=TRUE))

prfile = "/data2/3to5/I35/pr/PRISM/pr_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

#####
# get grid precipitation at location

nctest = nc_open(prfile)
lon = ncvar_get(nctest,"lon")-360
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

# create model grid

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

# get cells to use

pointarea = distfunc(loclon,loclat,modelgrid)

# get precip data

obsnames = c("Daymet","Livneh","PRISM")

for(i in 1:3){
  nctest = nc_open(prfile_hist[i])
  PRCP = ncvar_get(nctest,"pr",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))*86400
  nc_close(nctest)
  
  dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

  if(length(dates)>length(PRCP)){
    datesin = dates[-which(substr(dates,6,10)=="01-29")]
  } else {
    datesin = dates
  }
    
  if(i==1){
    outputframe = data.frame(datesin,PRCP)
    names(outputframe)[2] = obsnames[i]
  } else {
    tmpframe = data.frame(datesin,PRCP)
    names(tmpframe)[2] = obsnames[i]
    outputframe = merge(outputframe,tmpframe,by="datesin")
  }
}

#####

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
fulldata = fulldata[,c(1:3,10,11,which(names(fulldata)=="streamflow_mean"),which(names(fulldata)=="PRCP"))]

fulldata$streamflow_mean=as.numeric(fulldata$streamflow_mean)
fulldata$PRCP=as.numeric(fulldata$PRCP)

fulldata$DATE = as.Date(fulldata$DATE)

names(outputframe)[1] ="DATE"
fulldata = merge(fulldata,outputframe,by="DATE")

fulldata$PRCP = fulldata$PRCP*25.4

fulldata$streamflow_mean3d = c(NA,rollapply(fulldata$streamflow_mean,3,mean,na.rm=TRUE),NA)

fulldata$PRCP_1M = c(rep(NA,29),rollapply(fulldata$PRCP,width=30,FUN=sum,na.rm=TRUE))
fulldata$Daymet_1M = c(rep(NA,29),rollapply(fulldata$Daymet,width=30,FUN=sum,na.rm=TRUE))
fulldata$Livneh_1M = c(rep(NA,29),rollapply(fulldata$Livneh,width=30,FUN=sum,na.rm=TRUE))
fulldata$PRISM_1M = c(rep(NA,29),rollapply(fulldata$PRISM,width=30,FUN=sum,na.rm=TRUE))

fulldata$SFM_1M = c(rep(NA,29),rollapply(fulldata$streamflow_mean,width=30,FUN=mean,na.rm=TRUE))

####
# find whiplash events

### station precipitation
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


### Daymet precipitation
percentile_bottom = quantile(fulldata$Daymet_1M,probs=0.25,na.rm=TRUE)
percentile_top = quantile(fulldata$Daymet_1M,probs=0.75,na.rm=TRUE)

precip_hardcode = fulldata$Daymet_1M
precip_hardcode=ifelse(fulldata$Daymet_1M<=percentile_bottom | fulldata$Daymet_1M>=percentile_top,ifelse(fulldata$Daymet_1M<=percentile_bottom,0,2),1)

tmp = precip_hardcode[1:(length(fulldata$Daymet_1M)-10)]-precip_hardcode[11:length(fulldata$Daymet_1M)]
tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,10))
tmppd = c(ifelse(tmp==2,1,0),rep(NA,10))
PDcount_Daymet = sum(tmppd,na.rm=TRUE)
DPcount_Daymet = sum(tmpdp,na.rm=TRUE)

fulldata$Daymet_DP = tmpdp
fulldata$Daymet_PD = tmppd
fulldata$Daymet_WP = c(tmp,rep(NA,10))



### Livneh precipitation
percentile_bottom = quantile(fulldata$Livneh_1M,probs=0.25,na.rm=TRUE)
percentile_top = quantile(fulldata$Livneh_1M,probs=0.75,na.rm=TRUE)

precip_hardcode = fulldata$Livneh_1M
precip_hardcode=ifelse(fulldata$Livneh_1M<=percentile_bottom | fulldata$Livneh_1M>=percentile_top,ifelse(fulldata$Livneh_1M<=percentile_bottom,0,2),1)

tmp = precip_hardcode[1:(length(fulldata$Livneh_1M)-10)]-precip_hardcode[11:length(fulldata$Livneh_1M)]
tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,10))
tmppd = c(ifelse(tmp==2,1,0),rep(NA,10))
PDcount_Livneh = sum(tmppd,na.rm=TRUE)
DPcount_Livneh = sum(tmpdp,na.rm=TRUE)

fulldata$Livneh_DP = tmpdp
fulldata$Livneh_PD = tmppd
fulldata$Livneh_WP = c(tmp,rep(NA,10))


### PRISM precipitation
percentile_bottom = quantile(fulldata$PRISM_1M,probs=0.25,na.rm=TRUE)
percentile_top = quantile(fulldata$PRISM_1M,probs=0.75,na.rm=TRUE)

precip_hardcode = fulldata$PRISM_1M
precip_hardcode=ifelse(fulldata$PRISM_1M<=percentile_bottom | fulldata$PRISM_1M>=percentile_top,ifelse(fulldata$PRISM_1M<=percentile_bottom,0,2),1)

tmp = precip_hardcode[1:(length(fulldata$PRISM_1M)-10)]-precip_hardcode[11:length(fulldata$PRISM_1M)]
tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,10))
tmppd = c(ifelse(tmp==2,1,0),rep(NA,10))
PDcount_PRISM = sum(tmppd,na.rm=TRUE)
DPcount_PRISM = sum(tmpdp,na.rm=TRUE)

fulldata$PRISM_DP = tmpdp
fulldata$PRISM_PD = tmppd
fulldata$PRISM_WP = c(tmp,rep(NA,10))


### streamflow 
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
# write out file

save(list=c("fulldata"),file=paste("/home/woot0002/whiplashinfo_",streamflowID,".Rdata",sep=""))




