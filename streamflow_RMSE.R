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

streamflow = read.table(paste("/home/woot0002/streamflow_",streamflowID,sep=""),header=TRUE,sep="\t",fill=TRUE) # watch it with reordered tables based on gauge ID

if(streamflowID!="08144500"){
  names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")
  streamidx = 4
} else {
  names(streamflow) = c("agency","site","DATE","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC","streamflow_mean","streamflow_mean_QC")
  streamidx = 8
}

streamflow$streamflowm = c(NA,rollapply(streamflow$streamflow_mean,3,mean,na.rm=TRUE),NA)

####

load(file=paste("/home/woot0002/CPREP_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))

datesin = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

for(i in 1:27){
  
  tmp = histresults[[i]][[1]]
  
  if(length(tmp)<length(datesin)){
    datesuse = datesin[-which(substr(datesin,6,10)=="02-29")]
  } else {
    datesuse = datesin
  }
  
  modname = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")
  
  if(i==1){
    histdattab = data.frame(datesuse,tmp)
    names(histdattab) = c("date",modname)
  } else {
    tmpframe = data.frame(datesuse,tmp)
    names(tmpframe) = c("date",modname)
    histdattab = merge(histdattab,tmpframe,by="date")
  }
}

######

streamdat = streamflow[,c(3,10)]
names(streamdat) = c("date","obs")

streamdat$date = as.character(streamdat$date)
histdattab$date = as.character(histdattab$date)

histdattab = merge(streamdat,histdattab,by="date")

RMSEs = c()

for(i in 3:29){
  
  diff = histdattab[,i]-histdattab[,2]
  diffsq = diff^2
  
  RMSEs[(i-2)] = sqrt(mean(diffsq,na.rm=TRUE))
}

mean(RMSEs)

####

histdattabyearly = aggregate(histdattab[,2:29],by=list(year=as.numeric(substr(histdattab$date,1,4))),mean,na.rm=TRUE)

RMSEs = c()

for(i in 3:29){
  
  diff = histdattabyearly[,i]-histdattabyearly[,2]
  diffsq = diff^2
  
  RMSEs[(i-2)] = sqrt(mean(diffsq,na.rm=TRUE))
}

mean(RMSEs)

####

# climatology of yearly average streamflow ensemble RMSE
climotab = apply(histdattabyearly[,2:ncol(histdattabyearly)],2,mean,na.rm=TRUE)

diff = climotab[2:length(climotab)]-climotab[1]
diffsq = diff^2

RMSEval = sqrt(mean(diffsq,na.rm=TRUE))
RMSEval

sqrt(diff[1:3]^2)