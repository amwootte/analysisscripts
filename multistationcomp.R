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

gauges=read.table("/home/woot0002/usgsgaugelocations.csv",sep=",",header=TRUE)
gauges$river = c(rep("San Saba",3),rep("Llano",4))
streamflow_all = NULL
streamflowyear_all = NULL

for(i in 1:nrow(gauges)){

  gaugenum = paste("0",gauges[i,2],sep="")
  sitename = gauges[i,1]
  river = gauges[i,7]
  streamflow = read.table(paste("/home/woot0002/streamflow_",gaugenum,sep=""),header=TRUE,sep="\t",fill=TRUE)

  if(gaugenum!="08151500" & gaugenum!="08144500" & gaugenum!="08148500"){
    names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")
    streamidx = 4
  } else {
    names(streamflow) = c("agency","site","DATE","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC","streamflow_mean","streamflow_mean_QC")
    streamidx = 8
  }

  streamflow = streamflow[,c(1,2,3,streamidx)]
  streamflow$DATE = as.character(streamflow$DATE)
  streamflow$streamflow_mean=as.numeric(streamflow$streamflow_mean)
  tmp = c(NA,rollapply(streamflow$streamflow_mean,3,mean,na.rm=TRUE),NA)
  streamflow$streamflow_mean=tmp
  streamflow$YEAR = as.numeric(substr(streamflow$DATE,1,4))
  streamflow = subset(streamflow,YEAR>=1981 & YEAR<=2005)
  streamflow$sitename = sitename
  streamflow$river = river

  streamflow_all = rbind(streamflow_all,streamflow)
  
  streamflow_hist_year =  aggregate(streamflow[,4],by=list(year=streamflow$YEAR),mean,na.rm=TRUE)
  streamflow_hist_year$site = gaugenum
  streamflow_hist_year$sitename = sitename
  streamflow_hist_year$river = river
  names(streamflow_hist_year)[2]="streamflow_mean"
  streamflowyear_all = rbind(streamflowyear_all,streamflow_hist_year)
}


streamflowyear_all$sitename <- factor(streamflowyear_all$sitename,levels = gauges[,1],ordered = TRUE)
streamflowyear_all$river  <- factor(streamflowyear_all$river,levels = unique(gauges[,7]),ordered = TRUE)

streamflow_all$sitename <- factor(streamflow_all$sitename,levels = gauges[,1],ordered = TRUE)
streamflow_all$river  <- factor(streamflow_all$river,levels = unique(gauges[,7]),ordered = TRUE)


library(ggplot2)

ggplot(streamflowyear_all, aes(x=sitename, y=streamflow_mean,fill=river)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Yearly Average Streamflow (1981-2005)")+xlab("")+ylab("Average Streamflow (ft^3/s)") 

ggplot(streamflow_all, aes(x=sitename, y=streamflow_mean,fill=river)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Daily Average Streamflow (1981-2005)")+xlab("")+ylab("Average Streamflow (ft^3/s)") 

