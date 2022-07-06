
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

streamflowID = "boisdearcdaily"
climateID = "2444385"

streamflow = read.table(paste("/home/woot0002/",streamflowID,".csv",sep=""),header=TRUE,sep=",",fill=TRUE) # watch it with reordered tables based on gauge ID
climate = read.table(paste("/home/woot0002/",climateID,".csv",sep=""),header=TRUE,sep=",")

daily_stream=streamflow

idxout = which(daily_stream[,2]== -Inf)
if(length(idxout)>0){
daily_stream=daily_stream[-idxout,]
}
names(daily_stream) = c("date","tmax","tmin")

daily_stream$date = as.character(as.Date(as.character(daily_stream$date),"%m/%d/%y"))

names(climate)[3] = "date"
climate$TMAX = (climate$TMAX-32)*5/9
climate$TMIN = (climate$TMIN-32)*5/9
climate$TOBS = (climate$TOBS-32)*5/9

climate$date = as.character(climate$date)

alldata = merge(daily_stream,climate,by="date")

alldata$tavg = (alldata$tmax+alldata$tmin)/2

quants = quantile(alldata$tavg,probs=seq(0.0001,0.9999,by=0.0001),na.rm=TRUE)
thres = 30.8
quantdiff = quants-thres
minidx = which(abs(quantdiff)==min(abs(quantdiff)))

alldata$tavg4d = filter(alldata$tavg,rep(1/4,4),sides=2)
quants = quantile(alldata$tavg4d,probs=seq(0.0001,0.9999,by=0.0001),na.rm=TRUE)
thres = 32.5
quantdiff = quants-thres
minidx = which(abs(quantdiff)==min(abs(quantdiff)))

alldata$g24LT05 = ifelse(alldata$tavg>=30.8,1,0)
alldata$j24LT05 = ifelse(alldata$tavg>=32.5,1,0)

yearlycalcs = aggregate(alldata[,2:3],by=list(year=substr(alldata$date,1,4)),mean,na.rm=TRUE)

yearlycalcs = aggregate(alldata[,18:19],by=list(year=substr(alldata$date,1,4)),sum,na.rm=TRUE)


cor(alldata$TMAX,alldata$tmax,use="pairwise.complete.obs")

cor(alldata$TMIN,alldata$tmin,use="pairwise.complete.obs")

plot(alldata$tmax~alldata$TMAX,type="p",xlim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylab="Streamflow high temperature",xlab="Air High Temperature")
abline(coef=c(0,1),lty=2)

plot(alldata$tmin~alldata$TMIN,type="p",xlim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylab="Streamflow Low temperature",xlab="Air Low Temperature")
abline(coef=c(0,1),lty=2)

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylab="Temperature Deg C",main="High Temperatures Stream vs. Air")
lines(alldata$TMAX,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylab="Temperature Deg C",main="Low Temperatures Stream vs. Air")
lines(alldata$TMIN,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)
