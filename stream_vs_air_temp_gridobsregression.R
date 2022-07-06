
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

streamflowID = "charliesdaily"
loclon = -99.51830566
loclat = 30.91122222
dataset = "DAYMETv4"

streamflow = read.table(paste("/home/woot0002/",streamflowID,".csv",sep=""),header=TRUE,sep=",",fill=TRUE) # watch it with reordered tables based on gauge ID
daily_stream=streamflow

idxout = which(daily_stream[,2]== -Inf)
if(length(idxout)>0){
  daily_stream=daily_stream[-idxout,]
}
names(daily_stream) = c("date","tmax","tmin")

daily_stream$date = as.character(as.Date(as.character(daily_stream$date),"%m/%d/%y"))

#####
# determine which grid cell to pull

filecheck = "/data4/data/OBS/DAYMETv4/na/tasmax/daymet_v4_daily_na_tmax_1980.nc"

nctest = nc_open(filecheck)
lon = ncvar_get(nctest,"lon")
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

###
# create model grid

LON = as.vector(lon)
LAT = as.vector(lat)
R= rep(1:dim(lon)[1],each=dim(lon)[2])
C= rep(1:dim(lon)[2],dim(lon)[1])

modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")

###
# get cells to use

pointarea = distfunc(loclon,loclat,modelgrid)

###
# pull tmax and tmin data

filesintmax = system("ls /data4/data/OBS/DAYMETv4/na/tasmax/*.nc",intern=TRUE)
filesintmin = system("ls /data4/data/OBS/DAYMETv4/na/tasmin/*.nc",intern=TRUE)

tmaxid = c(grep("2017",filesintmax),grep("2018",filesintmax),grep("2019",filesintmax),grep("2020",filesintmax),grep("2021",filesintmax))
tminid = c(grep("2017",filesintmin),grep("2018",filesintmin),grep("2019",filesintmin),grep("2020",filesintmin),grep("2021",filesintmin))

filesintmax = filesintmax[tmaxid]
filesintmin = filesintmin[tminid]

TMAX = TMIN = c()

for(i in 1:5){
  nctest = nc_open(filesintmax[i])
  TMAX = c(TMAX,ncvar_get(nctest,"tmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1)))
  nc_close(nctest)
  
  nctest = nc_open(filesintmin[i])
  TMIN = c(TMIN,ncvar_get(nctest,"tmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1)))
  nc_close(nctest)
}

dates = seq(as.Date("2017-01-01"),as.Date("2021-12-31"),by="day")
if(length(TMAX)<length(dates)){
  dates = dates[-which(substr(dates,6,10)=="02-29")]
}

climdat = data.frame(dates,TMAX,TMIN)

names(climdat) = c("date","TMAX","TMIN")

climdat$date = as.character(climdat$date)
daily_stream$date = as.character(daily_stream$date)

#which(as.character(climdat$date) %in% as.character(daily_stream$date))

####
# merge gridded obs with stream temperature data

alldata = merge(daily_stream,climdat,by="date")

####
# analyze results

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


plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$tmin),na.rm=TRUE),ylab="Temperature Deg C",main="Stream Temperatures")
lines(alldata$tmin,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

plot(alldata$tmax-alldata$tmin,type="l",col="blue",ylim=range(alldata$tmax-alldata$tmin,na.rm=TRUE),ylab="Temperature difference",main="Stream Temperatures")

####
# analyze results

alldata$TMAX3D = filter(alldata$TMAX,rep(1/3,3),sides=2)
alldata$TMIN3D = filter(alldata$TMIN,rep(1/3,3),sides=2)
alldata$TMAX5D = filter(alldata$TMAX,rep(1/5,5),sides=2)
alldata$TMIN5D = filter(alldata$TMIN,rep(1/5,5),sides=2)
alldata$TMAX7D = filter(alldata$TMAX,rep(1/7,7),sides=2)
alldata$TMIN7D = filter(alldata$TMIN,rep(1/7,7),sides=2)
alldata$TMAX10D = filter(alldata$TMAX,rep(1/10,10),sides=2)
alldata$TMIN10D = filter(alldata$TMIN,rep(1/10,10),sides=2)
alldata$TMAX14D = filter(alldata$TMAX,rep(1/14,14),sides=2)
alldata$TMIN14D = filter(alldata$TMIN,rep(1/14,14),sides=2)

cor(alldata$TMAX3D,alldata$tmax,use="pairwise.complete.obs")
cor(alldata$TMIN3D,alldata$tmin,use="pairwise.complete.obs")

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylab="Temperature Deg C",main="High Temperatures Stream vs. Air")
lines(alldata$TMAX3D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylab="Temperature Deg C",main="Low Temperatures Stream vs. Air")
lines(alldata$TMIN3D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

####

cor(alldata$TMAX5D,alldata$tmax,use="pairwise.complete.obs")
cor(alldata$TMIN5D,alldata$tmin,use="pairwise.complete.obs")

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylab="Temperature Deg C",main="High Temperatures Stream vs. Air")
lines(alldata$TMAX5D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylab="Temperature Deg C",main="Low Temperatures Stream vs. Air")
lines(alldata$TMIN5D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

###
cor(alldata$TMAX7D,alldata$tmax,use="pairwise.complete.obs")
cor(alldata$TMIN7D,alldata$tmin,use="pairwise.complete.obs")

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylab="Temperature Deg C",main="High Temperatures Stream vs. Air")
lines(alldata$TMAX7D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylab="Temperature Deg C",main="Low Temperatures Stream vs. Air")
lines(alldata$TMIN7D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

###

cor(alldata$TMAX10D,alldata$tmax,use="pairwise.complete.obs")
cor(alldata$TMIN10D,alldata$tmin,use="pairwise.complete.obs")

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylab="Temperature Deg C",main="High Temperatures Stream vs. Air")
lines(alldata$TMAX10D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylab="Temperature Deg C",main="Low Temperatures Stream vs. Air")
lines(alldata$TMIN10D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

###

cor(alldata$TMAX14D,alldata$tmax,use="pairwise.complete.obs")
cor(alldata$TMIN14D,alldata$tmin,use="pairwise.complete.obs")

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$TMAX),na.rm=TRUE),ylab="Temperature Deg C",main="High Temperatures Stream vs. Air")
lines(alldata$TMAX14D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$TMIN),na.rm=TRUE),ylab="Temperature Deg C",main="Low Temperatures Stream vs. Air")
lines(alldata$TMIN14D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Stream","Air"),col=c("blue","black"),lty=1)

####
# linear regression

lmfit = lm(tmax~TMAX,data=alldata)
alldata$tmaxfit = fitted.values(lmfit)

lmfit = lm(tmin~TMIN,data=alldata)
alldata$tminfit = fitted.values(lmfit)

lmfit = lm(tmax~TMAX3D,data=alldata)
alldata$tmaxfit3D = c(rep(NA,2),fitted.values(lmfit))

lmfit = lm(tmin~TMIN3D,data=alldata)
alldata$tminfit3D = c(rep(NA,2),fitted.values(lmfit))

lmfit = lm(tmax~TMAX5D,data=alldata)
alldata$tmaxfit5D = c(rep(NA,4),fitted.values(lmfit))
save(list="lmfit",file="TMAX5daylmfit_charlies.Rdata")

lmfit = lm(tmin~TMIN5D,data=alldata)
alldata$tminfit5D = c(rep(NA,4),fitted.values(lmfit))
save(list="lmfit",file="TMIN5daylmfit_charlies.Rdata")


lmfit = lm(tmax~TMAX7D,data=alldata)
alldata$tmaxfit7D = c(rep(NA,6),fitted.values(lmfit))

lmfit = lm(tmin~TMIN7D,data=alldata)
alldata$tminfit7D = c(rep(NA,6),fitted.values(lmfit))

lmfit = lm(tmax~TMAX10D,data=alldata)
alldata$tmaxfit10D = c(rep(NA,9),fitted.values(lmfit))

lmfit = lm(tmin~TMIN10D,data=alldata)
alldata$tminfit10D = c(rep(NA,9),fitted.values(lmfit))

lmfit = lm(tmax~TMAX14D,data=alldata)
alldata$tmaxfit14D = c(rep(NA,13),fitted.values(lmfit))

lmfit = lm(tmin~TMIN14D,data=alldata)
alldata$tminfit14D = c(rep(NA,13),fitted.values(lmfit))

####
# plot results

par(mfrow=c(2,1))

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$tmaxfit),na.rm=TRUE),ylab="Temperature Deg C",main="High Stream Temperatures - obs vs. fitted")
lines(alldata$tmaxfit,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tmaxfit-alldata$tmax,type="l",col="blue",ylab="Temperature Diff Deg C",main="High Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmax,alldata$tmaxfit,use="pairwise.complete.obs")

###

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$tminfit),na.rm=TRUE),ylab="Temperature Deg C",main="Low Stream Temperatures - obs vs. fitted")
lines(alldata$tminfit,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tminfit-alldata$tmin,type="l",col="blue",ylab="Temperature Diff Deg C",main="Low Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmin,alldata$tminfit,use="pairwise.complete.obs")

####

par(mfrow=c(2,1))

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$tmaxfit3D),na.rm=TRUE),ylab="Temperature Deg C",main="High Stream Temperatures - obs vs. fitted")
lines(alldata$tmaxfit3D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tmaxfit3D-alldata$tmax,type="l",col="blue",ylab="Temperature Diff Deg C",main="High Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmax,alldata$tmaxfit3D,use="pairwise.complete.obs")

###

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$tminfit3D),na.rm=TRUE),ylab="Temperature Deg C",main="Low Stream Temperatures - obs vs. fitted")
lines(alldata$tminfit3D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tminfit3D-alldata$tmin,type="l",col="blue",ylab="Temperature Diff Deg C",main="Low Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmin,alldata$tminfit3D,use="pairwise.complete.obs")

####

par(mfrow=c(2,1))

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$tmaxfit5D),na.rm=TRUE),ylab="Temperature Deg C",main="High Stream Temperatures - obs vs. fitted")
lines(alldata$tmaxfit5D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tmaxfit5D-alldata$tmax,type="l",col="blue",ylab="Temperature Diff Deg C",main="High Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmax,alldata$tmaxfit5D,use="pairwise.complete.obs")

###

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$tminfit5D),na.rm=TRUE),ylab="Temperature Deg C",main="Low Stream Temperatures - obs vs. fitted")
lines(alldata$tminfit5D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tminfit5D-alldata$tmin,type="l",col="blue",ylab="Temperature Diff Deg C",main="Low Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmin,alldata$tminfit5D,use="pairwise.complete.obs")

####

par(mfrow=c(2,1))

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$tmaxfit7D),na.rm=TRUE),ylab="Temperature Deg C",main="High Stream Temperatures - obs vs. fitted")
lines(alldata$tmaxfit7D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tmaxfit7D-alldata$tmax,type="l",col="blue",ylab="Temperature Diff Deg C",main="High Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmax,alldata$tmaxfit7D,use="pairwise.complete.obs")

###

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$tminfit7D),na.rm=TRUE),ylab="Temperature Deg C",main="Low Stream Temperatures - obs vs. fitted")
lines(alldata$tminfit7D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tminfit7D-alldata$tmin,type="l",col="blue",ylab="Temperature Diff Deg C",main="Low Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmin,alldata$tminfit7D,use="pairwise.complete.obs")

####

par(mfrow=c(2,1))

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$tmaxfit10D),na.rm=TRUE),ylab="Temperature Deg C",main="High Stream Temperatures - obs vs. fitted")
lines(alldata$tmaxfit10D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tmaxfit10D-alldata$tmax,type="l",col="blue",ylab="Temperature Diff Deg C",main="High Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmax,alldata$tmaxfit10D,use="pairwise.complete.obs")

###

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$tminfit10D),na.rm=TRUE),ylab="Temperature Deg C",main="Low Stream Temperatures - obs vs. fitted")
lines(alldata$tminfit10D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tminfit10D-alldata$tmin,type="l",col="blue",ylab="Temperature Diff Deg C",main="Low Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmin,alldata$tminfit10D,use="pairwise.complete.obs")

####

par(mfrow=c(2,1))

plot(alldata$tmax,type="l",col="blue",ylim=range(c(alldata$tmax,alldata$tmaxfit14D),na.rm=TRUE),ylab="Temperature Deg C",main="High Stream Temperatures - obs vs. fitted")
lines(alldata$tmaxfit14D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tmaxfit14D-alldata$tmax,type="l",col="blue",ylab="Temperature Diff Deg C",main="High Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmax,alldata$tmaxfit14D,use="pairwise.complete.obs")

###

plot(alldata$tmin,type="l",col="blue",ylim=range(c(alldata$tmin,alldata$tminfit14D),na.rm=TRUE),ylab="Temperature Deg C",main="Low Stream Temperatures - obs vs. fitted")
lines(alldata$tminfit14D,type="l",col="black")
abline(h=0,lty=2)
legend("bottomleft",legend=c("Obs.","Fitted"),col=c("blue","black"),lty=1)

plot(alldata$tminfit14D-alldata$tmin,type="l",col="blue",ylab="Temperature Diff Deg C",main="Low Stream Temperatures - residuals")
abline(h=0,lty=2)

cor(alldata$tmin,alldata$tminfit14D,use="pairwise.complete.obs")

