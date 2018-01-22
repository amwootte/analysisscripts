################
#
# Analysis for Variance Check

############
# library and function load

############
# load libraries

###############3
# 
# Variance Check

library(ncdf)
library(fields)
library(maps)

############
# load calculator functions
source("componentfunctions.R")

###
# set variable name
varname = "tmax"

###
# set extra options
printPDF = TRUE # if you want some plots printed to PDF set this to TRUE

##############
# Read in components from netcdf file

compfile = paste("uncertainty/SE/analysisid/",varname,"_components_fix.nc",sep="")
# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file
test = open.ncdf(compfile)
lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
time = 1950:2099
V = get.var.ncdf(test,"V")
M = get.var.ncdf(test,"M")
D = get.var.ncdf(test,"D")
S = get.var.ncdf(test,"S")
G = get.var.ncdf(test,"G")
T = get.var.ncdf(test,"T")
close.ncdf(test)

################
# Subset to future period only

timeidx = which(time>=2006)
M=M[,,timeidx]
D=D[,,timeidx]
S=S[,,timeidx]
G=G[,,timeidx]
T=T[,,timeidx]

################
# Get multiple domain tseries

# Central Florida= -84 -> -80, 27 -> 29
# Southern Appalachians= -85.5 -> -81, 34.5 -> 37
# Southwest Missouri= -95.5 -> -93, 35.5 -> 39

lonidx1 = which(lon>= -84 & lon<= -80)
latidx1 = which(lat>= 27 & lat>= 29)

lonidx2 = which(lon>= -85.5 & lon<= -81)
latidx2 = which(lat>= 34.5 & lat>= 34)

lonidx3 = which(lon>= -95.5 & lon<= -93)
latidx3 = which(lat>= 35.5 & lat>= 39)

Vts1 = mean(V[lonidx1,latidx1],na.rm=TRUE) # central Florida time series
Mts1 = apply(M[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Dts1 = apply(D[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Sts1 = apply(S[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Gts1 = apply(G[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Tts1 = apply(T[lonidx1,latidx1,],3,mean,na.rm=TRUE)

Vts2 = mean(V[lonidx2,latidx2],na.rm=TRUE) # Southern Appalachians time series
Mts2 = apply(M[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Dts2 = apply(D[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Sts2 = apply(S[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Gts2 = apply(G[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Tts2 = apply(T[lonidx2,latidx2,],3,mean,na.rm=TRUE)

Vts3 = mean(V[lonidx3,latidx3],na.rm=TRUE) # Southwest Missouri time series
Mts3 = apply(M[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Dts3 = apply(D[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Sts3 = apply(S[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Gts3 = apply(G[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Tts3 = apply(T[lonidx3,latidx3,],3,mean,na.rm=TRUE)

###################
###################
# Get anomaly values

files = system(paste("ls uncertainty/SE/anomdatid/",varname,"*.nc",sep=""),intern=TRUE)

if(varname=="pr" & length(grep("pr25",files))>0){
files = files[-grep("pr25",files)]
}

if(varname=="tmax" & length(grep("tmax95",files))>0){
files = files[-grep("tmax95",files)]
}

if(varname=="tmin" & length(grep("tmin32",files))>0){
files = files[-grep("tmin32",files)]
}

files = files[-grep("PRISM",files)]

CFts = list()
SAts = list()
SMts = list()

for(f in 1:length(files)){

test = open.ncdf(files[f])

if(f==1){
lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
}

vardata = get.var.ncdf(test,paste(varname,"_anom",sep=""))

spunits = test$var[[3]]$units

time = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
close.ncdf(test)

#if(datasetname=="maca" |  datasetname=="cmip5.bcca"){
##vardata2 = vardata
#vardata[,1,] = NA
#vardata[1,,]=NA
#} 

#if(datasetname=="dcp"){
#mask = ifelse(is.na(vardata[,,100])==TRUE,0,1)
#for(i in 1:56){
#vardata[,,i] = ifelse(mask==0,NA,vardata[,,i])
#}
#vardata[,,11] = NA
#}

#if(datasetname=="ccr" & (varname=="tmin" | varname=="tmax")){
#vardata = ifelse(vardata < -5,NA,vardata)
#}

CFts[[f]]=apply(vardata[lonidx1,latidx1,timeidx],3,mean,na.rm=TRUE)
SAts[[f]]=apply(vardata[lonidx2,latidx2,timeidx],3,mean,na.rm=TRUE)
SMts[[f]]=apply(vardata[lonidx3,latidx3,timeidx],3,mean,na.rm=TRUE)

message("Finished time series grabs for file ",f," / ",length(files))

}

#CFts2 = lapply(CFts,filter,sides=2,filter=rep(1,10)/10)
#SAts2 = lapply(SAts,filter,sides=2,filter=rep(1,10)/10)
#SMts2 = lapply(SMts,filter,sides=2,filter=rep(1,10)/10)

#Gts1 = ifelse(Gts1==0,NA,Gts1)
#Gts2 = ifelse(Gts2==0,NA,Gts2)
#Gts3 = ifelse(Gts3==0,NA,Gts2)

######################
###################
# Plotting results

if(varname=="tmax"){
  titlevar = "Average High Temperature"
  titlevarunits = "(C)"
}

if(varname=="tmax95"){
  titlevar = "Average Number of Days over 95F"
  titlevarunits = "days"
}

if(varname=="tmin"){
  titlevar = "Average Low Temperature"
  titlevarunits = "(C)"
}

if(varname=="tmin32"){
  titlevar = "Average Number of Days under 32F"
  titlevarunits = "days"
}


if(varname=="pr"){
  titlevar = "Average Total Precipitation"
  titlevarunits = "(%)"
}

if(varname=="pr25"){
  titlevar = "Average Number of Days with Rainfall > 1 in"
  titlevarunits = "days"
}

pdf(paste("uncertainty/SE/analysisid/",varname,"_variancecheck_fix.pdf",sep=""),height=6,width=6,onefile=TRUE)

############
# Central Florida

###########
# Fractional Uncertainty

year = 2006:2099

###
# Central Florida

####
# ens. sd overlaid with estimated sd

CFtab = do.call(rbind,CFts)
CFmean = apply(CFtab,2,mean,na.rm=TRUE)
CFsd = apply(CFtab,2,sd,na.rm=TRUE)
CFub = CFmean+2*CFsd
CFlb = CFmean-2*CFsd
CFsub = Gts1+2*sqrt(Tts1)
CFslb = Gts1-2*sqrt(Tts1)

for(f in 1:length(files)){
if(f==1){
plot(CFts[[f]]~year,col="gray",lwd=0.5,main=paste("Central Florida Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(CFtab,na.rm=TRUE),max(CFtab,na.rm=TRUE)))
} else {
lines(CFts[[f]]~year,col="gray",lwd=0.5)
}
}
lines(CFmean~year,col="black",lwd=2)
lines(CFub~year,col="black",lwd=2,lty=2)
lines(CFlb~year,col="black",lwd=2,lty=2)
lines(Gts1~year,col="purple",lwd=2)
lines(CFsub~year,col="purple",lwd=2,lty=2)
lines(CFslb~year,col="purple",lwd=2,lty=2)
legend("topleft",c("Ensemble members","Ens. Mean","Ens. Upper / Lower bounds","Derived Mean","Derived Upper / Lower Bounds"),col=c("gray","black","black","purple","purple"),lwd=c(0.5,2,2,2,2),lty=c(1,1,2,1,2))

####
# SD plot

plot(sqrt(Tts1)~CFsd,main=paste("Central Florida Derived vs Ensemble SD \nfor ",titlevar,sep=""),pch=19,ylab="Derived SD",xlab="Ensemble SD",ylim=c(0,max(CFsd,sqrt(Tts1),na.rm=TRUE)),xlim=c(0,max(CFsd,sqrt(Tts1),na.rm=TRUE)))
abline(coef=c(0,1))

####
# contribution overlaid

CFsub1 = Gts1+2*sqrt(Vts1)
CFslb1 = Gts1-2*sqrt(Vts1)
CFsub2 = Gts1+2*(sqrt(Vts1+Mts1))
CFslb2 = Gts1-2*(sqrt(Vts1+Mts1))
CFsub3 = Gts1+2*(sqrt(Vts1+Mts1+Dts1))
CFslb3 = Gts1-2*(sqrt(Vts1+Mts1+Dts1))
CFsub4 = Gts1+2*(sqrt(Vts1+Mts1+Dts1+Sts1))
CFslb4 = Gts1-2*(sqrt(Vts1+Mts1+Dts1+Sts1))

for(f in 1:length(files)){
if(f==1){
plot(CFts[[f]]~year,col="gray",lwd=0.5,main=paste("Central Florida Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(CFtab,na.rm=TRUE),max(CFtab,na.rm=TRUE)))
} else {
lines(CFts[[f]]~year,col="gray",lwd=0.5)
}
}

lines(Gts1~year,col="purple",lwd=2)
lines(CFsub1~year,col="orange",lwd=2,lty=2)
lines(CFslb1~year,col="orange",lwd=2,lty=2)
lines(CFsub2~year,col="blue",lwd=2,lty=2)
lines(CFslb2~year,col="blue",lwd=2,lty=2)
lines(CFsub3~year,col="red",lwd=2,lty=2)
lines(CFslb3~year,col="red",lwd=2,lty=2)
lines(CFsub4~year,col="darkgreen",lwd=2,lty=2)
lines(CFslb4~year,col="darkgreen",lwd=2,lty=2)

legend("topleft",c("Ensemble members","Derived Mean","Natural Variability","GCM Model Uncertainty","Downscaling Uncertainty","Scenario Uncertainty"),col=c("gray","purple","orange","blue","red","darkgreen"),lwd=c(0.5,2,2,2,2,2),lty=c(1,2,2,2,2,2))


###
# Southern Appalachians

####
# ens. sd overlaid with estimated sd

SAtab = do.call(rbind,SAts)
SAmean = apply(SAtab,2,mean,na.rm=TRUE)
SAsd = apply(SAtab,2,sd,na.rm=TRUE)
SAub = SAmean+2*SAsd
SAlb = SAmean-2*SAsd
SAsub = Gts2+2*sqrt(Tts2)
SAslb = Gts2-2*sqrt(Tts2)

for(f in 1:length(files)){
if(f==1){
plot(SAts[[f]]~year,col="gray",lwd=0.5,main=paste("Southern Appalachians Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SAtab,na.rm=TRUE),max(SAtab,na.rm=TRUE)))
} else {
lines(SAts[[f]]~year,col="gray",lwd=0.5)
}
}
lines(SAmean~year,col="black",lwd=2)
lines(SAub~year,col="black",lwd=2,lty=2)
lines(SAlb~year,col="black",lwd=2,lty=2)
lines(Gts2~year,col="purple",lwd=2)
lines(SAsub~year,col="purple",lwd=2,lty=2)
lines(SAslb~year,col="purple",lwd=2,lty=2)
legend("topleft",c("Ensemble members","Ens. Mean","Ens. Upper / Lower bounds","Derived Mean","Derived Upper / Lower Bounds"),col=c("gray","black","black","purple","purple"),lwd=c(0.5,2,2,2,2),lty=c(1,1,2,1,2))

####
# SD plot

plot(sqrt(Tts2)~SAsd,main=paste("Southern Appalachians Derived vs Ensemble SD \nfor ",titlevar,sep=""),pch=19,ylab="Derived SD",xlab="Ensemble SD",ylim=c(0,max(SAsd,sqrt(Tts2),na.rm=TRUE)),xlim=c(0,max(SAsd,sqrt(Tts2),na.rm=TRUE)))
abline(coef=c(0,1))

####
# contribution overlaid

SAsub1 = Gts2+2*sqrt(Vts2)
SAslb1 = Gts2-2*sqrt(Vts2)
SAsub2 = Gts2+2*(sqrt(Vts2+Mts2))
SAslb2 = Gts2-2*(sqrt(Vts2+Mts2))
SAsub3 = Gts2+2*(sqrt(Vts2+Mts2+Dts2))
SAslb3 = Gts2-2*(sqrt(Vts2+Mts2+Dts2))
SAsub4 = Gts2+2*(sqrt(Vts2+Mts2+Dts2+Sts2))
SAslb4 = Gts2-2*(sqrt(Vts2+Mts2+Dts2+Sts2))

for(f in 1:length(files)){
if(f==1){
plot(SAts[[f]]~year,col="gray",lwd=0.5,main=paste("Central Florida Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SAtab,na.rm=TRUE),max(SAtab,na.rm=TRUE)))
} else {
lines(SAts[[f]]~year,col="gray",lwd=0.5)
}
}

lines(Gts2~year,col="purple",lwd=2)
lines(SAsub1~year,col="orange",lwd=2,lty=2)
lines(SAslb1~year,col="orange",lwd=2,lty=2)
lines(SAsub2~year,col="blue",lwd=2,lty=2)
lines(SAslb2~year,col="blue",lwd=2,lty=2)
lines(SAsub3~year,col="red",lwd=2,lty=2)
lines(SAslb3~year,col="red",lwd=2,lty=2)
lines(SAsub4~year,col="darkgreen",lwd=2,lty=2)
lines(SAslb4~year,col="darkgreen",lwd=2,lty=2)

legend("topleft",c("Ensemble members","Derived Mean","Natural Variability","GCM Model Uncertainty","Downscaling Uncertainty","Scenario Uncertainty"),col=c("gray","purple","orange","blue","red","darkgreen"),lwd=c(0.5,2,2,2,2,2),lty=c(1,2,2,2,2,2))

###
# Southwest Missouri

####
# ens. sd overlaid with estimated sd

SMtab = do.call(rbind,SMts)
SMmean = apply(SMtab,2,mean,na.rm=TRUE)
SMsd = apply(SMtab,2,sd,na.rm=TRUE)
SMub = SMmean+2*SMsd
SMlb = SMmean-2*SMsd
SMsub = Gts3+2*sqrt(Tts3)
SMslb = Gts3-2*sqrt(Tts3)

for(f in 1:length(files)){
if(f==1){
plot(SMts[[f]]~year,col="gray",lwd=0.5,main=paste("Southwest Missouri Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SMtab,na.rm=TRUE),max(SMtab,na.rm=TRUE)))
} else {
lines(SMts[[f]]~year,col="gray",lwd=0.5)
}
}
lines(SMmean~year,col="black",lwd=2)
lines(SMub~year,col="black",lwd=2,lty=2)
lines(SMlb~year,col="black",lwd=2,lty=2)
lines(Gts3~year,col="purple",lwd=2)
lines(SMsub~year,col="purple",lwd=2,lty=2)
lines(SMslb~year,col="purple",lwd=2,lty=2)
legend("topleft",c("Ensemble members","Ens. Mean","Ens. Upper / Lower bounds","Derived Mean","Derived Upper / Lower Bounds"),col=c("gray","black","black","purple","purple"),lwd=c(0.5,2,2,2,2),lty=c(1,1,2,1,2))

####
# SD plot

plot(sqrt(Tts3)~SMsd,main=paste("Southwest Missouri Derived vs Ensemble SD \nfor ",titlevar,sep=""),pch=19,ylab="Derived SD",xlab="Ensemble SD",ylim=c(0,max(SMsd,sqrt(Tts3),na.rm=TRUE)),xlim=c(0,max(SMsd,sqrt(Tts3),na.rm=TRUE)))
abline(coef=c(0,1))

####
# contribution overlaid

SMsub1 = Gts3+2*sqrt(Vts3)
SMslb1 = Gts3-2*sqrt(Vts3)
SMsub2 = Gts3+2*(sqrt(Vts3+Mts3))
SMslb2 = Gts3-2*(sqrt(Vts3+Mts3))
SMsub3 = Gts3+2*(sqrt(Vts3+Mts3+Dts3))
SMslb3 = Gts3-2*(sqrt(Vts3+Mts3+Dts3))
SMsub4 = Gts3+2*(sqrt(Vts3+Mts3+Dts3+Sts3))
SMslb4 = Gts3-2*(sqrt(Vts3+Mts3+Dts3+Sts3))

for(f in 1:length(files)){
if(f==1){
plot(SMts[[f]]~year,col="gray",lwd=0.5,main=paste("Southwest Missouri Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SMtab,na.rm=TRUE),max(SMtab,na.rm=TRUE)))
} else {
lines(SMts[[f]]~year,col="gray",lwd=0.5)
}
}

lines(Gts3~year,col="purple",lwd=2)
lines(SMsub1~year,col="orange",lwd=2,lty=2)
lines(SMslb1~year,col="orange",lwd=2,lty=2)
lines(SMsub2~year,col="blue",lwd=2,lty=2)
lines(SMslb2~year,col="blue",lwd=2,lty=2)
lines(SMsub3~year,col="red",lwd=2,lty=2)
lines(SMslb3~year,col="red",lwd=2,lty=2)
lines(SMsub4~year,col="darkgreen",lwd=2,lty=2)
lines(SMslb4~year,col="darkgreen",lwd=2,lty=2)

legend("topleft",c("Ensemble members","Derived Mean","Natural Variability","GCM Model Uncertainty","Downscaling Uncertainty","Scenario Uncertainty"),col=c("gray","purple","orange","blue","red","darkgreen"),lwd=c(0.5,2,2,2,2,2),lty=c(1,2,2,2,2,2))



dev.off()




