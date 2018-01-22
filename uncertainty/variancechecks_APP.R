################
#
# Appalachians Analysis for Variance Check

############
# library and function load

############
# load libraries

###############3
# 
# Variance Check

library(sp)
library(fields)
library(raster)
library(rasterVis)
library(maptools)
library(maps)
library(ncdf)

############
# load calculator functions and maps needed
source("componentfunctions.R")

test = readShapePoly("physio_shp/physio")
projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

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
# Get mask of points for the Appalachians

test.sub <- test[as.character(test@data$FENCODE) %in% c("5a","5b","6b","8d"), ]
  #test.sub <- test[as.character(test@data$FENCODE) %in% c("5b"), ]
tmpV = V

modrasV = raster(t(tmpV)[length(lat):1,])
extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
mod.subV <- crop(modrasV, extent(test.sub))
mod.subV <- mask(modrasV, test.sub)

testin = matrix(getValues(mod.subV),nrow=length(lon),ncol=length(lat))
testin = testin[,length(lat):1]

maskuse = array(NA,dim=c(length(lon),length(lat),94))
for(i in 1:94) maskuse[,,i] = ifelse(is.na(testin)==FALSE,1,0)


################
# Get domain tseries

V = ifelse(maskuse==1,V,NA)
M = ifelse(maskuse==1,M,NA)
D = ifelse(maskuse==1,D,NA)
S = ifelse(maskuse==1,S,NA)
G = ifelse(maskuse==1,G,NA)
T = ifelse(maskuse==1,T,NA)

Vts1 = mean(V,na.rm=TRUE) 
Mts1 = apply(M,3,mean,na.rm=TRUE)
Dts1 = apply(D,3,mean,na.rm=TRUE)
Sts1 = apply(S,3,mean,na.rm=TRUE)
Gts1 = apply(G,3,mean,na.rm=TRUE)
Tts1 = apply(T,3,mean,na.rm=TRUE)

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

SAts = list()

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

vardata = vardata[,,timeidx]
vardata = ifelse(maskuse==1,vardata,NA)

SAts[[f]]=apply(vardata,3,mean,na.rm=TRUE)

message("Finished time series grabs for file ",f," / ",length(files))

}

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

year =  2006:2099

pdf(paste("uncertainty/SE/analysisid/",varname,"_variancecheck_fix_APP.pdf",sep=""),height=6,width=6,onefile=TRUE)

###
# Southern Appalachians

####
# ens. sd overlaid with estimated sd

SAtab = do.call(rbind,SAts)
SAmean = apply(SAtab,2,mean,na.rm=TRUE)
SAsd = apply(SAtab,2,sd,na.rm=TRUE)
SAub = SAmean+2*SAsd
SAlb = SAmean-2*SAsd
SAsub = Gts1+2*sqrt(Tts1)
SAslb = Gts1-2*sqrt(Tts1)

for(f in 1:length(files)){
if(f==1){
plot(SAts[[f]]~year,col="gray",lwd=0.5,main=paste("Appalachians Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SAtab,na.rm=TRUE),max(SAtab,na.rm=TRUE)))
} else {
lines(SAts[[f]]~year,col="gray",lwd=0.5)
}
}
lines(SAmean~year,col="black",lwd=2)
lines(SAub~year,col="black",lwd=2,lty=2)
lines(SAlb~year,col="black",lwd=2,lty=2)
lines(Gts1~year,col="purple",lwd=2)
lines(SAsub~year,col="purple",lwd=2,lty=2)
lines(SAslb~year,col="purple",lwd=2,lty=2)
legend("topleft",c("Ensemble members","Ens. Mean","Ens. Upper / Lower bounds","Derived Mean","Derived Upper / Lower Bounds"),col=c("gray","black","black","purple","purple"),lwd=c(0.5,2,2,2,2),lty=c(1,1,2,1,2))

####
# SD plot

plot(sqrt(Tts1)~SAsd,main=paste("Appalachians Derived vs Ensemble SD \nfor ",titlevar,sep=""),pch=19,ylab="Derived SD",xlab="Ensemble SD",ylim=c(0,max(SAsd,sqrt(Tts2),na.rm=TRUE)),xlim=c(0,max(SAsd,sqrt(Tts2),na.rm=TRUE)))
abline(coef=c(0,1))

####
# contribution overlaid

SAsub1 = Gts1+2*sqrt(Vts1)
SAslb1 = Gts1-2*sqrt(Vts1)
SAsub2 = Gts1+2*(sqrt(Vts1+Mts1))
SAslb2 = Gts1-2*(sqrt(Vts1+Mts1))
SAsub3 = Gts1+2*(sqrt(Vts1+Mts1+Dts1))
SAslb3 = Gts1-2*(sqrt(Vts1+Mts1+Dts1))
SAsub4 = Gts1+2*(sqrt(Vts1+Mts1+Dts1+Sts1))
SAslb4 = Gts1-2*(sqrt(Vts1+Mts1+Dts1+Sts1))

for(f in 1:length(files)){
if(f==1){
plot(SAts[[f]]~year,col="gray",lwd=0.5,main=paste("Appalachians Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SAtab,na.rm=TRUE),max(SAtab,na.rm=TRUE)))
} else {
lines(SAts[[f]]~year,col="gray",lwd=0.5)
}
}

lines(Gts1~year,col="purple",lwd=2,lty=2)
lines(SAsub1~year,col="orange",lwd=2,lty=2)
lines(SAslb1~year,col="orange",lwd=2,lty=2)
lines(SAsub2~year,col="blue",lwd=2,lty=2)
lines(SAslb2~year,col="blue",lwd=2,lty=2)
lines(SAsub3~year,col="red",lwd=2,lty=2)
lines(SAslb3~year,col="red",lwd=2,lty=2)
lines(SAsub4~year,col="darkgreen",lwd=2,lty=2)
lines(SAslb4~year,col="darkgreen",lwd=2,lty=2)

legend("topleft",c("Ensemble members","Derived Mean","Natural Variability","GCM Model Uncertainty","Downscaling Uncertainty","Scenario Uncertainty"),col=c("gray","purple","orange","blue","red","darkgreen"),lwd=c(0.5,2,2,2,2,2),lty=c(1,2,2,2,2,2))

dev.off()




