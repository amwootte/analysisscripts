################
#
# Spaghetti Plotter

############
# library and function load

############
# load libraries

library(ncdf)
library(fields)
library(maps)

############
# load calculator functions
source("componentfunctions.R")

###
# set variable name
varname = "tmax"
anomfilepath = "anomdat"
anfilepath = "analysis"

###
# set extra options
printPDF = TRUE # if you want some plots printed to PDF set this to TRUE

###################
###################
# Get anomaly values

files = system(paste("ls uncertainty/SE/",anomfilepath,"/",varname,"*.nc",sep=""),intern=TRUE)

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

lonidx1 = which(lon>= -83.5 & lon<= -80)
latidx1 = which(lat>= 27 & lat>= 29)

lonidx2 = which(lon>= -85.5 & lon<= -81)
latidx2 = which(lat>= 34.5 & lat>= 34)

lonidx3 = which(lon>= -95.5 & lon<= -93)
latidx3 = which(lat>= 35.5 & lat>= 39)
}

vardata = get.var.ncdf(test,paste(varname,"_anom",sep=""))

spunits = test$var[[3]]$units

time = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
close.ncdf(test)

if(datasetname=="maca" |  datasetname=="cmip5.bcca"){
#vardata2 = vardata
vardata[,1,] = NA
vardata[1,,]=NA
} 

if(datasetname=="dcp"){
mask = ifelse(is.na(vardata[,,100])==TRUE,0,1)
for(i in 1:56){
vardata[,,i] = ifelse(mask==0,NA,vardata[,,i])
}
vardata[,,11] = NA
}

if(datasetname=="ccr" & (varname=="tmin" | varname=="tmax")){
vardata = ifelse(vardata < -5,NA,vardata)
}

CFts[[f]]=apply(vardata[lonidx1,latidx1,timeidx],3,mean,na.rm=TRUE)
SAts[[f]]=apply(vardata[lonidx2,latidx2,timeidx],3,mean,na.rm=TRUE)
SMts[[f]]=apply(vardata[lonidx3,latidx3,timeidx],3,mean,na.rm=TRUE)

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

pdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_spaghettiplot_fix.pdf",sep=""),height=6,width=6,onefile=TRUE)

year = 2006:2099

###
# Central Florida

####
# ens. sd overlaid with estimated sd

CFtab = do.call(rbind,CFts)
CFtab[38,1:ncol(CFtab)] = rep(NA,ncol(CFtab))
CFtab[193,9] = NA
CFmean = apply(CFtab,2,mean,na.rm=TRUE)


for(f in 1:length(files)){
if(f==1){
plot(CFts[[f]]~year,col="gray",lwd=0.5,main=paste("Central Florida Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(CFtab,na.rm=TRUE),max(CFtab,na.rm=TRUE)))
} else {
lines(CFts[[f]]~year,col="gray",lwd=0.5)
}
}
lines(CFmean~year,col="black",lwd=2)
legend("topleft",c("Individual Members","Ens. Mean"),col=c("gray","black"),lwd=c(0.5,2),lty=c(1,1))


###
# Southern Appalachians

####
# ens. sd overlaid with estimated sd

SAtab = do.call(rbind,SAts)
SAtab[38,1:ncol(SAtab)] = rep(NA,ncol(SAtab))
SAtab[193,9] = NA
SAmean = apply(SAtab,2,mean,na.rm=TRUE)

for(f in 1:length(files)){
if(f==1){
plot(SAts[[f]]~year,col="gray",lwd=0.5,main=paste("Southern Appalachians Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SAtab,na.rm=TRUE),max(SAtab,na.rm=TRUE)))
} else {
lines(SAts[[f]]~year,col="gray",lwd=0.5)
}
}
lines(SAmean~year,col="black",lwd=2)
legend("topleft",c("Individual members","Ens. Mean"),col=c("gray","black"),lwd=c(0.5,2),lty=c(1,1))

###
# Southwest Missouri

####
# ens. sd overlaid with estimated sd

SMtab = do.call(rbind,SMts)
SMtab[38,1:ncol(SMtab)] = rep(NA,ncol(SMtab))
SMtab[193,9] = NA
SMmean = apply(SMtab,2,mean,na.rm=TRUE)

for(f in 1:length(files)){
if(f==1){
plot(SMts[[f]]~year,col="gray",lwd=0.5,main=paste("Southwest Missouri Spaghetti Plot \nfor ",titlevar,sep=""),ylab=paste(varname," anomaly ", titlevarunits,sep=""),xlab="Year",type="l",ylim=c(min(SMtab,na.rm=TRUE),max(SMtab,na.rm=TRUE)))
} else {
lines(SMts[[f]]~year,col="gray",lwd=0.5)
}
}
lines(SMmean~year,col="black",lwd=2)
legend("topleft",c("Individual members","Ens. Mean"),col=c("gray","black"),lwd=c(0.5,2),lty=c(1,1))

dev.off()




