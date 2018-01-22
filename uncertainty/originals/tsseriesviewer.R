#####################
# 
# domain average anomaly time series viewer

###
# Arguments set

varname = "tmax" # common name of variable

###
# load libraries and master list files
library(ncdf)
library(maps) # these two just to check plotting issues if need be
library(fields)

files = system(paste("ls uncertainty/regfits/",varname,"*.nc",sep=""),intern=TRUE)
files = files[-grep("PRISM",files)]

if(varname=="pr" & length(grep("pr25",files))>0){
files = files[-grep("pr25",files)]
}

if(varname=="tmax" & length(grep("tmax95",files))>0){
files = files[-grep("tmax95",files)]
}

if(varname=="tmin" & length(grep("tmin32",files))>0){
files = files[-grep("tmin32",files)]
}

####
# begin calcs

filenames = do.call(rbind,strsplit(files,"/"))

for(f in 1:length(files)){

filesplit = strsplit(filenames[f,3],"_")
GCM = filesplit[[1]][2]
DS = filesplit[[1]][3]
scen = filesplit[[1]][4]

test = open.ncdf(files[f])

if(f==1){
lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
}

#vardata = get.var.ncdf(test,paste(varname,"_anom",sep=""))
vardata = get.var.ncdf(test,"fittedvalues")

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

tseries = data.frame(matrix(vardata[17,12,],nrow=1,ncol=dim(vardata)[3]))
tmpframe = data.frame(scen,GCM,DS)
tmpframe = cbind(tmpframe,tseries)

if(f==1){
outputframe = tmpframe
} else {
outputframe = rbind(outputframe,tmpframe)
}

message("Finished file ",f," / ",length(files))

}

years=1950:2099

####
# scen var calc
outputframe$scenGCM = paste(outputframe$scen,outputframe$GCM,sep="-")
scenGCMavg = aggregate(outputframe[,4:153],by=list(scen=outputframe$scenGCM),sum,na.rm=TRUE)
scen = do.call(rbind,strsplit(as.character(scenGCMavg$scen),"-"))[,1]
scenGCMavg$scen = scen
scenavg = aggregate(scenGCMavg[,2:ncol(scenGCMavg)],by=list(scen=scenGCMavg$scen),mean,na.rm=TRUE)
S = apply(scenavg[,2:151],2,var,na.rm=TRUE)

#####
# GCM var calc
outputframe$scenDS = paste(outputframe$scen,outputframe$DS,sep="-")
scenDSvar = aggregate(outputframe[,4:153],by=list(scen=outputframe$scenDS),var,na.rm=TRUE)
scenDSvar$scen = do.call(rbind,strsplit(as.character(scenDSvar$scen),"-"))[,1]
DSsum = aggregate(scenDSvar[,2:151],by=list(scen=as.character(scenDSvar$scen)),sum,na.rm=TRUE)
M = apply(DSsum[,2:151],2,sum,na.rm=TRUE)*(1/7)*(1/6)

#####
# DS var calc
scenGCMvar = aggregate(outputframe[,4:153],by=list(scen=outputframe$scenGCM),var,na.rm=TRUE)
scenGCMvar$scen = do.call(rbind,strsplit(as.character(scenGCMvar$scen),"-"))[,1]
GCMsum = aggregate(scenGCMvar[,2:151],by=list(scen=as.character(scenGCMvar$scen)),sum,na.rm=TRUE)
D = apply(GCMsum[,2:151],2,sum,na.rm=TRUE)*(1/7)*(1/46)

#####
# plot values

plot(as.vector(S)~years,col="darkgreen",lwd=2,type="l",xlim=c(2000,2100))
lines(as.vector(M)~years,col="darkblue",lwd=2)
lines(as.vector(D)~years,col="red",lwd=2)



scenavg = aggregate(outputframe[,4:ncol(outputframe)],by=list(scen=outputframe$scen),mean,na.rm=TRUE)
GCMavg = aggregate(outputframe[,4:ncol(outputframe)],by=list(GCM=outputframe$GCM),mean,na.rm=TRUE)
DSavg = aggregate(outputframe[,4:ncol(outputframe)],by=list(DS=outputframe$DS),mean,na.rm=TRUE)

years=1950:2099

scens = as.matrix(scenavg[,2:ncol(scenavg)],nrow=nrow(scenavg),ncol=length(years))

plot(as.vector(scens[1,])~years,type="l",col="blue",ylim=c(-2,5.5),xlim=c(2000,2100))
lines(as.vector(scens[2,])~years,col="red")
lines(as.vector(scens[3,])~years,col="green")
lines(as.vector(scens[4,])~years,col="red",lty=2)
lines(as.vector(scens[5,])~years,col="darkblue",lty=2)
lines(as.vector(scens[6,])~years,col="darkred",lty=2)
lines(as.vector(scens[7,])~years,col="blue",lty=2)

scensd = apply(scens,2,sd,na.rm=TRUE)
plot(as.vector(scensd)~years,lwd=2,type="l",xlim=c(2000,2100))

GCMs = as.matrix(GCMavg[,2:ncol(GCMavg)],nrow=nrow(GCMavg),ncol=length(years))
GCMsd = apply(GCMs,2,sd,na.rm=TRUE)
lines(as.vector(GCMsd)~years,col="darkblue",lwd=2)

DSs = as.matrix(DSavg[,2:ncol(DSavg)],nrow=nrow(DSavg),ncol=length(years))
DSsd = apply(DSs,2,sd,na.rm=TRUE)
lines(as.vector(DSsd)~years,col="red",lwd=2)

















