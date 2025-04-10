setwd("pathto/data")

library(ncdf4)
library(maps)
library(mapdata)
library(maptools)
library(fields)
library(sp)
library(raster)
library(ggplot2)

distfunc = function(lon,lat,modelgrid){ # functions helps find nearest neighbors to desired location
  #model grid must have the lats and lons in the same format as the location lat and lon. modelgrid should also have r and c per grid point
  dist=3963 * acos(sin(lat/57.2958)*sin(modelgrid$lat/57.2958)+cos(lat/57.2958)*cos(modelgrid$lat/57.2958)*cos((modelgrid$lon/57.2958)-(lon/57.2958)))
  minidx = which(dist==min(dist,na.rm=TRUE))
  modelgrid[minidx,]
}

#####
# Version of this script to run parallel but also one domain as a time.

#statelist = c("louisiana","oklahoma","north carolina","florida","kentucky","SGP-NCA","SE-NCA")
enssize = 3
mask = "D"

#####
# masks

#for(s in 1:length(statelist)){
if(mask!="CONUS" & mask!="A" & mask!="B" & mask!="C" & mask!="D" & mask!="E" & mask!="F"){
  test = nc_open(paste("masks/",mask,"_mask.nc",sep=""))
  regionmask = ncvar_get(test,"mask")
  nc_close(test)
} else {
  if(mask=="A"){ loclon = -85.90625; loclat=31.78125;}
  if(mask=="B"){ loclon = -76.90625; loclat=41.65625;}
  if(mask=="C"){ loclon = -97.65625; loclat=32.78125;}
  if(mask=="D"){ loclon = -97.03125; loclat=42.28125;}
  if(mask=="E"){ loclon = -114.6562; loclat=33.28125;}
  if(mask=="F"){ loclon = -121.9062; loclat=45.78125;}
}

##########
#LOCA data pulls

histfile_LOCA = "pr_LOCA_climo_CONUS.nc"
projfile_LOCA = "pr_LOCA_climo_rcp85_CONUS.nc"

####
# pull historical data

nctest = nc_open(histfile_LOCA)
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")

varnamesprh = c()
prhist_LOCA = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamesprh[i]=nctest$var[[i]]$name
  prhist_LOCA[,,i]=ncvar_get(nctest,varnamesprh[i])
}

nc_close(nctest)

####
# pull projected data - pr

nctest = nc_open(projfile_LOCA)
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")

varnamesprp = c()
prproj_LOCA = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamesprp[i]=nctest$var[[i]]$name
  prproj_LOCA[,,i]=ncvar_get(nctest,varnamesprp[i])
}

nc_close(nctest)

##########
#GCM data pulls

histfile_GCM = "pr_GCM_climo_CONUS.nc"
projfile_GCM = "pr_GCM_climo_rcp85_CONUS.nc"

####
# pull historical data

nctest = nc_open(histfile_GCM)
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")

varnamesprh = c()
prhist_GCM = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamesprh[i]=nctest$var[[i]]$name
  prhist_GCM[,,i]=ncvar_get(nctest,varnamesprh[i])
}

nc_close(nctest)

####
# pull projected data - pr

nctest = nc_open(projfile_GCM)
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")

varnamesprp = c()
prproj_GCM = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamesprp[i]=nctest$var[[i]]$name
  prproj_GCM[,,i]=ncvar_get(nctest,varnamesprp[i])
}

nc_close(nctest)

######
# Calculate projected change

probs = prhist_LOCA[,,which(varnamesprh=="LIVNEH")]

varidx = which(varnamesprh!="LIVNEH" & varnamesprh!="ENS")
varnamesprh = varnamesprh[varidx]
prhist_LOCA = prhist_LOCA[,,varidx]
prhist_GCM = prhist_GCM[,,varidx]

varidx = which(varnamesprp!="LIVNEH" & varnamesprp!="ENS")
varnamesprp = varnamesprp[varidx]
prproj_LOCA = prproj_LOCA[,,varidx]
prproj_GCM = prproj_GCM[,,varidx]

prchange_LOCA = prproj_LOCA-prhist_LOCA
prchange_GCM = prproj_GCM-prhist_GCM

######
# mask everything

if(mask=="A" | mask=="B" | mask=="C" | mask=="D" | mask=="E" | mask=="F"){
  LON = rep(lon,each=length(lat))
  LAT = rep(lat,length(lon))
  R = rep(1:length(lon),each=length(lat))
  C = rep(1:length(lat),length(lon))
  modelgrid = data.frame(R,C,LON,LAT)
  names(modelgrid) = c("R","C","lon","lat")
  pointpick = distfunc(loclon,loclat,modelgrid)
  regionmask = matrix(NA,nrow=length(lon),ncol=length(lat))
  regionmask[pointpick$R[1],pointpick$C[1]]=1
}



for(i in 1:length(varnamesprh)){
  if(mask!="CONUS") {
    prchange_GCM[,,i] = ifelse(regionmask==1,prchange_GCM[,,i],NA)
    prhist_GCM[,,i] = ifelse(regionmask==1,prhist_GCM[,,i],NA)
    prproj_GCM[,,i] = ifelse(regionmask==1,prproj_GCM[,,i],NA)
    prchange_LOCA[,,i] = ifelse(regionmask==1,prchange_LOCA[,,i],NA)
    prhist_LOCA[,,i] = ifelse(regionmask==1,prhist_LOCA[,,i],NA)
    prproj_LOCA[,,i] = ifelse(regionmask==1,prproj_LOCA[,,i],NA)
  } 
}

if(mask!="CONUS") {
  probs = ifelse(regionmask==1,probs,NA)
}


####
# calculate means of all models

probsmean = mean(probs,na.rm=TRUE)

prhistmean_LOCA = apply(prhist_LOCA,3,mean,na.rm=TRUE)
prhistmean_GCM = apply(prhist_GCM,3,mean,na.rm=TRUE)

prprojmean_LOCA = apply(prproj_LOCA,3,mean,na.rm=TRUE)
prprojmean_GCM = apply(prproj_GCM,3,mean,na.rm=TRUE)

prchangemean_LOCA = apply(prchange_LOCA,3,mean,na.rm=TRUE)
prchangemean_GCM = apply(prchange_GCM,3,mean,na.rm=TRUE)

#####
# Determine percentiles - precip historical
tau <- seq(0,1,length=length(varnamesprh))
quant.target <- quantile(prhistmean_LOCA,tau,type=7,na.rm=TRUE)
perc.prhist_LOCA <- approx(quant.target,tau,prhistmean_LOCA,rule=2)$y

quant.target <- quantile(prhistmean_GCM,tau,type=7,na.rm=TRUE)
perc.prhist_GCM <- approx(quant.target,tau,prhistmean_GCM,rule=2)$y

prhisttab <- data.frame(varnamesprh,prhistmean_LOCA,prhistmean_GCM,perc.prhist_LOCA,perc.prhist_GCM)

#####
# Determine percentiles - precip projected
tau <- seq(0,1,length=length(varnamesprp))
quant.target <- quantile(prprojmean_LOCA,tau,type=7,na.rm=TRUE)
perc.prproj_LOCA <- approx(quant.target,tau,prprojmean_LOCA,rule=2)$y

quant.target <- quantile(prprojmean_GCM,tau,type=7,na.rm=TRUE)
perc.prproj_GCM <- approx(quant.target,tau,prprojmean_GCM,rule=2)$y

prprojtab <- data.frame(varnamesprp,prprojmean_LOCA,prprojmean_GCM,perc.prproj_LOCA,perc.prproj_GCM)


#####
# Determine percentiles - precip change
tau <- seq(0,1,length=length(varnamesprh))
quant.target <- quantile(prchangemean_LOCA,tau,type=7,na.rm=TRUE)
perc.prchange_LOCA <- approx(quant.target,tau,prchangemean_LOCA,rule=2)$y

quant.target <- quantile(prchangemean_GCM,tau,type=7,na.rm=TRUE)
perc.prchange_GCM <- approx(quant.target,tau,prchangemean_GCM,rule=2)$y

prchangetab <- data.frame(varnamesprh,prchangemean_LOCA,prchangemean_GCM,perc.prchange_LOCA,perc.prchange_GCM)

######
# RMSE calculation and plots for skill

combosused = combn(1:29,3,simplify=TRUE)
combonum = 1:ncol(combosused)
RMSEvals_GCM = RMSEvals_LOCA = c()

for(i in 1:ncol(combosused)){
  GCMvals = prhisttab[combosused[,i],3]
  LOCAvals = prhisttab[combosused[,i],2]
  
  RMSEvals_GCM[i] = sqrt(mean((GCMvals-probsmean)^2))
  RMSEvals_LOCA[i] = sqrt(mean((LOCAvals-probsmean)^2))
  message("Finished RMSE calcs for ",i," / ",ncol(combosused))
}

NRMSEvals_GCM = RMSEvals_GCM/sd(probs,na.rm=TRUE)
NRMSEvals_LOCA = RMSEvals_LOCA/sd(probs,na.rm=TRUE)

NRMSEtab = data.frame(combonum,RMSEvals_GCM,RMSEvals_LOCA,NRMSEvals_GCM,NRMSEvals_LOCA)
GCMminidx = which(NRMSEtab[,2]==min(NRMSEtab[,2])) # while NRMSE is calculated, RMSE version used for paper.
LOCAminidx = which(NRMSEtab[,3]==min(NRMSEtab[,3])) # switch to NRMSE makes no difference to GCM selection

## Plot when selecting from pre-downscaled ensemble, plotting historical period precipitation

pdf(paste("GCMselections_",mask,".pdf",sep=""),onefile=TRUE,width=18,height=6)
par(mfrow=c(1,3))

plot(prhistmean_GCM~perc.prhist_GCM,data=prhisttab,ylim=range(c(prhistmean_GCM,prhistmean_LOCA)),lwd=2,pch=19,ylab="Historical Precip (mm)",xlab="percentile",cex=2)
points(prhistmean_LOCA~perc.prhist_LOCA,data=prhisttab,lwd=2,col="blue",cex=2)
abline(h=probsmean,lty=3)

prhisttab_sub = prhisttab[combosused[,GCMminidx],]
mean(prhisttab_sub[,2])-probsmean
mean(prhisttab_sub[,3])-probsmean

GCMsin = as.character(prhisttab[combosused[,GCMminidx],1]) # grabs the names of selected GCMs

#abline(h=mean(prhisttab_sub[,2]),lty=2,col="blue")
#abline(h=mean(prhisttab_sub[,3]),lty=2,col="black")

points(prhistmean_GCM[1]~perc.prhist_GCM[1],data=prhisttab_sub,col="purple",lwd=2,pch=15,cex=2)
points(prhistmean_GCM[2]~perc.prhist_GCM[2],data=prhisttab_sub,col="green",lwd=2,pch=18,cex=2)
points(prhistmean_GCM[3]~perc.prhist_GCM[3],data=prhisttab_sub,col="red",lwd=2,pch=17,cex=2)

points(prhistmean_LOCA[1]~perc.prhist_LOCA[1],data=prhisttab_sub,col="purple",lwd=2,pch=0,cex=2)
points(prhistmean_LOCA[2]~perc.prhist_LOCA[2],data=prhisttab_sub,col="green",lwd=2,pch=5,cex=2)
points(prhistmean_LOCA[3]~perc.prhist_LOCA[3],data=prhisttab_sub,col="red",lwd=2,pch=2,cex=2)


legendtitles = c("CMIP5","LOCA",
                 paste(GCMsin[1],"CMIP5",sep=" "),paste(GCMsin[2],"CMIP5",sep=" "),paste(GCMsin[3],"CMIP5",sep=" "),
                 paste(GCMsin[1],"LOCA",sep=" "),paste(GCMsin[2],"LOCA",sep=" "),paste(GCMsin[3],"LOCA",sep=" "))

legend("topleft",legendtitles,
       pch=c(19,1,15,18,17,0,5,2),
       col=c("black","blue","purple","green","red","purple","green","red"),cex=1)


## Plot when selecting from pre-downscaled ensemble, plotting future period precipitation

plot(prprojmean_GCM~perc.prproj_GCM,data=prprojtab,ylim=range(c(prprojmean_GCM,prprojmean_LOCA)),lwd=2,pch=19,ylab="Future Precip (mm)",xlab="percentile",cex=2)
points(prprojmean_LOCA~perc.prproj_LOCA,data=prprojtab,lwd=2,col="blue",cex=2)

prprojtab_sub = prprojtab[combosused[,GCMminidx],]

#abline(h=mean(prhisttab_sub[,2]),lty=2,col="blue")
#abline(h=mean(prhisttab_sub[,3]),lty=2,col="black")

points(prprojmean_GCM[1]~perc.prproj_GCM[1],data=prprojtab_sub,col="purple",lwd=2,pch=15,cex=2)
points(prprojmean_GCM[2]~perc.prproj_GCM[2],data=prprojtab_sub,col="green",lwd=2,pch=18,cex=2)
points(prprojmean_GCM[3]~perc.prproj_GCM[3],data=prprojtab_sub,col="red",lwd=2,pch=17,cex=2)

points(prprojmean_LOCA[1]~perc.prproj_LOCA[1],data=prprojtab_sub,col="purple",lwd=2,pch=0,cex=2)
points(prprojmean_LOCA[2]~perc.prproj_LOCA[2],data=prprojtab_sub,col="green",lwd=2,pch=5,cex=2)
points(prprojmean_LOCA[3]~perc.prproj_LOCA[3],data=prprojtab_sub,col="red",lwd=2,pch=2,cex=2)

legendtitles = c("CMIP5","LOCA",
                 paste(GCMsin[1],"CMIP5",sep=" "),paste(GCMsin[2],"CMIP5",sep=" "),paste(GCMsin[3],"CMIP5",sep=" "),
                 paste(GCMsin[1],"LOCA",sep=" "),paste(GCMsin[2],"LOCA",sep=" "),paste(GCMsin[3],"LOCA",sep=" "))

legend("topleft",legendtitles,
       pch=c(19,1,15,18,17,0,5,2),
       col=c("black","blue","purple","green","red","purple","green","red"),cex=1)

## Plot when selecting from pre-downscaled ensemble, plotting projected change

plot(prchangemean_GCM~perc.prchange_GCM,data=prchangetab,ylim=range(c(prchangemean_GCM,prchangemean_LOCA)),lwd=2,pch=19,ylab="Projected Change (mm)",xlab="percentile",cex=2)
points(prchangemean_LOCA~perc.prchange_LOCA,data=prchangetab,lwd=2,col="blue",cex=2)
abline(h=0,lty=3)

prchangetab_sub = prchangetab[combosused[,GCMminidx],]

#abline(h=mean(prhisttab_sub[,2]),lty=2,col="blue")
#abline(h=mean(prhisttab_sub[,3]),lty=2,col="black")

points(prchangemean_GCM[1]~perc.prchange_GCM[1],data=prchangetab_sub,col="purple",lwd=2,pch=15,cex=2)
points(prchangemean_GCM[2]~perc.prchange_GCM[2],data=prchangetab_sub,col="green",lwd=2,pch=18,cex=2)
points(prchangemean_GCM[3]~perc.prchange_GCM[3],data=prchangetab_sub,col="red",lwd=2,pch=17,cex=2)

points(prchangemean_LOCA[1]~perc.prchange_LOCA[1],data=prchangetab_sub,col="purple",lwd=2,pch=0,cex=2)
points(prchangemean_LOCA[2]~perc.prchange_LOCA[2],data=prchangetab_sub,col="green",lwd=2,pch=5,cex=2)
points(prchangemean_LOCA[3]~perc.prchange_LOCA[3],data=prchangetab_sub,col="red",lwd=2,pch=2,cex=2)

legendtitles = c("CMIP5","LOCA",
                 paste(GCMsin[1],"CMIP5",sep=" "),paste(GCMsin[2],"CMIP5",sep=" "),paste(GCMsin[3],"CMIP5",sep=" "),
                 paste(GCMsin[1],"LOCA",sep=" "),paste(GCMsin[2],"LOCA",sep=" "),paste(GCMsin[3],"LOCA",sep=" "))

legend("topleft",legendtitles,
       pch=c(19,1,15,18,17,0,5,2),
       col=c("black","blue","purple","green","red","purple","green","red"),cex=1)

dev.off()


#############
# repeat plots for post-downscaled ensemble figures

## Plot when selecting from post-downscaled ensemble, plotting historical period precipitation

pdf(paste("LOCAselections_",mask,".pdf",sep=""),onefile=TRUE,width=18,height=6)
par(mfrow=c(1,3))

plot(prhistmean_GCM~perc.prhist_GCM,data=prhisttab,ylim=range(c(prhistmean_GCM,prhistmean_LOCA)),lwd=2,pch=19,ylab="Historical Precip (mm)",xlab="percentile",cex=2)
points(prhistmean_LOCA~perc.prhist_LOCA,data=prhisttab,lwd=2,col="blue",cex=2)
abline(h=probsmean,lty=3)

prhisttab_sub = prhisttab[combosused[,LOCAminidx],] # switch from GCMminidx to LOCAminidx to reverse the experiment
mean(prhisttab_sub[,2])-probsmean
mean(prhisttab_sub[,3])-probsmean

GCMsin = as.character(prhisttab[combosused[,LOCAminidx],1]) # grabs the names of selected GCMs

#abline(h=mean(prhisttab_sub[,2]),lty=2,col="blue")
#abline(h=mean(prhisttab_sub[,3]),lty=2,col="black")

points(prhistmean_GCM[1]~perc.prhist_GCM[1],data=prhisttab_sub,col="purple",lwd=2,pch=15,cex=2)
points(prhistmean_GCM[2]~perc.prhist_GCM[2],data=prhisttab_sub,col="green",lwd=2,pch=18,cex=2)
points(prhistmean_GCM[3]~perc.prhist_GCM[3],data=prhisttab_sub,col="red",lwd=2,pch=17,cex=2)

points(prhistmean_LOCA[1]~perc.prhist_LOCA[1],data=prhisttab_sub,col="purple",lwd=2,pch=0,cex=2)
points(prhistmean_LOCA[2]~perc.prhist_LOCA[2],data=prhisttab_sub,col="green",lwd=2,pch=5,cex=2)
points(prhistmean_LOCA[3]~perc.prhist_LOCA[3],data=prhisttab_sub,col="red",lwd=2,pch=2,cex=2)

legendtitles = c("CMIP5","LOCA",
                 paste(GCMsin[1],"CMIP5",sep=" "),paste(GCMsin[2],"CMIP5",sep=" "),paste(GCMsin[3],"CMIP5",sep=" "),
                 paste(GCMsin[1],"LOCA",sep=" "),paste(GCMsin[2],"LOCA",sep=" "),paste(GCMsin[3],"LOCA",sep=" "))

legend("topleft",legendtitles,
       pch=c(19,1,15,18,17,0,5,2),
       col=c("black","blue","purple","green","red","purple","green","red"),cex=1)

## Plot when selecting from post-downscaled ensemble, plotting future period precipitation

plot(prprojmean_GCM~perc.prproj_GCM,data=prprojtab,ylim=range(c(prprojmean_GCM,prprojmean_LOCA)),lwd=2,pch=19,ylab="Future Precip (mm)",xlab="percentile",cex=2)
points(prprojmean_LOCA~perc.prproj_LOCA,data=prprojtab,lwd=2,col="blue",cex=2)

prprojtab_sub = prprojtab[combosused[,LOCAminidx],]

#abline(h=mean(prhisttab_sub[,2]),lty=2,col="blue")
#abline(h=mean(prhisttab_sub[,3]),lty=2,col="black")

points(prprojmean_GCM[1]~perc.prproj_GCM[1],data=prprojtab_sub,col="purple",lwd=2,pch=15,cex=2)
points(prprojmean_GCM[2]~perc.prproj_GCM[2],data=prprojtab_sub,col="green",lwd=2,pch=18,cex=2)
points(prprojmean_GCM[3]~perc.prproj_GCM[3],data=prprojtab_sub,col="red",lwd=2,pch=17,cex=2)

points(prprojmean_LOCA[1]~perc.prproj_LOCA[1],data=prprojtab_sub,col="purple",lwd=2,pch=0,cex=2)
points(prprojmean_LOCA[2]~perc.prproj_LOCA[2],data=prprojtab_sub,col="green",lwd=2,pch=5,cex=2)
points(prprojmean_LOCA[3]~perc.prproj_LOCA[3],data=prprojtab_sub,col="red",lwd=2,pch=2,cex=2)

legendtitles = c("CMIP5","LOCA",
                 paste(GCMsin[1],"CMIP5",sep=" "),paste(GCMsin[2],"CMIP5",sep=" "),paste(GCMsin[3],"CMIP5",sep=" "),
                 paste(GCMsin[1],"LOCA",sep=" "),paste(GCMsin[2],"LOCA",sep=" "),paste(GCMsin[3],"LOCA",sep=" "))

legend("topleft",legendtitles,
       pch=c(19,1,15,18,17,0,5,2),
       col=c("black","blue","purple","green","red","purple","green","red"),cex=1)

## Plot when selecting from post-downscaled ensemble, plotting projected change

plot(prchangemean_GCM~perc.prchange_GCM,data=prchangetab,ylim=range(c(prchangemean_GCM,prchangemean_LOCA)),lwd=2,pch=19,ylab="Projected Change (mm)",xlab="percentile",cex=2)
points(prchangemean_LOCA~perc.prchange_LOCA,data=prchangetab,lwd=2,col="blue",cex=2)
abline(h=0,lty=3)

prchangetab_sub = prchangetab[combosused[,LOCAminidx],]

#abline(h=mean(prhisttab_sub[,2]),lty=2,col="blue")
#abline(h=mean(prhisttab_sub[,3]),lty=2,col="black")

points(prchangemean_GCM[1]~perc.prchange_GCM[1],data=prchangetab_sub,col="purple",lwd=2,pch=15,cex=2)
points(prchangemean_GCM[2]~perc.prchange_GCM[2],data=prchangetab_sub,col="green",lwd=2,pch=18,cex=2)
points(prchangemean_GCM[3]~perc.prchange_GCM[3],data=prchangetab_sub,col="red",lwd=2,pch=17,cex=2)

points(prchangemean_LOCA[1]~perc.prchange_LOCA[1],data=prchangetab_sub,col="purple",lwd=2,pch=0,cex=2)
points(prchangemean_LOCA[2]~perc.prchange_LOCA[2],data=prchangetab_sub,col="green",lwd=2,pch=5,cex=2)
points(prchangemean_LOCA[3]~perc.prchange_LOCA[3],data=prchangetab_sub,col="red",lwd=2,pch=2,cex=2)

legendtitles = c("CMIP5","LOCA",
                 paste(GCMsin[1],"CMIP5",sep=" "),paste(GCMsin[2],"CMIP5",sep=" "),paste(GCMsin[3],"CMIP5",sep=" "),
                 paste(GCMsin[1],"LOCA",sep=" "),paste(GCMsin[2],"LOCA",sep=" "),paste(GCMsin[3],"LOCA",sep=" "))

legend("topleft",legendtitles,
       pch=c(19,1,15,18,17,0,5,2),
       col=c("black","blue","purple","green","red","purple","green","red"),cex=1)

dev.off()
