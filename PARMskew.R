source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

setwd("/home/woot0002/")

PARMdtrfile = "dtr_day_I35dtrdetrp1-PARM-B10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
PARMtmaxfile = "tasmax_day_I35txdetrp1-PARM-B10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
PARMtminfile = "tasmin_day_I35tndetrp1-PARM-Bd10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"

EDQMtminfile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc" 
EDQMtmaxfile = "/data2/3to5/I35/tasmax/EDQM/tasmax_day_I35txp1-EDQM-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"

PRISMtminfile = "/data2/3to5/I35/tasmin/PRISM/tasmin_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
PRISMtmaxfile = "/data2/3to5/I35/tasmax/PRISM/tasmax_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

test = nc_open(PARMdtrfile)
PARMdtr = ncvar_get(test,"dtr")
lat = ncvar_get(test,"lat")
lon = ncvar_get(test,"lon")
nc_close(test)

test = nc_open(PARMtmaxfile)
PARMtmax = ncvar_get(test,"tasmax")
nc_close(test)

test = nc_open(PARMtminfile)
PARMtmin = ncvar_get(test,"tasmin")
nc_close(test)

histdates= seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

if(length(histdates)>dim(PARMdtr)[3]){
  histdates2 = histdates[-which(substr(histdates,6,10)=="02-29")]
}

test = nc_open(EDQMtmaxfile)
EDQMtmax = ncvar_get(test,"tasmax")
nc_close(test)

test = nc_open(EDQMtminfile)
EDQMtmin = ncvar_get(test,"tasmin")
nc_close(test)

EDQMdtr = EDQMtmax-EDQMtmin

test = nc_open(PRISMtmaxfile)
PRISMtmax = ncvar_get(test,"tasmax")
nc_close(test)

test = nc_open(PRISMtminfile)
PRISMtmin = ncvar_get(test,"tasmin")
nc_close(test)

PRISMdtr = PRISMtmax-PRISMtmin

PARMdtrcheck = PARMtmax - PARMtmin

PARMdtr2 = ifelse(PARMdtr<0.01,0.01,PARMdtr)
PARMtmincheck = PARMtmax - PARMdtr2

range(PARMdtrcheck-PARMdtr2,na.rm=TRUE)



library(e1071)

PARMdtrskew = PARMdtr2skew = EDQMdtrskew = PRISMdtrskew = PARMtmaxskew = EDQMtmaxskew = PRISMtmaxskew = PARMtminskew = EDQMtminskew = PRISMtminskew = matrix(NA,nrow=length(lon),ncol=length(lat))

for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    
    if(is.na(PARMtmax[r,c,1])==FALSE){
      PARMdtrskew[r,c] = skewness(PARMdtr[r,c,],na.rm=TRUE)
      PARMdtr2skew[r,c] = skewness(PARMdtr2[r,c,],na.rm=TRUE)
      EDQMdtrskew[r,c] = skewness(EDQMdtr[r,c,],na.rm=TRUE)
      PRISMdtrskew[r,c] = skewness(PRISMdtr[r,c,],na.rm=TRUE)
      
      PARMtmaxskew[r,c] = skewness(PARMtmax[r,c,],na.rm=TRUE)
      EDQMtmaxskew[r,c] = skewness(EDQMtmax[r,c,],na.rm=TRUE)
      PRISMtmaxskew[r,c] = skewness(PRISMtmax[r,c,],na.rm=TRUE)
      
      PARMtminskew[r,c] = skewness(PARMtmin[r,c,],na.rm=TRUE)
      EDQMtminskew[r,c] = skewness(EDQMtmin[r,c,],na.rm=TRUE)
      PRISMtminskew[r,c] = skewness(PRISMtmin[r,c,],na.rm=TRUE)
      
    }
  }
  message("Finished calcs for row ",r," / ",length(lon))
}

######

range(c(PARMdtrskew,PARMdtr2skew,EDQMdtrskew,PRISMdtrskew,PARMtminskew,EDQMtminskew,PRISMtminskew,PARMtmaxskew,EDQMtmaxskew,PRISMtmaxskew),na.rm=TRUE)

rawcolorbar = colorramp(c(PARMdtrskew,PARMdtr2skew,EDQMdtrskew,PRISMdtrskew,PARMtminskew,EDQMtminskew,PRISMtminskew,PARMtmaxskew,EDQMtmaxskew,PRISMtmaxskew),colorchoice="bluetored",type="difference",Blimit=50,use_fixed_scale = TRUE, fixed_scale = c(-1.3,1.3))

#testsfc1 = list(x=(lon-360),y=lat,z=PARMdtrskew)
#surface(testsfc1,type="I",main="PARM dtr skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
#map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMdtr2skew)
surface(testsfc1,type="I",main="PARM dtr (post-processed) skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMdtrskew)
surface(testsfc1,type="I",main="EDQM dtr skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMdtrskew)
surface(testsfc1,type="I",main="PRISM dtr skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMtmaxskew)
surface(testsfc1,type="I",main="PARM tmax skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMtmaxskew)
surface(testsfc1,type="I",main="EDQM tmax skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMtmaxskew)
surface(testsfc1,type="I",main="PRISM tmax skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMtminskew)
surface(testsfc1,type="I",main="PARM tmin skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMtminskew)
surface(testsfc1,type="I",main="EDQM tmin skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMtminskew)
surface(testsfc1,type="I",main="PRISM tmin skew",zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)



mask = ifelse(PARMtminskew<0 & PARMtmaxskew<0,1,0)
PARMskewdiff = ifelse(mask==1,PARMtminskew-PARMtmaxskew,NA)

testsfc1 = list(x=(lon-360),y=lat,z=mask)
surface(testsfc1,type="I",main="PARM skew mask")
map("state",add=TRUE)

diffcolorbar = colorramp(c(PARMskewdiff),colorchoice="bluetored",type="difference",Blimit=50,use_fixed_scale = TRUE,fixed_scale = c(-0.55,0.3))

testsfc1 = list(x=(lon-360),y=lat,z=PARMskewdiff)
surface(testsfc1,type="I",main="PARM skew difference (tmin-tmax)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

which(PARMskewdiff==min(PARMskewdiff,na.rm=TRUE),arr.ind = TRUE)

nrow(which(PARMskewdiff<=0,arr.ind = TRUE))
nrow(which(is.na(PARMskewdiff)==FALSE,arr.ind = TRUE))

i=5
j=117

PARMdtr2vals = PARMdtr2[i,j,] 
PARMtmaxvals = PARMtmax[i,j,]
PARMtminvals = PARMtmin[i,j,]
PRISMtmaxvals = PRISMtmax[i,j,]
PRISMtminvals = PRISMtmin[i,j,]
PRISMdtrvals = PRISMdtr[i,j,]
EDQMtmaxvals = EDQMtmax[i,j,]
EDQMtminvals = EDQMtmin[i,j,]
EDQMdtrvals = EDQMdtr[i,j,]


PARMtmaxquant = PARMtminquant = PARMdtrquant = c()
EDQMtmaxquant = EDQMtminquant = EDQMdtrquant = c()
PRISMtmaxquant = PRISMtminquant = PRISMdtrquant = c()
probsin = seq(0.01,0.99,by=0.01)

for(p in 1:length(probsin)){
  
  PARMtmaxquant[p] = quantile(PARMtmaxvals,probs=probsin[p],na.rm=TRUE) 
  PARMtminquant[p] = quantile(PARMtminvals,probs=probsin[p],na.rm=TRUE) 
  PARMdtrquant[p] = quantile(PARMdtr2vals,probs=probsin[p],na.rm=TRUE) 
  
  EDQMtmaxquant[p] = quantile(EDQMtmaxvals,probs=probsin[p],na.rm=TRUE) 
  EDQMtminquant[p] = quantile(EDQMtminvals,probs=probsin[p],na.rm=TRUE) 
  EDQMdtrquant[p] = quantile(EDQMdtrvals,probs=probsin[p],na.rm=TRUE) 
  
  PRISMtmaxquant[p] = quantile(PRISMtmaxvals,probs=probsin[p],na.rm=TRUE) 
  PRISMtminquant[p] = quantile(PRISMtminvals,probs=probsin[p],na.rm=TRUE) 
  PRISMdtrquant[p] = quantile(PRISMdtrvals,probs=probsin[p],na.rm=TRUE) 
  
}

PARMtmaxquantanom = PARMtmaxquant-mean(PARMtmaxquant)
PARMtminquantanom = PARMtminquant-mean(PARMtminquant)
PARMdtrquantanom = PARMdtrquant-mean(PARMdtrquant)

plot(probsin~PARMtmaxquantanom,xlim=range(c(PARMtmaxquantanom,PARMtminquantanom,PARMdtrquantanom)),lwd=2,type="l",col="red")
lines(probsin~PARMtminquantanom,lwd=2,col="blue")
lines(probsin~PARMdtrquantanom,lwd=2,col="black")

plot(probsin~PARMdtrquant,lwd=2,type="l")

##########

PARMdtr2skew = EDQMdtrskew = PRISMdtrskew = PARMtmaxskew = EDQMtmaxskew = PRISMtmaxskew = PARMtminskew = EDQMtminskew = PRISMtminskew = array(NA,dim=c(length(lon),length(lat),12))

for(m in 1:12){
  
  monidx = which(as.numeric(substr(histdates2,6,7))==m)

for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    
    if(is.na(PARMtmax[r,c,1])==FALSE){
      PARMdtr2skew[r,c,m] = skewness(PARMdtr2[r,c,monidx],na.rm=TRUE)
      EDQMdtrskew[r,c,m] = skewness(EDQMdtr[r,c,monidx],na.rm=TRUE)
      PRISMdtrskew[r,c,m] = skewness(PRISMdtr[r,c,monidx],na.rm=TRUE)
      
      PARMtmaxskew[r,c,m] = skewness(PARMtmax[r,c,monidx],na.rm=TRUE)
      EDQMtmaxskew[r,c,m] = skewness(EDQMtmax[r,c,monidx],na.rm=TRUE)
      PRISMtmaxskew[r,c,m] = skewness(PRISMtmax[r,c,monidx],na.rm=TRUE)
      
      PARMtminskew[r,c,m] = skewness(PARMtmin[r,c,monidx],na.rm=TRUE)
      EDQMtminskew[r,c,m] = skewness(EDQMtmin[r,c,monidx],na.rm=TRUE)
      PRISMtminskew[r,c,m] = skewness(PRISMtmin[r,c,monidx],na.rm=TRUE)
      
    }
  }
  message("Finished calcs for row ",r," / ",length(lon))
}
  message("Finished calcs for month ",m)
}

####

rawcolorbar = colorramp(c(PARMdtr2skew,EDQMdtrskew,PRISMdtrskew,PARMtminskew,EDQMtminskew,PRISMtminskew,PARMtmaxskew,EDQMtmaxskew,PRISMtmaxskew),colorchoice="bluetored",type="difference",Blimit=30)

pdf("monthlyskewness.pdf",onefile=TRUE,width=15,height=5)
par(mfrow=c(1,3))

for(m in 1:12){

testsfc1 = list(x=(lon-360),y=lat,z=PARMdtr2skew[,,m])
surface(testsfc1,type="I",main=paste("PARM dtr (post-processed) skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMdtrskew[,,m])
surface(testsfc1,type="I",main=paste("EDQM dtr skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMdtrskew[,,m])
surface(testsfc1,type="I",main=paste("PRISM dtr skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMtmaxskew[,,m])
surface(testsfc1,type="I",main=paste("PARM tmax skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMtmaxskew[,,m])
surface(testsfc1,type="I",main=paste("EDQM tmax skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMtmaxskew[,,m])
surface(testsfc1,type="I",main=paste("PRISM tmax skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMtminskew[,,m])
surface(testsfc1,type="I",main=paste("PARM tmin skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMtminskew[,,m])
surface(testsfc1,type="I",main=paste("EDQM tmin skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMtminskew[,,m])
surface(testsfc1,type="I",main=paste("PRISM tmin skew for month: ",m,sep=""),zlim=rawcolorbar[[1]],col=rawcolorbar[[3]],breaks=rawcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

}

dev.off()

j=5
i=110

PARMdtr2skew[i,j,]
PARMtmaxskew[i,j,]
PARMtminskew[i,j,]

probsin = seq(0.01,0.99,by=0.01)
PARMtmaxquant = PARMtminquant = PARMdtrquant = matrix(NA,nrow=length(probsin),ncol=12)
EDQMtmaxquant = EDQMtminquant = EDQMdtrquant = matrix(NA,nrow=length(probsin),ncol=12)
PRISMtmaxquant = PRISMtminquant = PRISMdtrquant = matrix(NA,nrow=length(probsin),ncol=12)

for(m in 1:12){
  monidx = which(as.numeric(substr(histdates2,6,7))==m)
for(p in 1:length(probsin)){
  
  PARMtmaxquant[p,m] = quantile(PARMtmax[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  PARMtminquant[p,m] = quantile(PARMtmin[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  PARMdtrquant[p,m] = quantile(PARMdtr2[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  
  EDQMtmaxquant[p,m] = quantile(EDQMtmax[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  EDQMtminquant[p,m] = quantile(EDQMtmin[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  EDQMdtrquant[p,m] = quantile(EDQMdtr[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  
  PRISMtmaxquant[p,m] = quantile(PRISMtmax[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  PRISMtminquant[p,m] = quantile(PRISMtmin[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  PRISMdtrquant[p,m] = quantile(PRISMdtr[i,j,monidx],probs=probsin[p],na.rm=TRUE) 
  
}

}

plot(probsin~PARMtmaxquant[,5],xlim=range(c(PARMtmaxquant[,c(5,9)],PARMtminquant[,c(5,9)])),lwd=2,type="l",col="red")
lines(probsin~PARMtminquant[,5],lwd=2,col="blue")

lines(probsin~PARMtmaxquant[,9],lwd=2,col="red",lty=2)
lines(probsin~PARMtminquant[,9],lwd=2,col="blue",lty=2)

plot(probsin~PARMdtrquant[,5],xlim=range(c(PARMdtrquant[,c(5,9)])),lwd=2,type="l",col="black")
lines(probsin~PARMdtrquant[,9],lwd=2,col="black",lty=2)

m=5
monidx = which(as.numeric(substr(histdates2,6,7))==m)
hist(PARMtmax[i,j,monidx],breaks=20)
hist(PARMtmin[i,j,monidx],breaks=20)
hist(PARMdtr2[i,j,monidx],breaks=20)

m=9
monidx = which(as.numeric(substr(histdates2,6,7))==m)
hist(PARMtmax[i,j,monidx],breaks=20)
hist(PARMtmin[i,j,monidx],breaks=20)
hist(PARMdtr2[i,j,monidx],breaks=20)



