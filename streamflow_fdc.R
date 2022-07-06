
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
library(hydroTSM)

gaugenum = "08151500"
futureperiod = c(2070,2099)

load(paste("/home/woot0002/streamflowreg_3daymovingmean_",gaugenum,".Rdata",sep=""))

# load DS files

prfile_hist = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_*historical*.nc",intern=TRUE))

tasmaxfile_hist = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_*historical*.nc",intern=TRUE))

tasminfile_hist = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_*historical*.nc",intern=TRUE))

prfile_proj = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_day*rcp*.nc",intern=TRUE))

tasmaxfile_proj = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_day*rcp*.nc",intern=TRUE))

tasminfile_proj = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_day*rcp*.nc",intern=TRUE))

#####
# Create file breakdown tables

filebreakdown = do.call(rbind,strsplit(prfile_proj,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(prfile_hist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

################
# load data

load(file=paste("/home/woot0002/CPREP_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))

#################
# Flow Duration Curve - function

fdcurve = function(x){
  
  if(length(which(x<=0))<length(x)){
  if(length(which(x<=0))>0){
    x[which(x<=0)]=NA
  }
  if(length(which(is.na(x)==TRUE))>0){
    x = x[which(is.na(x)==FALSE)]
  }
  
  ordobs = order(x,decreasing=TRUE)
  
  xobs = x[ordobs]
  ranks = 1:length(xobs)
  n = length(xobs)
  
  P = 100*(ranks/(n+1))
  
  outframe = data.frame(xobs,P)
  } else {
    xobs = rep(NA,length(x))
    ranks = 1:length(xobs)
    n = length(xobs)
    P = 100*(ranks/(n+1))
    outframe = data.frame(xobs,P)
  }
  outframe
  
}

#################
# obs regression fdc

obsfdc = fdcurve(fulldata2$streamflowmod)

#################
# hist DS fdc

histfdc = list()

for(i in 1:length(histresults)){
  histtmp = histresults[[i]][[1]]
  histfdc[[i]] = fdcurve(histtmp)
}

#################
# proj DS fdc

projfdc = list()

for(i in 1:length(projresults)){
  projtmp = projresults[[i]][[1]]
  projfdc[[i]] = fdcurve(projtmp)
}

#################
# Figure out scales for y axis

histonly_y = c(obsfdc[,1])
for(i in 1:length(histfdc)){
  histonly_y=c(histonly_y,histfdc[[i]][,1])
}

histonly_ylog = log(histonly_y,10)
histonly_ylog = histonly_ylog[which(histonly_ylog!= -Inf)]
histylim = range(histonly_ylog)

projonly_y = c()
for(i in 1:length(histfdc)){
  projonly_y=c(projonly_y,histfdc[[i]][,1])
}

for(i in 1:length(projfdc)){
  projonly_y=c(projonly_y,projfdc[[i]][,1])
}

projonly_ylog = log(projonly_y,10)
projonly_ylog = projonly_ylog[which(projonly_ylog!= -Inf)]
projylim = range(projonly_ylog)

#################
# plotting fdc curves

#plot(log(xobs,10)~P,data=obsfdc,axes=FALSE,type="l",col="black",ylim=histylim,xlab="Exceedance Probability (%)",ylab="Streamflow (cfs)")
#axis(side = 1)
#ticks <- c(0.01,0.1,1,10,100,1000,10000)
#locs <- c(-2,-1,0,1,2,3,4)
#axis(side = 2, at = locs, labels = ticks)

# historical only
#for(i in 1:28){
#  if(i==1){
#    plot(log(xobs,10)~P,data=obsfdc,type="l",col="black",ylim=histylim,xlab="Exceedance Probability (%)",ylab="Log(streamflow)")
#  } else {
#    lines(log(xobs,10)~P,data=histfdc[[i-1]],type="l",lwd=0.5,col="gray")
#  }
#}

# historical + projected
for(i in 1:27){
  if(i==1){
    plot(log(xobs,10)~P,data=histfdc[[i]],type="l",axes=FALSE,col="black",xlim=c(0,100),ylim=projylim,xlab="Exceedance Probability (%)",ylab="Streamflow (cfs)",main=paste("Flow Duration Curves: ",gaugenum,sep=""))
  } else {
    lines(log(xobs,10)~P,data=histfdc[[i]],type="l",lwd=0.5,col="gray")
  }
}

for(i in 1:81){
  
  if(projfilebreakdown$scen[i]=="rcp26"){
    colin ="green"
  }
  
  if(projfilebreakdown$scen[i]=="rcp45"){
    colin ="blue"
  }
  
  if(projfilebreakdown$scen[i]=="rcp85"){
    colin ="red"
  }
  
  lines(log(xobs,10)~P,data=projfdc[[i]],type="l",lwd=0.5,col=colin)
  
}

axis(side = 1)
ticks <- c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
locs <- c(-4,-3,-2,-1,0,1,2,3,4)
axis(side = 2, at = locs, labels = ticks)


legend("bottomleft",
       legend=c("Historical (1981-2005)",paste("RCP 2.6 (",futureperiod[1],"-",futureperiod[2],")",sep=""),
       paste("RCP 4.5 (",futureperiod[1],"-",futureperiod[2],")",sep=""),paste("RCP 8.5 (",futureperiod[1],"-",futureperiod[2],")",sep="")),
       col=c("gray","green","blue","red"),
       lwd=0.5)

######
# calculate ensemble means

rcp26idx = which(projfilebreakdown$scen=="rcp26")
rcp45idx = which(projfilebreakdown$scen=="rcp45")
rcp85idx = which(projfilebreakdown$scen=="rcp85")

probsused = seq(0.01,99.99,by=0.01)

projdates = seq(as.Date(paste(futureperiod[1],"-01-01",sep="")),as.Date(paste(futureperiod[2],"-12-31",sep="")),by="day")

for(i in 1:length(rcp26idx)){
  
  projtmp = projfdc[[rcp26idx[i]]]
  
  xobstmp = c()
  for(j in 1:length(probsused)){
    diffP = projtmp[,2]-probsused[j]
    useidx = which(abs(diffP)==min(abs(diffP)))
    if(length(useidx)==1){
      xobstmp[j] = projtmp[useidx,1]
    }
    if(length(useidx)>1){
      xobstmp[j] = projtmp[useidx[1],1]
    }
    if(length(useidx)==0){
      xobstmp[j] = NA
    }
  }
  
  if(i==1){
    rcp26table = data.frame(probsused,xobstmp)
    names(rcp26table) = c("P",paste(projfilebreakdown$GCM[rcp26idx[i]],projfilebreakdown$DS[rcp26idx[i]],projfilebreakdown$obs[rcp26idx[i]],sep="_"))
  } else {
    rcp26table = cbind(rcp26table,xobstmp)
    names(rcp26table)[ncol(rcp26table)] = paste(projfilebreakdown$GCM[rcp26idx[i]],projfilebreakdown$DS[rcp26idx[i]],projfilebreakdown$obs[rcp26idx[i]],sep="_")
    }
}

rcp26table$ensmean = apply(rcp26table[,2:28],1,mean,na.rm=TRUE)
rcp26table$ensmin = apply(rcp26table[,2:28],1,min,na.rm=TRUE)
rcp26table$ensmax = apply(rcp26table[,2:28],1,max,na.rm=TRUE)

###

for(i in 1:length(rcp45idx)){
  
  projtmp = projfdc[[rcp45idx[i]]]
  
  xobstmp = c()
  for(j in 1:length(probsused)){
    diffP = projtmp[,2]-probsused[j]
    useidx = which(abs(diffP)==min(abs(diffP)))
    if(length(useidx)==1){
      xobstmp[j] = projtmp[useidx,1]
    }
    if(length(useidx)>1){
      xobstmp[j] = projtmp[useidx[1],1]
    }
    if(length(useidx)==0){
      xobstmp[j] = NA
    }
  }
  
  if(i==1){
    rcp45table = data.frame(probsused,xobstmp)
    names(rcp45table) = c("P",paste(projfilebreakdown$GCM[rcp45idx[i]],projfilebreakdown$DS[rcp45idx[i]],projfilebreakdown$obs[rcp45idx[i]],sep="_"))
  } else {
    rcp45table = cbind(rcp45table,xobstmp)
    names(rcp45table)[ncol(rcp45table)] = paste(projfilebreakdown$GCM[rcp45idx[i]],projfilebreakdown$DS[rcp45idx[i]],projfilebreakdown$obs[rcp45idx[i]],sep="_")
  }
}

rcp45table$ensmean = apply(rcp45table[,2:28],1,mean,na.rm=TRUE)
rcp45table$ensmin = apply(rcp45table[,2:28],1,min,na.rm=TRUE)
rcp45table$ensmax = apply(rcp45table[,2:28],1,max,na.rm=TRUE)

###

for(i in 1:length(rcp85idx)){
  
  projtmp = projfdc[[rcp85idx[i]]]
  
  xobstmp = c()
  for(j in 1:length(probsused)){
    diffP = projtmp[,2]-probsused[j]
    useidx = which(abs(diffP)==min(abs(diffP)))
    if(length(useidx)==1){
      xobstmp[j] = projtmp[useidx,1]
    }
    if(length(useidx)>1){
      xobstmp[j] = projtmp[useidx[1],1]
    }
    if(length(useidx)==0){
      xobstmp[j] = NA
    }
  }
  
  if(i==1){
    rcp85table = data.frame(probsused,xobstmp)
    names(rcp85table) = c("P",paste(projfilebreakdown$GCM[rcp85idx[i]],projfilebreakdown$DS[rcp85idx[i]],projfilebreakdown$obs[rcp85idx[i]],sep="_"))
  } else {
    rcp85table = cbind(rcp85table,xobstmp)
    names(rcp85table)[ncol(rcp85table)] = paste(projfilebreakdown$GCM[rcp85idx[i]],projfilebreakdown$DS[rcp85idx[i]],projfilebreakdown$obs[rcp85idx[i]],sep="_")
  }
}

rcp85table$ensmean = apply(rcp85table[,2:28],1,mean,na.rm=TRUE)
rcp85table$ensmin = apply(rcp85table[,2:28],1,min,na.rm=TRUE)
rcp85table$ensmax = apply(rcp85table[,2:28],1,max,na.rm=TRUE)

###

for(i in 1:nrow(histfilebreakdown)){
  
  projtmp = histfdc[[i]]
  
  xobstmp = c()
  for(j in 1:length(probsused)){
    diffP = projtmp[,2]-probsused[j]
    useidx = which(abs(diffP)==min(abs(diffP)))
    if(length(useidx)==1){
      xobstmp[j] = projtmp[useidx,1]
    }
    if(length(useidx)>1){
      xobstmp[j] = projtmp[useidx[1],1]
    }
    if(length(useidx)==0){
      xobstmp[j] = NA
    }
  }
  
  if(i==1){
    histtable = data.frame(probsused,xobstmp)
    names(histtable) = c("P",paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_"))
  } else {
    histtable = cbind(histtable,xobstmp)
    names(histtable)[ncol(histtable)] = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")
  }
}

histtable$ensmean = apply(histtable[,2:28],1,mean,na.rm=TRUE)
histtable$ensmin = apply(histtable[,2:28],1,min,na.rm=TRUE)
histtable$ensmax = apply(histtable[,2:28],1,max,na.rm=TRUE)

###########
# plot ens means

# historical only
#for(i in 1:28){
#  if(i==1){
#    plot(log(histtable$ensmean,10)~histtable$P,type="l",lwd=2,col="black",ylim=histylim,xlab="Exceedance Probability (%)",ylab="Log(streamflow)")
#  } else {
#    lines(log(histtable[,i],10)~histtable$P,type="l",lwd=0.5,col="gray")
#  }
#}

plot(log(histtable$ensmean,10)~histtable$P,axes=FALSE,type="l",lwd=1,lty=2,col="black",xlab="Exceedance Probability (%)",ylab="Streamflow (cfs)",main=paste("Flow Duration Curves: ",gaugenum,sep=""))
lines(log(rcp26table$ensmean,10)~rcp26table$P,type="l",lwd=1,col="green")
lines(log(rcp45table$ensmean,10)~rcp45table$P,type="l",lwd=1,col="blue")
lines(log(rcp85table$ensmean,10)~rcp85table$P,type="l",lwd=1,col="red")


axis(side = 1)
ticks <- c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
locs <- c(-4,-3,-2,-1,0,1,2,3,4)
axis(side = 2, at = locs, labels = ticks)


legend("bottomleft",
       legend=c("Historical (1981-2005)",paste("RCP 2.6 (",futureperiod[1],"-",futureperiod[2],")",sep=""),
                paste("RCP 4.5 (",futureperiod[1],"-",futureperiod[2],")",sep=""),paste("RCP 8.5 (",futureperiod[1],"-",futureperiod[2],")",sep="")),
       col=c("black","green","blue","red"),
       lwd=2,lty=c(2,1,1,1))

###

plot(log(histtable$ensmean,10)~histtable$P,axes=FALSE,type="l",lwd=1,lty=1,col="black",xlab="Exceedance Probability (%)",ylab="Streamflow (cfs)",main=paste("Flow Duration Curves: ",gaugenum,sep=""))
lines(log(rcp26table$ensmean,10)~rcp26table$P,type="l",lwd=1,col="green")
lines(log(rcp45table$ensmean,10)~rcp45table$P,type="l",lwd=1,col="blue")
lines(log(rcp85table$ensmean,10)~rcp85table$P,type="l",lwd=1,col="red")

lines(log(rcp26table$ensmin,10)~rcp26table$P,type="l",lwd=0.5,col="green",lty=2)
lines(log(rcp45table$ensmin,10)~rcp45table$P,type="l",lwd=0.5,col="blue",lty=2)
lines(log(rcp85table$ensmin,10)~rcp85table$P,type="l",lwd=0.5,col="red",lty=2)

lines(log(rcp26table$ensmax,10)~rcp26table$P,type="l",lwd=0.5,col="green",lty=2)
lines(log(rcp45table$ensmax,10)~rcp45table$P,type="l",lwd=0.5,col="blue",lty=2)
lines(log(rcp85table$ensmax,10)~rcp85table$P,type="l",lwd=0.5,col="red",lty=2)

lines(log(histtable$ensmin,10)~histtable$P,type="l",lwd=0.5,col="black",lty=2)
lines(log(histtable$ensmax,10)~histtable$P,type="l",lwd=0.5,col="black",lty=2)


axis(side = 1)
ticks <- c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
locs <- c(-4,-3,-2,-1,0,1,2,3,4)
axis(side = 2, at = locs, labels = ticks)

legend("bottomleft",
       legend=c("Historical (1981-2005)",paste("RCP 2.6 (",futureperiod[1],"-",futureperiod[2],")",sep=""),
                paste("RCP 4.5 (",futureperiod[1],"-",futureperiod[2],")",sep=""),paste("RCP 8.5 (",futureperiod[1],"-",futureperiod[2],")",sep="")),
       col=c("black","green","blue","red"),
       lwd=2,lty=1)

###########
# filled polygons version

plot(log(histtable$ensmean,10)~histtable$P,axes=FALSE,type="l",lwd=2,lty=2,col="black",xlab="Exceedance Probability (%)",ylab="Streamflow (cfs)",main=paste("Flow Duration Curves: ",gaugenum,sep=""))
polygon(c(rev(histtable$P), histtable$P), c(rev(log(histtable$ensmin,10)), log(histtable$ensmax,10)),col = rgb(224/255,224/255,224/255,0.5), border = NA)

lines(log(rcp26table$ensmean,10)~rcp26table$P,type="l",lwd=2,lty=2,col="green")
polygon(c(rev(rcp26table$P), rcp26table$P), c(rev(log(rcp26table$ensmin,10)), log(rcp26table$ensmax,10)),col = rgb(0/255,255/255,0/255,0.25), border = NA)

lines(log(rcp45table$ensmean,10)~rcp45table$P,type="l",lwd=2,lty=2,col="blue")
polygon(c(rev(rcp45table$P), rcp45table$P), c(rev(log(rcp45table$ensmin,10)), log(rcp45table$ensmax,10)),col = rgb(0/255,0/255,255/255,0.25), border = NA)

lines(log(rcp85table$ensmean,10)~rcp85table$P,type="l",lwd=2,lty=2,col="red")
polygon(c(rev(rcp85table$P), rcp85table$P), c(rev(log(rcp85table$ensmin,10)), log(rcp85table$ensmax,10)),col = rgb(255/255,0/255,0/255,0.25), border = NA)

axis(side = 1)
ticks <- c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
locs <- c(-4,-3,-2,-1,0,1,2,3,4)
axis(side = 2, at = locs, labels = ticks)

legend("bottomleft",
       legend=c("Historical (1981-2005)",paste("RCP 2.6 (",futureperiod[1],"-",futureperiod[2],")",sep=""),
                paste("RCP 4.5 (",futureperiod[1],"-",futureperiod[2],")",sep=""),paste("RCP 8.5 (",futureperiod[1],"-",futureperiod[2],")",sep="")),
       col=c("black","green","blue","red"),
       lwd=2,lty=1)


####
# limiting range


plot(log(histtable$ensmean,10)~histtable$P,ylim=c(-2.5,1.7),axes=FALSE,type="l",lwd=1,lty=1,col="black",xlab="Exceedance Probability (%)",ylab="Streamflow (cfs)",main=paste("Flow Duration Curves: ",gaugenum,sep=""))
lines(log(rcp26table$ensmean,10)~rcp26table$P,type="l",lwd=1,col="green")
lines(log(rcp45table$ensmean,10)~rcp45table$P,type="l",lwd=1,col="blue")
lines(log(rcp85table$ensmean,10)~rcp85table$P,type="l",lwd=1,col="red")

lines(log(rcp26table$ensmin,10)~rcp26table$P,type="l",lwd=0.5,col="green",lty=2)
lines(log(rcp45table$ensmin,10)~rcp45table$P,type="l",lwd=0.5,col="blue",lty=2)
lines(log(rcp85table$ensmin,10)~rcp85table$P,type="l",lwd=0.5,col="red",lty=2)

lines(log(rcp26table$ensmax,10)~rcp26table$P,type="l",lwd=0.5,col="green",lty=2)
lines(log(rcp45table$ensmax,10)~rcp45table$P,type="l",lwd=0.5,col="blue",lty=2)
lines(log(rcp85table$ensmax,10)~rcp85table$P,type="l",lwd=0.5,col="red",lty=2)

lines(log(histtable$ensmin,10)~histtable$P,type="l",lwd=0.5,col="black",lty=2)
lines(log(histtable$ensmax,10)~histtable$P,type="l",lwd=0.5,col="black",lty=2)


axis(side = 1)
ticks <- c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
locs <- c(-4,-3,-2,-1,0,1,2,3,4)
axis(side = 2, at = locs, labels = ticks)

legend("bottomleft",
       legend=c("Historical (1981-2005)",paste("RCP 2.6 (",futureperiod[1],"-",futureperiod[2],")",sep=""),
                paste("RCP 4.5 (",futureperiod[1],"-",futureperiod[2],")",sep=""),paste("RCP 8.5 (",futureperiod[1],"-",futureperiod[2],")",sep="")),
       col=c("black","green","blue","red"),
       lwd=2,lty=1)




