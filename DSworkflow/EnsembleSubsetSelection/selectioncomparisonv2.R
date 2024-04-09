source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(mapdata)
library(maptools)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(ggplot2)
library(modi)
library(parallel)
library(foreach)
library(doParallel)


#####
# Version of this script to run parallel but also one domain as a time.

enssize = 5
mask = "SGP-NCA"
ensgroup="CMIP6"

#####
# load RMSEs

load(paste("/home/woot0002/RCMES/allvars_",ensgroup,"_optselect_enssize",enssize,"_dom",mask,".Rdata",sep=""))

pr_optim = alldat2pr
tasmax_optim = alldat2tx
tasmin_optim = alldat2tn
mv_optim = alldat2mv

histfile1 = paste("/home/woot0002/RCMES/pr_",ensgroup,"_climo_",mask,"_annual.nc",sep="")
histfile2 = paste("/home/woot0002/RCMES/tasmax_",ensgroup,"_climo_",mask,"_annual.nc",sep="")
histfile3 = paste("/home/woot0002/RCMES/tasmin_",ensgroup,"_climo_",mask,"_annual.nc",sep="")

#####
# pull historical data - pr

nctest = nc_open(histfile1)
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")

varnamesprh = c()
prhist = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamesprh[i]=nctest$var[[i]]$name
  prhist[,,i]=ncvar_get(nctest,varnamesprh[i])
}

nc_close(nctest)

#####
# pull historical data - tasmax

nctest = nc_open(histfile2)
lathist = ncvar_get(nctest,"lat")
lonhist=ncvar_get(nctest,"lon")
varnamestxh = c()
tasmaxhist = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamestxh[i]=nctest$var[[i]]$name
  tasmaxhist[,,i]=ncvar_get(nctest,varnamestxh[i])
}
nc_close(nctest)

#####
# pull historical data - tasmin

nctest = nc_open(histfile3)
lathist = ncvar_get(nctest,"lat")
lonhist=ncvar_get(nctest,"lon")
varnamestnh = c()
tasminhist = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamestnh[i]=nctest$var[[i]]$name
  tasminhist[,,i]=ncvar_get(nctest,varnamestnh[i])
}
nc_close(nctest)



######
# pull projected values

scenario = "rcp85"

projfile1 = paste("/home/woot0002/RCMES/pr_",ensgroup,"_climo_",scenario,"_",mask,"_annual.nc",sep="")
projfile2 = paste("/home/woot0002/RCMES/tasmax_",ensgroup,"_climo_",scenario,"_",mask,"_annual.nc",sep="")
projfile3 = paste("/home/woot0002/RCMES/tasmin_",ensgroup,"_climo_",scenario,"_",mask,"_annual.nc",sep="")

#####
# pull projected data - pr

nctest = nc_open(projfile1)
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")

varnamesprp = c()
prproj = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamesprp[i]=nctest$var[[i]]$name
  prproj[,,i]=ncvar_get(nctest,varnamesprp[i])
}

nc_close(nctest)

#####
# pull projected data - tasmax

nctest = nc_open(projfile2)
varnamestxp = c()
tasmaxproj = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamestxp[i]=nctest$var[[i]]$name
  tasmaxproj[,,i]=ncvar_get(nctest,varnamestxp[i])
}
nc_close(nctest)

#####
# pull projected data - tasmin

nctest = nc_open(projfile3)
varnamestnp = c()
tasminproj = array(NA,dim=c(length(lon),length(lat),length(nctest$var)))
for(i in 1:length(nctest$var)){
  varnamestnp[i]=nctest$var[[i]]$name
  tasminproj[,,i]=ncvar_get(nctest,varnamestnp[i])
}
nc_close(nctest)

######
# Calculate projected change and mask

probs = prhist[,,which(varnamesprh=="OBS")]
tasmaxobs = tasmaxhist[,,which(varnamestxh=="OBS")]
tasminobs = tasminhist[,,which(varnamestnh=="OBS")]


varidx = which(varnamesprh!="OBS" & varnamesprh!="ENS")
varnamesprh = varnamesprh[varidx]
prhist = prhist[,,varidx]

varidx = which(varnamesprp!="OBS" & varnamesprp!="ENS")
varnamesprp = varnamesprp[varidx]
prproj = prproj[,,varidx]

varidx = which(varnamestxh!="OBS" & varnamestxh!="ENS")
varnamestxh = varnamestxh[varidx]
tasmaxhist = tasmaxhist[,,varidx]

varidx = which(varnamestxp!="OBS" & varnamestxp!="ENS")
varnamestxp = varnamestxp[varidx]
tasmaxproj = tasmaxproj[,,varidx]

varidx = which(varnamestnh!="OBS" & varnamestnh!="ENS")
varnamestnh = varnamestnh[varidx]
tasminhist = tasminhist[,,varidx]

varidx = which(varnamestnp!="OBS" & varnamestnp!="ENS")
varnamestnp = varnamestnp[varidx]
tasminproj = tasminproj[,,varidx]

tasminchange = tasminproj-tasminhist
tasmaxchange = tasmaxproj-tasmaxhist
prchange = prproj-prhist

#####
# calc mean changes

pro_GCMs=txo_GCMs=tno_GCMs=mvo_GCMs=c()
for(g in 8:12){
  pro_GCMs[g-7]=as.character(pr_optim[1,g])
  txo_GCMs[g-7]=as.character(tasmax_optim[1,g])
  tno_GCMs[g-7]=as.character(tasmin_optim[1,g])
  mvo_GCMs[g-7]=as.character(mv_optim[1,(g-3)])
}

varidx_pro = which(varnamesprp %in% pro_GCMs)
varidx_txo = which(varnamestxp %in% txo_GCMs)
varidx_tno = which(varnamestnp %in% tno_GCMs)
varidx_mvo = which(varnamesprp %in% mvo_GCMs)

prmeanall = apply(prchange,c(1,2),mean,na.rm=TRUE)
txmeanall = apply(tasmaxchange,c(1,2),mean,na.rm=TRUE)
tnmeanall = apply(tasminchange,c(1,2),mean,na.rm=TRUE)

prmean_pro = apply(prchange[,,varidx_pro],c(1,2),mean,na.rm=TRUE)
prmean_txo = apply(prchange[,,varidx_txo],c(1,2),mean,na.rm=TRUE)
prmean_tno = apply(prchange[,,varidx_tno],c(1,2),mean,na.rm=TRUE)
prmean_mvo = apply(prchange[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)

txmean_pro = apply(tasmaxchange[,,varidx_pro],c(1,2),mean,na.rm=TRUE)
txmean_txo = apply(tasmaxchange[,,varidx_txo],c(1,2),mean,na.rm=TRUE)
txmean_tno = apply(tasmaxchange[,,varidx_tno],c(1,2),mean,na.rm=TRUE)
txmean_mvo = apply(tasmaxchange[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)

tnmean_pro = apply(tasminchange[,,varidx_pro],c(1,2),mean,na.rm=TRUE)
tnmean_txo = apply(tasminchange[,,varidx_txo],c(1,2),mean,na.rm=TRUE)
tnmean_tno = apply(tasminchange[,,varidx_tno],c(1,2),mean,na.rm=TRUE)
tnmean_mvo = apply(tasminchange[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)

#####
# error calcs

prmeanall_error = apply(prhist,c(1,2),mean,na.rm=TRUE)-probs
prmean_pro_error = apply(prhist[,,varidx_pro],c(1,2),mean,na.rm=TRUE)-probs
prmean_txo_error = apply(prhist[,,varidx_txo],c(1,2),mean,na.rm=TRUE)-probs
prmean_tno_error = apply(prhist[,,varidx_tno],c(1,2),mean,na.rm=TRUE)-probs
prmean_mvo_error = apply(prhist[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)-probs

txmeanall_error = apply(tasmaxhist,c(1,2),mean,na.rm=TRUE)-tasmaxobs
txmean_pro_error = apply(tasmaxhist[,,varidx_pro],c(1,2),mean,na.rm=TRUE)-tasmaxobs
txmean_txo_error = apply(tasmaxhist[,,varidx_txo],c(1,2),mean,na.rm=TRUE)-tasmaxobs
txmean_tno_error = apply(tasmaxhist[,,varidx_tno],c(1,2),mean,na.rm=TRUE)-tasmaxobs
txmean_mvo_error = apply(tasmaxhist[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)-tasmaxobs

tnmeanall_error = apply(tasminhist,c(1,2),mean,na.rm=TRUE)-tasminobs
tnmean_pro_error = apply(tasminhist[,,varidx_pro],c(1,2),mean,na.rm=TRUE)-tasminobs
tnmean_txo_error = apply(tasminhist[,,varidx_txo],c(1,2),mean,na.rm=TRUE)-tasminobs
tnmean_tno_error = apply(tasminhist[,,varidx_tno],c(1,2),mean,na.rm=TRUE)-tasminobs
tnmean_mvo_error = apply(tasminhist[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)-tasminobs


#####
# plot taylor diagram

save(list=c("probs","tasmaxobs","tasminobs","prhist","tasmaxhist","tasminhist","varidx_pro","varidx_txo","varidx_tno","varidx_mvo"),file=paste("/home/woot0002/RCMES/evaldata_",ensgroup,"_dom",mask,"_size",enssize,"_shortEAA.Rdata",sep=""))

library(plotrix)

plotfilename = paste("/home/woot0002/RCMES/all_taylordiagrams_",ensgroup,"_dom",mask,"_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=5,height=5)
taylor.diagram(ref=as.vector(probs),model=as.vector(apply(prhist,c(1,2),mean,na.rm=TRUE)),main=paste("Precip Taylor Diagram - ",mask,sep=""))
taylor.diagram(ref=as.vector(probs),model=as.vector(apply(prhist[,,varidx_pro],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="blue",pch=5)
taylor.diagram(ref=as.vector(probs),model=as.vector(apply(prhist[,,varidx_txo],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="darkgreen",pch=2)
taylor.diagram(ref=as.vector(probs),model=as.vector(apply(prhist[,,varidx_tno],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="yellow4",pch=2)
taylor.diagram(ref=as.vector(probs),model=as.vector(apply(prhist[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="purple",pch=1)
legend("topright",legend=c("Full","pr optim", "tasmax optim", "tasmin optim","mv optim"),pch=c(19,5,2,2,1),col=c("red","blue","darkgreen","yellow4","purple"),cex=0.75)

taylor.diagram(ref=as.vector(tasmaxobs),model=as.vector(apply(tasmaxhist,c(1,2),mean,na.rm=TRUE)),main=paste("tasmax Taylor Diagram - ",mask,sep=""))
taylor.diagram(ref=as.vector(tasmaxobs),model=as.vector(apply(tasmaxhist[,,varidx_pro],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="blue",pch=5)
taylor.diagram(ref=as.vector(tasmaxobs),model=as.vector(apply(tasmaxhist[,,varidx_txo],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="darkgreen",pch=2)
taylor.diagram(ref=as.vector(tasmaxobs),model=as.vector(apply(tasmaxhist[,,varidx_tno],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="yellow4",pch=2)
taylor.diagram(ref=as.vector(tasmaxobs),model=as.vector(apply(tasmaxhist[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="purple",pch=1)
legend("topright",legend=c("Full","pr optim", "tasmax optim", "tasmin optim","mv optim"),pch=c(19,5,2,2,1),col=c("red","blue","darkgreen","yellow4","purple"),cex=0.75)

taylor.diagram(ref=as.vector(tasminobs),model=as.vector(apply(tasminhist,c(1,2),mean,na.rm=TRUE)),main=paste("tasmin Taylor Diagram - ",mask,sep=""))
taylor.diagram(ref=as.vector(tasminobs),model=as.vector(apply(tasminhist[,,varidx_pro],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="blue",pch=5)
taylor.diagram(ref=as.vector(tasminobs),model=as.vector(apply(tasminhist[,,varidx_txo],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="darkgreen",pch=2)
taylor.diagram(ref=as.vector(tasminobs),model=as.vector(apply(tasminhist[,,varidx_tno],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="yellow4",pch=2)
taylor.diagram(ref=as.vector(tasminobs),model=as.vector(apply(tasminhist[,,varidx_mvo],c(1,2),mean,na.rm=TRUE)),add=TRUE,col="purple",pch=1)
legend("topright",legend=c("Full","pr optim", "tasmax optim", "tasmin optim","mv optim"),pch=c(19,5,2,2,1),col=c("red","blue","darkgreen","yellow4","purple"),cex=0.75)

dev.off()

#####
# plot 3 member blocks

if(mask!="full"){
  if(mask=="NE-NCA"){
    xlimin = c(-86,-64)
    ylimin = c(36,48)
    H=7
    W=10
  }
  if(mask=="SW-NCA"){
    xlimin = c(-125,-102)
    ylimin = c(31,42)
    H=7
    W=10
  }
  if(mask=="SGP-NCA"){
    xlimin = c(-107,-93)
    ylimin = c(25,40)
    H=10
    W=5
    
    Ha=20
    Wa=10
    textx = -106.5
    texty1 =27.45
    texty2 =26.45
    texty3 =25.45
  }
  if(mask=="SE-NCA"){
    xlimin = c(-95,-75)
    ylimin = c(24,40)
    H=6
    W=8
  }
  if(mask=="MW-NCA"){
    xlimin = c(-99,-80)
    ylimin = c(36,50)
    H=7
    W=6
  }
  if(mask=="NGP-NCA"){
    xlimin = c(-118,-95)
    ylimin = c(40,49)
    H=7
    W=6
  }
  if(mask=="NW-NCA"){
    xlimin = c(-125,-111)
    ylimin = c(42,49)
    H=7
    W=6
  }
} else {
  xlimin = range(lonhist)
  ylimin = range(lathist)
  H=7
  W=10
  
  Ha=14
  Wa=20
  textx = -125
  texty1 =28.45
  texty2 =26.45
  texty3 =24.45
}

#####

range(c(prmeanall,prmean_pro,prmean_txo,prmean_tno,prmean_mvo,prchange[,,varidx_pro],prchange[,,varidx_txo],prchange[,,varidx_tno],prchange[,,varidx_mvo]),na.rm=TRUE)

colorchoicediff="browntogreen"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/pr_changeplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_block_shortEAA.pdf",sep="") 
diffcolorbar = colorramp(c(prmeanall,prmean_pro,prmean_txo,prmean_tno,prmean_mvo,prchange[,,varidx_pro],prchange[,,varidx_txo],prchange[,,varidx_tno],prchange[,,varidx_mvo]),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(-250,250))

pdf(plotfilename,onefile=TRUE,width=15,height=26)

par(mfrow=c(4,6))

for(i in 1:24){
  if(i<=5){
    testsfc1 = list(x=lonhist,y=lathist,z=prchange[,,varidx_pro[i]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=pro_GCMs[i],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==6){
    testsfc1 = list(x=lonhist,y=lathist,z=prmean_pro)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - pr optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  if(i>=7 & i<=11){
    testsfc1 = list(x=lonhist,y=lathist,z=prchange[,,varidx_txo[(i-6)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=txo_GCMs[(i-6)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==12){
    testsfc1 = list(x=lonhist,y=lathist,z=prmean_txo)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - tasmax optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  if(i>=13 & i<=17){
    testsfc1 = list(x=lonhist,y=lathist,z=prchange[,,varidx_tno[(i-12)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=tno_GCMs[(i-12)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==18){
    testsfc1 = list(x=lonhist,y=lathist,z=prmean_tno)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - tasmin optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=19 & i<=23){
    testsfc1 = list(x=lonhist,y=lathist,z=prchange[,,varidx_mvo[(i-18)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=mvo_GCMs[(i-18)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==24){
    testsfc1 = list(x=lonhist,y=lathist,z=prmean_mvo)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - mv optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(prmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(prmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
}

dev.off()

#####


range(c(txmeanall,txmean_pro,txmean_txo,txmean_tno,txmean_mvo,tasmaxchange[,,varidx_pro],tasmaxchange[,,varidx_txo],tasmaxchange[,,varidx_tno],tasmaxchange[,,varidx_mvo]),na.rm=TRUE)

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmax_changeplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_block_shortEAA.pdf",sep="") 
diffcolorbar = colorramp(c(txmeanall,txmean_pro,txmean_txo,txmean_tno,txmean_mvo,tasmaxchange[,,varidx_pro],tasmaxchange[,,varidx_txo],tasmaxchange[,,varidx_tno],tasmaxchange[,,varidx_mvo]),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,10))

pdf(plotfilename,onefile=TRUE,width=15,height=26)

par(mfrow=c(4,6))

for(i in 1:24){
  if(i<=5){
    testsfc1 = list(x=lonhist,y=lathist,z=tasmaxchange[,,varidx_pro[i]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=pro_GCMs[i],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasmaxchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasmaxchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasmaxchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==6){
    testsfc1 = list(x=lonhist,y=lathist,z=txmean_pro)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - pr optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(txmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(txmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  if(i>=7 & i<=11){
    testsfc1 = list(x=lonhist,y=lathist,z=tasmaxchange[,,varidx_txo[(i-6)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=txo_GCMs[(i-6)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasmaxchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasmaxchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasmaxchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==12){
    testsfc1 = list(x=lonhist,y=lathist,z=txmean_txo)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - tasmax optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(txmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(txmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  if(i>=13 & i<=17){
    testsfc1 = list(x=lonhist,y=lathist,z=tasmaxchange[,,varidx_tno[(i-12)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=tno_GCMs[(i-12)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasmaxchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasmaxchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasmaxchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==18){
    testsfc1 = list(x=lonhist,y=lathist,z=txmean_tno)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - tasmin optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(txmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(txmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=19 & i<=23){
    testsfc1 = list(x=lonhist,y=lathist,z=tasmaxchange[,,varidx_mvo[(i-18)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=mvo_GCMs[(i-18)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasmaxchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasmaxchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasmaxchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==24){
    testsfc1 = list(x=lonhist,y=lathist,z=txmean_mvo)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - mv optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(txmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(txmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
}

dev.off()

#####


range(c(tnmeanall,tnmean_pro,tnmean_txo,tnmean_tno,tnmean_mvo,tasminchange[,,varidx_pro],tasminchange[,,varidx_txo],tasminchange[,,varidx_tno],tasminchange[,,varidx_mvo]),na.rm=TRUE)

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmin_changeplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_block_shortEAA.pdf",sep="") 
diffcolorbar = colorramp(c(tnmeanall,tnmean_pro,tnmean_txo,tnmean_tno,tnmean_mvo,tasminchange[,,varidx_pro],tasminchange[,,varidx_txo],tasminchange[,,varidx_tno],tasminchange[,,varidx_mvo]),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,10))

pdf(plotfilename,onefile=TRUE,width=15,height=26)

par(mfrow=c(4,6))

for(i in 1:24){
  if(i<=5){
    testsfc1 = list(x=lonhist,y=lathist,z=tasminchange[,,varidx_pro[i]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=pro_GCMs[i],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasminchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasminchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasminchange[,,varidx_pro[i]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==6){
    testsfc1 = list(x=lonhist,y=lathist,z=tnmean_pro)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - pr optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  if(i>=7 & i<=11){
    testsfc1 = list(x=lonhist,y=lathist,z=tasminchange[,,varidx_txo[(i-6)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=txo_GCMs[(i-6)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasminchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasminchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasminchange[,,varidx_txo[(i-6)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==12){
    testsfc1 = list(x=lonhist,y=lathist,z=tnmean_txo)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - tasmax optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  if(i>=13 & i<=17){
    testsfc1 = list(x=lonhist,y=lathist,z=tasminchange[,,varidx_tno[(i-12)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=tno_GCMs[(i-12)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasminchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasminchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasminchange[,,varidx_tno[(i-12)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==18){
    testsfc1 = list(x=lonhist,y=lathist,z=tnmean_tno)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - tasmin optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=19 & i<=23){
    testsfc1 = list(x=lonhist,y=lathist,z=tasminchange[,,varidx_mvo[(i-18)]])
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main=mvo_GCMs[(i-18)],zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tasminchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tasminchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tasminchange[,,varidx_mvo[(i-18)]],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    
  }
  if(i==24){
    testsfc1 = list(x=lonhist,y=lathist,z=tnmean_mvo)
    surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="mean - mv optim",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
    map("worldHires",add=TRUE)
    text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
}

dev.off()


#####
# change spread

prmeanall_spread = apply(prchange,c(1,2),max,na.rm=FALSE)-apply(prchange,c(1,2),min,na.rm=FALSE)
prmean_pro_spread = apply(prchange[,,varidx_pro],c(1,2),max,na.rm=FALSE)-apply(prchange[,,varidx_pro],c(1,2),min,na.rm=FALSE)
prmean_txo_spread = apply(prchange[,,varidx_txo],c(1,2),max,na.rm=FALSE)-apply(prchange[,,varidx_txo],c(1,2),min,na.rm=FALSE)
prmean_tno_spread = apply(prchange[,,varidx_tno],c(1,2),max,na.rm=FALSE)-apply(prchange[,,varidx_tno],c(1,2),min,na.rm=FALSE)
prmean_mvo_spread = apply(prchange[,,varidx_mvo],c(1,2),max,na.rm=FALSE)-apply(prchange[,,varidx_mvo],c(1,2),min,na.rm=FALSE)

txmeanall_spread = apply(tasmaxchange,c(1,2),max,na.rm=FALSE)-apply(tasmaxchange,c(1,2),min,na.rm=FALSE)
txmean_pro_spread = apply(tasmaxchange[,,varidx_pro],c(1,2),max,na.rm=FALSE)-apply(tasmaxchange[,,varidx_pro],c(1,2),min,na.rm=FALSE)
txmean_txo_spread = apply(tasmaxchange[,,varidx_txo],c(1,2),max,na.rm=FALSE)-apply(tasmaxchange[,,varidx_txo],c(1,2),min,na.rm=FALSE)
txmean_tno_spread = apply(tasmaxchange[,,varidx_tno],c(1,2),max,na.rm=FALSE)-apply(tasmaxchange[,,varidx_tno],c(1,2),min,na.rm=FALSE)
txmean_mvo_spread = apply(tasmaxchange[,,varidx_mvo],c(1,2),max,na.rm=FALSE)-apply(tasmaxchange[,,varidx_mvo],c(1,2),min,na.rm=FALSE)

tnmeanall_spread = apply(tasminchange,c(1,2),max,na.rm=FALSE)-apply(tasminchange,c(1,2),min,na.rm=FALSE)
tnmean_pro_spread = apply(tasminchange[,,varidx_pro],c(1,2),max,na.rm=FALSE)-apply(tasminchange[,,varidx_pro],c(1,2),min,na.rm=FALSE)
tnmean_txo_spread = apply(tasminchange[,,varidx_txo],c(1,2),max,na.rm=FALSE)-apply(tasminchange[,,varidx_txo],c(1,2),min,na.rm=FALSE)
tnmean_tno_spread = apply(tasminchange[,,varidx_tno],c(1,2),max,na.rm=FALSE)-apply(tasminchange[,,varidx_tno],c(1,2),min,na.rm=FALSE)
tnmean_mvo_spread = apply(tasminchange[,,varidx_mvo],c(1,2),max,na.rm=FALSE)-apply(tasminchange[,,varidx_mvo],c(1,2),min,na.rm=FALSE)


prmean_pro_sr = prmean_pro_spread/prmeanall_spread
prmean_txo_sr = prmean_txo_spread/prmeanall_spread
prmean_tno_sr = prmean_tno_spread/prmeanall_spread
prmean_mvo_sr = prmean_mvo_spread/prmeanall_spread

txmean_pro_sr = txmean_pro_spread/txmeanall_spread
txmean_txo_sr = txmean_txo_spread/txmeanall_spread
txmean_tno_sr = txmean_tno_spread/txmeanall_spread
txmean_mvo_sr = txmean_mvo_spread/txmeanall_spread

tnmean_pro_sr = tnmean_pro_spread/tnmeanall_spread
tnmean_txo_sr = tnmean_txo_spread/tnmeanall_spread
tnmean_tno_sr = tnmean_tno_spread/tnmeanall_spread
tnmean_mvo_sr = tnmean_mvo_spread/tnmeanall_spread

#testsfc = list(x=lonhist,y=lathist,z=prmean_pro_sr)
#surface(testsfc,type="I",xlim=c(textx,-102),ylim=c(31,42))

#####
# plot results - pr

range(c(prmeanall,prmean_pro,prmean_txo,prmean_tno,prmean_mvo),na.rm=TRUE)

colorchoicediff="browntogreen"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/pr_changeplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(prmeanall,prmean_pro,prmean_txo,prmean_tno,prmean_mvo),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(-100,100))

testsfc1 = list(x=lonhist,y=lathist,z=prmeanall)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - full ensemble mean change",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_pro)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean change, pr optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_txo)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean change, tasmax optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_tno)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean change, tasmin optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=prmean_mvo)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean change, mv optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()

####
# plot results - pr error

range(c(prmeanall_error,prmean_pro_error,prmean_txo_error,prmean_tno_error,prmean_mvo_error),na.rm=TRUE)

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/pr_errorplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(prmeanall_error,prmean_pro_error,prmean_txo_error,prmean_tno_error,prmean_mvo_error),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(-700,700))

testsfc1 = list(x=lonhist,y=lathist,z=prmeanall_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - full ensemble mean error",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_pro_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean error, pr optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_txo_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean error, tasmax optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_tno_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean error, tasmin optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_mvo_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - ensemble mean error, mv optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()

####
# plot results - pr spread

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="ratio"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/pr_spreadplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,".pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(prmean_pro_sr,prmean_r5o_sr,prmean_mvo_sr),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,1))

testsfc1 = list(x=lonhist,y=lathist,z=prmean_pro_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - spread ratio (pr optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_txo_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - spread ratio (tasmax optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=prmean_tno_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - spread ratio (tasmin optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=prmean_mvo_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="pr - spread ratio (mv optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(prmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(prmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(prmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()

#####
# plot results - tasmax

range(c(txmeanall,txmean_pro,txmean_txo,txmean_tno,txmean_mvo),na.rm=TRUE)

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmax_changeplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(txmeanall,txmean_pro,txmean_txo,txmean_tno,txmean_mvo),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,8))

testsfc1 = list(x=lonhist,y=lathist,z=txmeanall)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - full ensemble mean change",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_pro)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean change, pr optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_txo)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean change, tasmax optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_tno)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean change, tasmin optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=txmean_mvo)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean change, mv optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()


####
# plot results - tasmax error

range(c(txmeanall_error,txmean_pro_error,txmean_txo_error,txmean_tno_error,txmean_mvo_error),na.rm=TRUE)

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmax_errorplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(txmeanall_error,txmean_pro_error,txmean_txo_error,txmean_tno_error,txmean_mvo_error),colorchoice=colorchoicediff,Blimit=BINLIMIT,use_fixed_scale = TRUE,fixed_scale = c(-10,10))

testsfc1 = list(x=lonhist,y=lathist,z=txmeanall_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - full ensemble mean error",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_pro_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean error, pr optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_txo_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean error, tasmax optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_tno_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean error, tasmin optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=txmean_mvo_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - ensemble mean error, mv optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()

#####
# tasmax spread plots

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="ratio"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmax_spreadplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(txmean_pro_sr,txmean_txo_sr,txmean_tno_sr,txmean_mvo_sr),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,1))

testsfc1 = list(x=lonhist,y=lathist,z=txmean_pro_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - spread ratio (pr optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_txo_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - spread ratio (tasmax optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=txmean_tno_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - spread ratio (tasmin optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=txmean_mvo_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmax - spread ratio (mv optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(txmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(txmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(txmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()

#####
# plot results - tasmin

range(c(tnmeanall,tnmean_pro,tnmean_txo,tnmean_tno,tnmean_mvo),na.rm=TRUE)

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmin_changeplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(tnmeanall,tnmean_pro,tnmean_txo,tnmean_tno,tnmean_mvo),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,8))

testsfc1 = list(x=lonhist,y=lathist,z=tnmeanall)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - full ensemble mean change",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmeanall,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_pro)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean change, pr optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_pro,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_txo)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean change, tasmax optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_txo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_tno)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean change, tasmin optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_tno,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=tnmean_mvo)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean change, mv optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_mvo,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()


####
# plot results - tasmin error

range(c(tnmeanall_error,tnmean_pro_error,tnmean_txo_error,tnmean_tno_error,tnmean_mvo_error),na.rm=TRUE)

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="difference"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmin_errorplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(tnmeanall_error,tnmean_pro_error,tnmean_txo_error,tnmean_tno_error,tnmean_mvo_error),colorchoice=colorchoicediff,Blimit=BINLIMIT,use_fixed_scale = TRUE,fixed_scale = c(-10,10))

testsfc1 = list(x=lonhist,y=lathist,z=tnmeanall_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - full ensemble mean error",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmeanall_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_pro_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean error, pr optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_pro_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_txo_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean error, tasmax optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_txo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_tno_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean error, tasmin optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_tno_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=tnmean_mvo_error)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - ensemble mean error, mv optimal subset",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_mvo_error,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()

#####
# tasmin spread plots

colorchoicediff="bluetored"
BINLIMIT=30
diffbartype="ratio"
#use_fixed_scale=TRUE
#fixed_scale = c(-120,120)

plotfilename = paste("/home/woot0002/RCMES/tasmin_spreadplots_",ensgroup,"_dom",mask,"_2070-2099_size",enssize,"_shortEAA.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=W,height=H)

diffcolorbar = colorramp(c(tnmean_pro_sr,tnmean_txo_sr,tnmean_tno_sr,tnmean_mvo_sr),colorchoice=colorchoicediff,Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,1))

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_pro_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - spread ratio (pr optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_pro_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_txo_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - spread ratio (tasmax optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_txo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lonhist,y=lathist,z=tnmean_tno_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - spread ratio (tasmin optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_tno_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)


testsfc1 = list(x=lonhist,y=lathist,z=tnmean_mvo_sr)
surface(testsfc1,type="I",xlim=xlimin,ylim=ylimin,main="tasmin - spread ratio (mv optimal / full)",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
map("worldHires",add=TRUE)
text(textx,texty1,labels=paste("MAX = ",round(max(tnmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty2,labels=paste("MEAN = ",round(mean(tnmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(textx,texty3,labels=paste("MIN = ",round(min(tnmean_mvo_sr,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()



#####
# ggplot of domain results
# boxplot variations
datavals = c(prmeanall,prmean_pro,prmean_txo,prmean_tno,prmean_mvo)
varname = rep("pr",length(prmeanall)*5)
group = rep(c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),each = length(prmeanall))

prtab = data.frame(varname,group,datavals)
prtab$group = factor(prtab$group,levels=c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),ordered=TRUE)


datavals = c(txmeanall,txmean_pro,txmean_txo,txmean_tno,txmean_mvo)
varname = rep("tasmax",length(txmeanall)*5)
group = rep(c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),each = length(txmeanall))

txtab = data.frame(varname,group,datavals)
txtab$group = factor(txtab$group,levels=c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),ordered=TRUE)


datavals = c(tnmeanall,tnmean_pro,tnmean_txo,tnmean_tno,tnmean_mvo)
varname = rep("tasmin",length(txmeanall)*5)
group = rep(c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),each = length(tnmeanall))

tntab = data.frame(varname,group,datavals)
tntab$group = factor(tntab$group,levels=c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),ordered=TRUE)


prtab = prtab[-which(is.na(prtab$datavals)==TRUE),]
txtab = txtab[-which(is.na(txtab$datavals)==TRUE),]
tntab = tntab[-which(is.na(tntab$datavals)==TRUE),]

####

library(ggplot2)

ggplot(prtab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("pr Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Change (mm)")

ggplot(txtab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("tasmax Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Change (C)")

ggplot(tntab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("tasmin Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Change (C)")


#####
# ggplot of domain results
# boxplot variations - error
datavals = c(prmeanall_error,prmean_pro_error,prmean_txo_error,prmean_tno_error,prmean_mvo_error)
varname = rep("pr",length(prmeanall_error)*5)
group = rep(c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),each = length(prmeanall))

prtab = data.frame(varname,group,datavals)
prtab$group = factor(prtab$group,levels=c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),ordered=TRUE)


datavals = c(txmeanall_error,txmean_pro_error,txmean_txo_error,txmean_tno_error,txmean_mvo_error)
varname = rep("tasmax",length(txmeanall_error)*5)
group = rep(c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),each = length(txmeanall))

txtab = data.frame(varname,group,datavals)
txtab$group = factor(txtab$group,levels=c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),ordered=TRUE)


datavals = c(tnmeanall_error,tnmean_pro_error,tnmean_txo_error,tnmean_tno_error,tnmean_mvo_error)
varname = rep("tasmin",length(tnmeanall_error)*5)
group = rep(c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),each = length(tnmeanall))

tntab = data.frame(varname,group,datavals)
tntab$group = factor(tntab$group,levels=c("full ens","pr optim ens","tasmax optim ens","tasmin optim ens","mv optim ens"),ordered=TRUE)

prtab = prtab[-which(is.na(prtab$datavals)==TRUE),]
txtab = txtab[-which(is.na(txtab$datavals)==TRUE),]
tntab = tntab[-which(is.na(tntab$datavals)==TRUE),]

####

library(ggplot2)

ggplot(prtab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("pr Error Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Error (mm)")

ggplot(txtab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("tasmax Error Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Error (C)")

ggplot(tntab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("tasmin Error Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Error (C)")

#####
# ggplot of domain results
# boxplot variations - frc
datavals = c(prmean_pro_sr,prmean_r5o_sr,prmean_mvo_sr)
varname = rep("pr",length(prmean_pro_sr)*3)
group = rep(c("pr optim ens","rx5day optim ens","mv optim ens"),each = length(prmeanall))

prtab = data.frame(varname,group,datavals)
prtab$group = factor(prtab$group,levels=c("pr optim ens","rx5day optim ens","mv optim ens"),ordered=TRUE)

datavals = c(r5mean_pro_sr,r5mean_r5o_sr,r5mean_mvo_sr)
varname = rep("rx5day",length(r5mean_pro_sr)*3)
group = rep(c("pr optim ens","rx5day optim ens","mv optim ens"),each = length(prmeanall))

r5tab = data.frame(varname,group,datavals)
r5tab$group = factor(r5tab$group,levels=c("pr optim ens","rx5day optim ens","mv optim ens"),ordered=TRUE)

prtab = prtab[-which(is.na(prtab$datavals)==TRUE),]
r5tab = r5tab[-which(is.na(r5tab$datavals)==TRUE),]

taball = rbind(prtab,r5tab)

####

library(ggplot2)

ggplot(taball, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("spread Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Spread Ratio") +facet_wrap(facets=vars(varname))

ggplot(prtab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("pr spread Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Spread Ratio")+ylim(0,1)

ggplot(r5tab, aes(x=group, y=datavals)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("rx5day spread Comparisons - full ensemble vs. subsets size=",enssize,sep=""))+xlab("Ensemble")+ylab("Spread Ratio")+ylim(0,1)

######
# plot all members in ensemble and subset selected and full ensemble averages.

prchangeavg = apply(prchange,3,mean,na.rm=TRUE)
tasmaxchangeavg = apply(tasmaxchange,3,mean,na.rm=TRUE)
tasminchangeavg = apply(tasminchange,3,mean,na.rm=TRUE)

prmean_mvo_avg = apply(prchange[,,varidx_mvo],3,mean,na.rm=TRUE)
txmean_mvo_avg = apply(tasmaxchange[,,varidx_mvo],3,mean,na.rm=TRUE)
tnmean_mvo_avg = apply(tasminchange[,,varidx_mvo],3,mean,na.rm=TRUE)

fullmean_pr = mean(prmeanall,na.rm=TRUE)
fullmean_tx = mean(txmeanall,na.rm=TRUE)
fullmean_tn = mean(tnmeanall,na.rm=TRUE)

mvomean_pr = mean(prmean_mvo,na.rm=TRUE)
mvomean_tx = mean(txmean_mvo,na.rm=TRUE)
mvomean_tn = mean(tnmean_mvo,na.rm=TRUE)

prerrs = prhist
txerrs = tasmaxhist
tnerrs = tasminhist

for(i in 1:dim(prhist)[3]){
  prerrs[,,i]=prhist[,,i]-probs
  txerrs[,,i]=tasmaxhist[,,i]-tasmaxobs
  tnerrs[,,i]=tasminhist[,,i]-tasminobs
}

prerrsavg = apply(prerrs,3,mean,na.rm=TRUE)
txerrsavg = apply(txerrs,3,mean,na.rm=TRUE)
tnerrsavg = apply(tnerrs,3,mean,na.rm=TRUE)

prerrsavg_mvo = apply(prerrs[,,varidx_mvo],3,mean,na.rm=TRUE)
txerrsavg_mvo = apply(txerrs[,,varidx_mvo],3,mean,na.rm=TRUE)
tnerrsavg_mvo = apply(tnerrs[,,varidx_mvo],3,mean,na.rm=TRUE)

####

boxplot(prchangeavg,main="Pr projected change",ylab="mm")
points(fullmean_pr,col="red",pch=19)
points(mvomean_pr,col="purple",pch=19)
abline(h=0,lty=3)
legend("topleft",legend=c("Full","subset"),pch=19,col=c("red","purple"))

boxplot(prchangeavg,main="Pr projected change",ylab="mm")
for(j in 1:5){
  points(prmean_mvo_avg[j],col="red",pch=19)
}
abline(h=0,lty=3)
legend("topleft",legend=c("subset models "),pch=19,col=c("red"),cex=0.7)

###

boxplot(tasmaxchangeavg,main="tasmax projected change",ylab="C")
points(fullmean_tx,col="red",pch=19)
points(mvomean_tx,col="purple",pch=19)
abline(h=0,lty=3)
legend("topleft",legend=c("Full","subset"),pch=19,col=c("red","purple"))


boxplot(tasmaxchangeavg,main="tasmax projected change",ylab="C")
for(j in 1:5){
  points(txmean_mvo_avg[j],col="red",pch=19)
}
abline(h=0,lty=3)
legend("topleft",legend=c("subset models "),pch=19,col=c("red"),cex=0.7)

###

boxplot(tasminchangeavg,main="tasmin projected change",ylab="C")
points(fullmean_tn,col="red",pch=19)
points(mvomean_tn,col="purple",pch=19)
abline(h=0,lty=3)
legend("topleft",legend=c("Full","subset"),pch=19,col=c("red","purple"))

boxplot(tasmaxchangeavg,main="tasmin projected change",ylab="C")
for(j in 1:5){
  points(tnmean_mvo_avg[j],col="red",pch=19)
}
abline(h=0,lty=3)
legend("topleft",legend=c("subset models "),pch=19,col=c("red"),cex=0.7)

####

boxplot(prerrsavg,main="Pr Error",ylab="mm")
points(mean(prerrsavg),col="red",pch=19)
points(mean(prerrsavg_mvo),col="purple",pch=19)
abline(h=0,lty=3)
legend("topleft",legend=c("Full","subset"),pch=19,col=c("red","purple"))

boxplot(prerrsavg,main="Pr Error",ylab="mm")
for(j in 1:5){
  points(prerrsavg_mvo[j],col="red",pch=19)
}
abline(h=0,lty=3)
legend("topleft",legend=c("subset models "),pch=19,col=c("red"),cex=0.7)

####

boxplot(txerrsavg,main="Tasmax Error",ylab="C")
points(mean(txerrsavg),col="red",pch=19)
points(mean(txerrsavg_mvo),col="purple",pch=19)
abline(h=0,lty=3)
legend("topleft",legend=c("Full","subset"),pch=19,col=c("red","purple"))

boxplot(txerrsavg,main="Tasmax Error",ylab="C")
for(j in 1:5){
  points(txerrsavg_mvo[j],col="red",pch=19)
}
abline(h=0,lty=3)
legend("topleft",legend=c("subset models "),pch=19,col=c("red"),cex=0.7)

####

boxplot(tnerrsavg,main="Tasmin Error",ylab="C")
points(mean(tnerrsavg),col="red",pch=19)
points(mean(tnerrsavg_mvo),col="purple",pch=19)
abline(h=0,lty=3)
legend("topleft",legend=c("Full","subset"),pch=19,col=c("red","purple"))

boxplot(tnerrsavg,main="Tasmin Error",ylab="C")
for(j in 1:5){
  points(tnerrsavg_mvo[j],col="red",pch=19)
}
abline(h=0,lty=3)
legend("topleft",legend=c("subset models "),pch=19,col=c("red"),cex=0.7)


