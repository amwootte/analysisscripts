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

setwd("/home/woot0002/DS_ind/")

load("valuesarrays_WUpr.Rdata")

GCM_bias_pr_WUpr = GCM_bias_pr
LOCA_bias_pr_WUpr = LOCA_bias_pr
GCM_change_pr_WUpr = GCM_change_pr
LOCA_change_pr_WUpr = LOCA_change_pr

GCM_bias_tmax_WUpr = GCM_bias_tmax
LOCA_bias_tmax_WUpr = LOCA_bias_tmax
GCM_change_tmax_WUpr = GCM_change_tmax
LOCA_change_tmax_WUpr = LOCA_change_tmax

load("valuesarrays_WUtmax.Rdata")

GCM_bias_pr_WUtmax = GCM_bias_pr
LOCA_bias_pr_WUtmax = LOCA_bias_pr
GCM_change_pr_WUtmax = GCM_change_pr
LOCA_change_pr_WUtmax = LOCA_change_pr

GCM_bias_tmax_WUtmax = GCM_bias_tmax
LOCA_bias_tmax_WUtmax = LOCA_bias_tmax
GCM_change_tmax_WUtmax = GCM_change_tmax
LOCA_change_tmax_WUtmax = LOCA_change_tmax

#####
# Create color bars

biasprbar = colorramp(c(GCM_bias_pr_WUpr,LOCA_bias_pr_WUpr,GCM_bias_pr_WUtmax,LOCA_bias_pr_WUtmax),colorchoice="bluetored",Blimit=30,type="difference",use_fixed_scale = TRUE,fixed_scale = c(-1200,1200))
biastmaxbar = colorramp(c(GCM_bias_tmax_WUpr,LOCA_bias_tmax_WUpr,GCM_bias_tmax_WUtmax,LOCA_bias_tmax_WUtmax),colorchoice="bluetored",Blimit=20,type="difference",use_fixed_scale = FALSE)

changeprbar = colorramp(c(GCM_change_pr_WUpr,LOCA_change_pr_WUpr,GCM_change_pr_WUtmax,LOCA_change_pr_WUtmax),colorchoice="browntogreen",Blimit=30,type="difference",use_fixed_scale = TRUE, fixed_scale = c(-250,250))
changetmaxbar = colorramp(c(GCM_change_tmax_WUpr,LOCA_change_tmax_WUpr,GCM_change_tmax_WUtmax,LOCA_change_tmax_WUtmax),colorchoice="bluetored",Blimit=20,type="difference",use_fixed_scale = TRUE,fixed_scale = c(0,8))

biastmaxbar_GCM = colorramp(c(GCM_bias_tmax_WUpr,GCM_bias_tmax_WUtmax),colorchoice="bluetored",Blimit=20,type="difference",use_fixed_scale = FALSE)

biasprbar_GCM = colorramp(c(GCM_bias_pr_WUpr,GCM_bias_pr_WUtmax),colorchoice="redtoblue",Blimit=30,type="difference",use_fixed_scale = FALSE)

#####
# value limiter

LOCA_bias_pr_WUpr2 = ifelse(LOCA_bias_pr_WUpr>250,250,LOCA_bias_pr_WUpr)
LOCA_bias_pr_WUpr2 = ifelse(LOCA_bias_pr_WUpr2< -250,-250,LOCA_bias_pr_WUpr2)

LOCA_bias_pr_WUtmax2 = ifelse(LOCA_bias_pr_WUtmax >250,250,LOCA_bias_pr_WUtmax)
LOCA_bias_pr_WUtmax2 = ifelse(LOCA_bias_pr_WUtmax2< -250,-250,LOCA_bias_pr_WUtmax2)

biasprbar_LOCA = colorramp(c(LOCA_bias_pr_WUpr2,LOCA_bias_pr_WUtmax2),colorchoice="redtoblue",Blimit=20,type="difference",use_fixed_scale = TRUE,fixed_scale=c(-250,250))


LOCA_bias_tmax_WUpr2 = ifelse(LOCA_bias_tmax_WUpr>2,2,LOCA_bias_tmax_WUpr)
LOCA_bias_tmax_WUpr2 = ifelse(LOCA_bias_tmax_WUpr2< -2,-2,LOCA_bias_tmax_WUpr2)

LOCA_bias_tmax_WUtmax2 = ifelse(LOCA_bias_tmax_WUtmax >2,2,LOCA_bias_tmax_WUtmax)
LOCA_bias_tmax_WUtmax2 = ifelse(LOCA_bias_tmax_WUtmax2< -2,-2,LOCA_bias_tmax_WUtmax2)

biastmaxbar_LOCA = colorramp(c(LOCA_bias_tmax_WUpr2,LOCA_bias_tmax_WUtmax2),colorchoice="bluetored",Blimit=20,type="difference",use_fixed_scale = TRUE,fixed_scale=c(-2,2))



#####
# Create plots - pr displayed

for(w in 1:3){
for(s in 1:3){
if(s==1){
  XLIM = range(lon)
  YLIM = range(lat)
  textlabx = -109 
  textlaby1 = 28.45
  textlaby2 = 27.45
  textlaby3 = 26.45
}

if(s==2){
  XLIM = c(-94.5,-89.5)
  YLIM = c(28,34)
  textlabx = -94 
  textlaby1 = 29.25
  textlaby2 = 29
  textlaby3 = 28.75
}

if(s==3){
  XLIM = c(-109.25,-102.75)
  YLIM = c(29,37.5)
  textlabx = -109 
  textlaby1 = 30.25
  textlaby2 = 30
  textlaby3 = 29.75
}

schemes = c("unweighted","skill","SI-h","SI-c","BMA")
pdf(paste("GCMprmaps_WU",WUs[w],"_SA",SAs[s],"_customscale2.pdf",sep=""),width=25,height=20,compress=TRUE)
par(mfrow=c(4,5))
for(i in 1:20){
  if(i>=1 & i<=5){
    testsfc = list(x=lon,y=lat,z=GCM_bias_pr_WUtmax[,,i,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[i]," pr bias, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biasprbar_GCM[[1]],col=biasprbar_GCM[[3]],breaks=biasprbar_GCM[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_bias_pr_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_bias_pr_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_bias_pr_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=6 & i<=10){
    I=i-5
    testsfc = list(x=lon,y=lat,z=GCM_change_pr_WUtmax[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," pr change, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changeprbar[[1]],col=changeprbar[[3]],breaks=changeprbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_change_pr_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_change_pr_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_change_pr_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i==11 | i==16){
    frame()
  }
  
  if(i>=12 & i<=15){
    I=i-10
    testsfc = list(x=lon,y=lat,z=GCM_bias_pr_WUpr[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," pr bias, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biasprbar_GCM[[1]],col=biasprbar_GCM[[3]],breaks=biasprbar_GCM[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_bias_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_bias_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_bias_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=17 & i<=20){
    I=i-15
    testsfc = list(x=lon,y=lat,z=GCM_change_pr_WUpr[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," pr change, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changeprbar[[1]],col=changeprbar[[3]],breaks=changeprbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_change_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_change_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_change_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  
}
dev.off()



pdf(paste("LOCAprmaps_WU",WUs[w],"_SA",SAs[s],"_customscale2.pdf",sep=""),width=25,height=20)
par(mfrow=c(4,5))
for(i in 1:20){
  if(i>=1 & i<=5){
    testsfc = list(x=lon,y=lat,z=LOCA_bias_pr_WUtmax2[,,i,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[i]," pr bias, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biasprbar_LOCA[[1]],col=biasprbar_LOCA[[3]],breaks=biasprbar_LOCA[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_bias_pr_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_bias_pr_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_bias_pr_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=6 & i<=10){
    I=i-5
    testsfc = list(x=lon,y=lat,z=LOCA_change_pr_WUtmax[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," pr change, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changeprbar[[1]],col=changeprbar[[3]],breaks=changeprbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_change_pr_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_change_pr_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_change_pr_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i==11 | i==16){
    frame()
  }
  
  if(i>=12 & i<=15){
    I=i-10
    testsfc = list(x=lon,y=lat,z=LOCA_bias_pr_WUpr2[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," pr bias, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biasprbar_LOCA[[1]],col=biasprbar_LOCA[[3]],breaks=biasprbar_LOCA[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_bias_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_bias_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_bias_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=17 & i<=20){
    I=i-15
    testsfc = list(x=lon,y=lat,z=LOCA_change_pr_WUpr[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," pr change, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changeprbar[[1]],col=changeprbar[[3]],breaks=changeprbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_change_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_change_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_change_pr_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  
}
dev.off()

###########

pdf(paste("GCMtmaxmaps_WU",WUs[w],"_SA",SAs[s],"_customscale.pdf",sep=""),width=25,height=20)
par(mfrow=c(4,5))
for(i in 1:20){
  if(i>=1 & i<=5){
    testsfc = list(x=lon,y=lat,z=GCM_bias_tmax_WUtmax[,,i,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[i]," tmax bias, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biastmaxbar_GCM[[1]],col=biastmaxbar_GCM[[3]],breaks=biastmaxbar_GCM[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_bias_tmax_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_bias_tmax_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_bias_tmax_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=6 & i<=10){
    I=i-5
    testsfc = list(x=lon,y=lat,z=GCM_change_tmax_WUtmax[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," tmax change, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changetmaxbar[[1]],col=changetmaxbar[[3]],breaks=changetmaxbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_change_tmax_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_change_tmax_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_change_tmax_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i==11 | i==16){
    frame()
  }
  
  if(i>=12 & i<=15){
    I=i-10
    testsfc = list(x=lon,y=lat,z=GCM_bias_tmax_WUpr[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," tmax bias, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biastmaxbar_GCM[[1]],col=biastmaxbar_GCM[[3]],breaks=biastmaxbar_GCM[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_bias_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_bias_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_bias_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=17 & i<=20){
    I=i-15
    testsfc = list(x=lon,y=lat,z=GCM_change_tmax_WUpr[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," tmax change, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changetmaxbar[[1]],col=changetmaxbar[[3]],breaks=changetmaxbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(GCM_change_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(GCM_change_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(GCM_change_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  
}
dev.off()



pdf(paste("LOCAtmaxmaps_WU",WUs[w],"_SA",SAs[s],"_customscale.pdf",sep=""),width=25,height=20)
par(mfrow=c(4,5))
for(i in 1:20){
  if(i>=1 & i<=5){
    testsfc = list(x=lon,y=lat,z=LOCA_bias_tmax_WUtmax2[,,i,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[i]," pr bias, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biastmaxbar_LOCA[[1]],col=biastmaxbar_LOCA[[3]],breaks=biastmaxbar_LOCA[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_bias_tmax_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_bias_tmax_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_bias_tmax_WUtmax[,,i,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=6 & i<=10){
    I=i-5
    testsfc = list(x=lon,y=lat,z=LOCA_change_tmax_WUtmax[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," tmax change, tmax weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changetmaxbar[[1]],col=changetmaxbar[[3]],breaks=changetmaxbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_change_tmax_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_change_tmax_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_change_tmax_WUtmax[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i==11 | i==16){
    frame()
  }
  
  if(i>=12 & i<=15){
    I=i-10
    testsfc = list(x=lon,y=lat,z=LOCA_bias_tmax_WUpr2[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," tmax bias, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=biastmaxbar_LOCA[[1]],col=biastmaxbar_LOCA[[3]],breaks=biastmaxbar_LOCA[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_bias_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_bias_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_bias_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(i>=17 & i<=20){
    I=i-15
    testsfc = list(x=lon,y=lat,z=LOCA_change_tmax_WUpr[,,I,w,s])   
    surface(testsfc,type="I",xlim=XLIM,ylim=YLIM,main=paste(schemes[I]," tmax change, pr weighted",sep=""),xlab="Longtitude",ylab="Latitude",zlim=changetmaxbar[[1]],col=changetmaxbar[[3]],breaks=changetmaxbar[[2]])
    map("state",add=TRUE)
    text(textlabx,textlaby1,labels=paste("MAX = ",round(max(LOCA_change_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby2,labels=paste("MEAN = ",round(mean(LOCA_change_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(textlabx,textlaby3,labels=paste("MIN = ",round(min(LOCA_change_tmax_WUpr[,,I,w,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  
}
dev.off()

}
}











