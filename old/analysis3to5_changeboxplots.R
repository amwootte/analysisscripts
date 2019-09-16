#####################
#
# Boxplot creation - projected change data

library(sp)
library(fields)
library(maps)
library(raster)
library(rasterVis)
library(maptools)
setwd("/home/woot0002/3to5/")

#####
# obs grab

load("MH3to5climo.Rdata")
load("MF3to5climo.Rdata")
load("Mask3to5.Rdata")

prcptot = rnnmm = rx1day = c()

for(i in 1:3){
  prcptot = c(prcptot,ifelse(mask==1,MFprcptot[[i]],NA)-ifelse(mask==1,MHprcptot[[i]],NA))
  rnnmm = c(rnnmm,ifelse(mask==1,MFrnnmm[[i]],NA)-ifelse(mask==1,MHrnnmm[[i]],NA))
  rx1day = c(rx1day,ifelse(mask==1,MFrx1day[[i]],NA)-ifelse(mask==1,MHrx1day[[i]],NA))
}

models = GCMs
model = rep(models,each=length(rx1day)/3)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = outputframe[-NAidx,]

###
# LOCA grab

load("LOCAH3to5climo.Rdata")
load("LOCAF3to5climo.Rdata")

test = readShapePoly("/home/woot0002/Shape/cb_2017_us_nation_5m") # appropriate filename format
projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

prcptot = rnnmm = rx1day = c()

for(i in 1:3){
  
  tmpV = LOCAFprcptot[[i]]-LOCAHprcptot[[i]] # need the 
  modrasV = raster(t(tmpV)[length(latLOCA):1,])
  if(all(lonLOCA>0)){
    extent(modrasV) = c(min(lonLOCA)-360,max(lonLOCA)-360,min(latLOCA),max(latLOCA))
  } else {
    extent(modrasV) = c(min(lonLOCA),max(lonLOCA),min(latLOCA),max(latLOCA))
  }
  mod.subV <- crop(modrasV, extent(test))
  mod.subV <- mask(modrasV, test)
  testin = matrix(getValues(mod.subV),nrow=length(lonLOCA),ncol=length(latLOCA))
  testin = testin[,length(latLOCA):1]
  prcptot = c(prcptot,testin)
  
  tmpV = LOCAFrnnmm[[i]]-LOCAHrnnmm[[i]] # need the 
  modrasV = raster(t(tmpV)[length(latLOCA):1,])
  if(all(lonLOCA>0)){
    extent(modrasV) = c(min(lonLOCA)-360,max(lonLOCA)-360,min(latLOCA),max(latLOCA))
  } else {
    extent(modrasV) = c(min(lonLOCA),max(lonLOCA),min(latLOCA),max(latLOCA))
  }
  mod.subV <- crop(modrasV, extent(test))
  mod.subV <- mask(modrasV, test)
  testin = matrix(getValues(mod.subV),nrow=length(lonLOCA),ncol=length(latLOCA))
  testin = testin[,length(latLOCA):1]
  rnnmm = c(rnnmm,testin)
  
  tmpV = LOCAFrx1day[[i]]-LOCAHrx1day[[i]] # need the 
  modrasV = raster(t(tmpV)[length(latLOCA):1,])
  if(all(lonLOCA>0)){
    extent(modrasV) = c(min(lonLOCA)-360,max(lonLOCA)-360,min(latLOCA),max(latLOCA))
  } else {
    extent(modrasV) = c(min(lonLOCA),max(lonLOCA),min(latLOCA),max(latLOCA))
  }
  mod.subV <- crop(modrasV, extent(test))
  mod.subV <- mask(modrasV, test)
  testin = matrix(getValues(mod.subV),nrow=length(lonLOCA),ncol=length(latLOCA))
  testin = testin[,length(latLOCA):1]
  rx1day = c(rx1day,testin)
}

models = paste("LOCA_",GCMs,sep="")
model = rep(models,each=length(rx1day)/3)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

####
# MACA grab

load("MACAH3to5climo.Rdata")
load("MACAF3to5climo.Rdata")

prcptot = rnnmm = rx1day = c()

for(i in 1:2){
  prcptot = c(prcptot,MACAFprcptot[[i]]-MACAHprcptot[[i]])
  rnnmm = c(rnnmm,MACAFrnnmm[[i]]-MACAHrnnmm[[i]])
  rx1day = c(rx1day,MACAFrx1day[[i]]-MACAHrx1day[[i]])
}

models = paste("MACA_",GCMs,sep="")
model = rep(models,each=length(rx1day)/2)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

###
# QDM grab

load("QDMH3to5climo.Rdata")
load("QDMF3to5climo.Rdata")
prcptot = rnnmm = rx1day = c()

for(i in 1:9){
  prcptot = c(prcptot,QDMFprcptot[[i]]-QDMHprcptot[[i]])
  rnnmm = c(rnnmm,QDMFrnnmm[[i]]-QDMHrnnmm[[i]])
  rx1day = c(rx1day,QDMFrx1day[[i]]-QDMHrx1day[[i]])
}

models = paste("QDM",QDMForder,sep="_")
model = rep(models,each=length(rx1day)/9)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

###
# DeltaSD grab

load("Obs3to5climo.Rdata")
load("DeltaSDF3to5climo.Rdata")

prcptot = rnnmm = rx1day = c()

for(i in 1:9){
  if(i==1 | i==4 | i==7) obsidx=1
  if(i==2 | i==5 | i==8) obsidx=2
  if(i==3 | i==6 | i==9) obsidx=3
  prcptot = c(prcptot,DeltaSDFprcptot[[i]]-obsprcptot[[obsidx]])
  rnnmm = c(rnnmm,DeltaSDFrnnmm[[i]]-obsrnnmm[[obsidx]])
  rx1day = c(rx1day,DeltaSDFrx1day[[i]]-obsrx1day[[obsidx]])
}

models = paste("DeltaSD",DeltaSDForder,sep="_")
model = rep(models,each=length(rx1day)/9)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

#######
# Create boxplots

alldat$plotorder = NA

alldat$plotorder[which(alldat$model=="CCSM4")]=1
alldat$plotorder[which(alldat$model=="MIROC5")]=2
alldat$plotorder[which(alldat$model=="MPI-ESM-LR")]=3

alldat$plotorder[which(alldat$model=="DeltaSD_CCSM4_daymet")]=4
alldat$plotorder[which(alldat$model=="DeltaSD_MIROC5_daymet")]=5
alldat$plotorder[which(alldat$model=="DeltaSD_MPI-ESM-LR_daymet")]=6

alldat$plotorder[which(alldat$model=="QDM_CCSM4_daymet")]=7
alldat$plotorder[which(alldat$model=="QDM_MIROC5_daymet")]=8
alldat$plotorder[which(alldat$model=="QDM_MPI-ESM-LR_daymet")]=9

alldat$plotorder[which(alldat$model=="DeltaSD_CCSM4_livneh")]=10
alldat$plotorder[which(alldat$model=="DeltaSD_MIROC5_livneh")]=11
alldat$plotorder[which(alldat$model=="DeltaSD_MPI-ESM-LR_livneh")]=12

alldat$plotorder[which(alldat$model=="LOCA_CCSM4")]=13
alldat$plotorder[which(alldat$model=="LOCA_MIROC5")]=14
alldat$plotorder[which(alldat$model=="LOCA_MPI-ESM-LR")]=15

alldat$plotorder[which(alldat$model=="MACA_CCSM4")]=16
alldat$plotorder[which(alldat$model=="MACA_MIROC5")]=17

alldat$plotorder[which(alldat$model=="QDM_CCSM4_livneh")]=18
alldat$plotorder[which(alldat$model=="QDM_MIROC5_livneh")]=19
alldat$plotorder[which(alldat$model=="QDM_MPI-ESM-LR_livneh")]=20

alldat$plotorder[which(alldat$model=="DeltaSD_CCSM4_prism")]=21
alldat$plotorder[which(alldat$model=="DeltaSD_MIROC5_prism")]=22
alldat$plotorder[which(alldat$model=="DeltaSD_MPI-ESM-LR_prism")]=23

alldat$plotorder[which(alldat$model=="QDM_CCSM4_prism")]=24
alldat$plotorder[which(alldat$model=="QDM_MIROC5_prism")]=25
alldat$plotorder[which(alldat$model=="QDM_MPI-ESM-LR_prism")]=26

labelsused = c("CCSM4","MIROC5","MPI-ESM-LR",
               "DeltaSD_CCSM4_daymet","DeltaSD_MIROC5_daymet","DeltaSD_MPI-ESM-LR_daymet",
               "QDM_CCSM4_daymet","QDM_MIROC5_daymet","QDM_MPI-ESM-LR_daymet",
               "DeltaSD_CCSM4_livneh","DeltaSD_MIROC5_livneh","DeltaSD_MPI-ESM-LR_livneh",
               "LOCA_CCSM4","LOCA_MIROC5","LOCA_MPI-ESM-LR",
               "MACA_CCSM4","MACA_MIROC5",
               "QDM_CCSM4_livneh","QDM_MIROC5_livneh","QDM_MPI-ESM-LR_livneh",
               "DeltaSD_CCSM4_prism","DeltaSD_MIROC5_prism","DeltaSD_MPI-ESM-LR_prism",
               "QDM_CCSM4_prism","QDM_MIROC5_prism","QDM_MPI-ESM-LR_prism")

locations = c(1,2,3,5,6,7,9,10,11,13,14,15,17,18,19,21,22,24,25,26,28,29,30,32,33,34)
colors=c(rep(c("brown","cyan","yellow"),5),"brown","cyan",rep(c("brown","cyan","yellow"),3))

boxplot(prcptot~plotorder,data=alldat,main="prcptot 2070-2099 minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",ylab="mm",at=locations,xaxt="n",col=colors)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsused, xpd = TRUE,cex=0.5)

boxplot(rnnmm~plotorder,data=alldat,main="rnnmm 2070-2099 minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",ylab="days",at=locations,xaxt="n",col=colors)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsused, xpd = TRUE,cex=0.5)

boxplot(rx1day~plotorder,data=alldat,main="rx1day 2070-2099 minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",ylab="mm",at=locations,xaxt="n",col=colors)
text(x =locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsused, xpd = TRUE,cex=0.5)


