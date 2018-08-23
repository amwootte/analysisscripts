#####################
#
# Boxplot creation - historical data

library(sp)
library(fields)
library(maps)
library(raster)
library(rasterVis)
library(maptools)
setwd("/home/woot0002/3to5/")

#####
# obs grab

load("Obs3to5climo.Rdata")
load("Mask3to5.Rdata")

prcptot = rnnmm = rx1day = c()

for(i in 1:3){
  prcptot = c(prcptot,ifelse(mask==1,obsprcptot[[i]],NA))
  rnnmm = c(rnnmm,ifelse(mask==1,obsrnnmm[[i]],NA))
  rx1day = c(rx1day,ifelse(mask==1,obsrx1day[[i]],NA))
}

models = obs
model = rep(models,each=length(rx1day)/3)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = outputframe[-NAidx,]

#####
# GCM grab

load("MH3to5climo.Rdata")
load("Domain3to5.Rdata")

prcptot = rnnmm = rx1day = c()

for(i in 1:3){
  prcptot = c(prcptot,ifelse(mask==1,MHprcptot[[i]],NA))
  rnnmm = c(rnnmm,ifelse(mask==1,MHrnnmm[[i]],NA))
  rx1day = c(rx1day,ifelse(mask==1,MHrx1day[[i]],NA))
}

models = GCMs
model = rep(models,each=length(rx1day)/3)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

#####
# LOCA grab

load("LOCAH3to5climo.Rdata")

test = readShapePoly("/home/woot0002/Shape/cb_2017_us_nation_5m") # appropriate filename format
projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

prcptot = rnnmm = rx1day = c()

for(i in 1:3){

tmpV = LOCAHprcptot[[i]] # need the 
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
LOCAHprcptot[[i]]=testin
prcptot = c(prcptot,testin)

tmpV = LOCAHrnnmm[[i]] # need the 
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
LOCAHrnnmm[[i]]=testin
rnnmm = c(rnnmm,testin)

tmpV = LOCAHrx1day[[i]] # need the 
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
LOCAHrx1day[[i]]=testin
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

prcptot = rnnmm = rx1day = c()

for(i in 1:2){
  prcptot = c(prcptot,MACAHprcptot[[i]])
  rnnmm = c(rnnmm,MACAHrnnmm[[i]])
  rx1day = c(rx1day,MACAHrx1day[[i]])
}

models = paste("MACA_",GCMs,sep="")
model = rep(models,each=length(rx1day)/2)
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

###
# QDM grab

load("QDMH3to5climo.Rdata")
QDMHorder2 = c("CCSM4_daymet","MIROC5_daymet","MPI-ESM-LR_daymet","CCSM4_livneh","MIROC5_livneh","MPI-ESM-LR_livneh","CCSM4_prism","MIROC5_prism","MPI-ESM-LR_prism")
prcptot = rnnmm = rx1day = c()

for(i in 1:9){
  idx = which(QDMHorder == QDMHorder2[i])
  prcptot = c(prcptot,QDMHprcptot[[idx]])
  rnnmm = c(rnnmm,QDMHrnnmm[[idx]])
  rx1day = c(rx1day,QDMHrx1day[[idx]])
}

models = paste("QDM",QDMHorder2,sep="_")
model = as.character(rep(models,each=length(rx1day)/9))
outputframe = data.frame(model,prcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

#######
# Create boxplots

#boxplot(prcptot~model,data=alldat,at=c(1,2,3,5,6,7,9,10,11,13,14,16,17,18),xaxt="n",col=c("red","purple","magenta","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","brown","cyan","yellow"))
#text(x =  c(1,2,3,5,6,7,9,10,11,13,14,16,17,18), y = par("usr")[3] - 1, srt = 45, adj = 1,
#     labels = unique(alldat$model), xpd = TRUE,cex=0.75)
#boxplot(prcptot~model,data=alldat,at=c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22),xaxt="n",col=c("red","purple","magenta","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","brown","cyan","yellow","brown","cyan","yellow"))
#text(x =  c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22), y = par("usr")[3] - 1, srt = 45, adj = 1,
#     labels = unique(alldat$model), xpd = TRUE,cex=0.75)

alldat$plotorder = NA
alldat$plotorder[which(alldat$model=="daymet")]=1
alldat$plotorder[which(alldat$model=="livneh")]=2
alldat$plotorder[which(alldat$model=="prism")]=3

alldat$plotorder[which(alldat$model=="CCSM4")]=4
alldat$plotorder[which(alldat$model=="MIROC5")]=5
alldat$plotorder[which(alldat$model=="MPI-ESM-LR")]=6

alldat$plotorder[which(alldat$model=="LOCA_CCSM4")]=7
alldat$plotorder[which(alldat$model=="LOCA_MIROC5")]=8
alldat$plotorder[which(alldat$model=="LOCA_MPI-ESM-LR")]=9

alldat$plotorder[which(alldat$model=="MACA_CCSM4")]=10
alldat$plotorder[which(alldat$model=="MACA_MIROC5")]=11

alldat$plotorder[which(alldat$model=="QDM_CCSM4_daymet")]=12
alldat$plotorder[which(alldat$model=="QDM_MIROC5_daymet")]=13
alldat$plotorder[which(alldat$model=="QDM_MPI-ESM-LR_daymet")]=14

alldat$plotorder[which(alldat$model=="QDM_CCSM4_livneh")]=15
alldat$plotorder[which(alldat$model=="QDM_MIROC5_livneh")]=16
alldat$plotorder[which(alldat$model=="QDM_MPI-ESM-LR_livneh")]=17

alldat$plotorder[which(alldat$model=="QDM_CCSM4_prism")]=18
alldat$plotorder[which(alldat$model=="QDM_MIROC5_prism")]=19
alldat$plotorder[which(alldat$model=="QDM_MPI-ESM-LR_prism")]=20

boxplot(prcptot~plotorder,data=alldat,main="prcptot Obs, GCM, and DS Values across I35 region grid points",ylab="mm",at=c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22,24,25,26),xaxt="n",col=c("red","purple","magenta","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","yellow"))
text(x =  c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22,24,25,26), y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = unique(alldat$model), xpd = TRUE,cex=0.5)

boxplot(rnnmm~plotorder,data=alldat,main="rnnmm Obs, GCM, and DS Values across I35 region grid points",ylab="days",at=c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22,24,25,26),xaxt="n",col=c("red","purple","magenta","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","yellow"))
text(x =  c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22,24,25,26), y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = unique(alldat$model), xpd = TRUE,cex=0.5)

boxplot(rx1day~plotorder,data=alldat,main="rx1day Obs, GCM, and DS Values across I35 region grid points",ylab="mm",at=c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22,24,25,26),xaxt="n",col=c("red","purple","magenta","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","brown","cyan","yellow","brown","cyan","yellow","brown","cyan","yellow"))
text(x =c(1,2,3,5,6,7,9,10,11,13,14,16,17,18,20,21,22,24,25,26), y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = unique(alldat$model), xpd = TRUE,cex=0.5)


sapply(MHprcptot,mean,na.rm=TRUE)



