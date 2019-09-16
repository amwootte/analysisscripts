##############
#
# Sanderson weight calculator
# plus consistency

RMSE = function(exp,obs){
  err = exp-obs
  errsq = err^2
  merr = mean(errsq,na.rm=TRUE)
  if(merr!=0){
    rmse = sqrt(merr)
  } else {
    rmse =0
  }
  rmse
}

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

load("/home/woot0002/RMSEmats_combo.Rdata")
#load("/home/woot0002/RMSEfiles/tasmin_RMSEmats.Rdata")

numobs = 5
numtotal = nrow(normRMSEmat)
nummodels = numtotal-numobs
obsidx = (nummodels+1):numtotal
obsnames = simnames[obsidx]

###
# skill check - each variable independently - tasmax, tasmin, pr

Duvec = seq(0.05,2,by=0.01)
Dc = 0.7
Dq = 0.7

##
# similarity weight calculation

Wumat = matrix(NA,nrow=nummodels,ncol=length(Duvec))
for(u in 1:length(Duvec)){

  S_top = -(normRMSEmat/Duvec[u])^2
  S = exp(S_top)

  Ru = c()
  for(i in 1:nummodels){
    Ru[i]=1+sum(S[1:nummodels,i],na.rm=TRUE)
  }

Wumat[,u] = Ru^(-1)
}

##
# quality weight calculation

Wq =c()
t_dist = c()
for(i in 1:nummodels){
  if(metadat$training[i]=="Daymet"){
    useidx = obsidx[1]
  }
  if(metadat$training[i]=="Livneh"){
    useidx = obsidx[2]
  }
  if(metadat$training[i]=="PRISM"){
    useidx = obsidx[3]
  }
  if(metadat$training[i]=="METDATA"){
    useidx = obsidx[4]
  }
  if(metadat$training[i]=="MAURER"){
    useidx = obsidx[5]
  }
  
  t_dist[i] = normRMSEmat[useidx,i]
  
  Q_top = -(normRMSEmat[useidx,i]/Dq)^2
  Wq[i] = exp(Q_top)
}

##
# consitency weight calculation
Wc = c()
c_dist = c()
for(i in 1:nummodels){
  if(metadat$training[i]=="Daymet"){
    useidx = obsidx[2:5]
  }
  if(metadat$training[i]=="Livneh"){
    useidx = obsidx[c(1,3:5)]
  }
  if(metadat$training[i]=="PRISM"){
    useidx = obsidx[c(1:2,4:5)]
  }
  if(metadat$training[i]=="METDATA"){
    useidx = obsidx[c(1:3,5)]
  }
  if(metadat$training[i]=="MAURER"){
    useidx = obsidx[1:4]
  }
  
  c_dist[i]= mean(normRMSEmat[useidx,i],na.rm=TRUE)
  
  C_top = -(mean(normRMSEmat[useidx,i],na.rm=TRUE)/Dc)^2
  Wc[i] = exp(C_top)
}

########
# Calculate final weight

W = WA = array(NA,dim=c(nummodels,length(Duvec)))

for(j in 1:length(Duvec)){
  W[,j] = Wumat[,j]*Wc*Wq
  WA[,j] = W[,j]/sum(W[,j])
}

#######
# get climos, weighted ensemble mean
var = "tasmax"

EDQMfiles = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*EDQM*historical*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",var,"_*histclimo.nc",sep=""),intern=TRUE)
MACALfiles = system(paste("ls /home/woot0002/MACA_LIVNEH/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)
MACAMfiles = system(paste("ls /home/woot0002/MACA_METDATA/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)
BCCAfiles = system(paste("ls /home/woot0002/BCCA/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)

DAYMETfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*daymet*.nc",sep=""),intern=TRUE)
LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*livneh*.nc",sep=""),intern=TRUE)
PRISMfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*prism*.nc",sep=""),intern=TRUE)
METDATAfile = system(paste("ls /data2/3to5/I35/METDATA/",var,"_histclimo_mon.nc",sep=""),intern=TRUE)
MAURERfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_maurer_histclimo_mon.nc",sep=""),intern=TRUE)

###
# create full filelist + metadata table

#EDQM
filelist1 = do.call("rbind",strsplit(EDQMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,5],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$DS = "EDQM"
filelist2$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filelist2$training = rep(c("Daymet","Livneh","PRISM"),3)
filelist2$exp = filelist2[,5]
EDQMdat = filelist2[,c(10,12,9,11)]

#LOCA
filelist1 = do.call("rbind",strsplit(LOCAfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCAdat = filelist2[,c(2,3,4,6)]
names(LOCAdat) = names(EDQMdat)

#MACA Livneh
filelist1 = do.call("rbind",strsplit(MACALfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
filelist2$DS = "MACA"
MACALdat = filelist2[,c(3,4,7,6)]
names(MACALdat) = names(EDQMdat)

#MACA METDATA
filelist1 = do.call("rbind",strsplit(MACAMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "METDATA"
filelist2$DS = "MACA"
MACAMdat = filelist2[,c(3,4,8,7)]
names(MACAMdat) = names(EDQMdat)

#BCCA
filelist1 = do.call("rbind",strsplit(BCCAfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "MAURER"
filelist2$DS = "BCCA"
BCCAdat = filelist2[,c(5,7,10,9)]
names(BCCAdat) = names(EDQMdat)

#All metadata
GCM = rep(NA,5)
exp = rep(NA,5)
DS = rep(NA,5)
training = c("DAYMET","LIVNEH","PRISM","METDATA","MAURER")
obsdat = data.frame(GCM,exp,DS,training)

metadat = rbind(EDQMdat,LOCAdat)
metadat = rbind(metadat,MACALdat)
metadat = rbind(metadat,MACAMdat)
metadat = rbind(metadat,BCCAdat)
metadat = rbind(metadat,obsdat)

# all files
allfiles = c(EDQMfiles,LOCAfiles,MACALfiles,MACAMfiles,BCCAfiles,DAYMETfile,LIVNEHfile,PRISMfile,METDATAfile,MAURERfile)

###
# Get data

tasmaxvardatalist = list()
for(i in 1:length(allfiles)){
  nctest = nc_open(allfiles[i])
  tasmaxvardatalist[[i]] = ncvar_get(nctest,nctest$var[[1]]$name)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

###
# Get weighted and unweighted ensemble mean
weightedmeanarray = array(0,dim=c(length(lon),length(lat),12,length(Duvec)))

for(m in 1:12){
  for(c in 1:length(Duvec)){
    for(i in 1:nummodels){
      tmp = tasmaxvardatalist[[i]][,,m]
      weightedmeanarray[,,m,c]=weightedmeanarray[,,m,c]+(tmp*WA[i,c])
    }
  }
  message("Finished weighted ensemble mean calculation for month: ",m)
}

unweightedmeanarray = array(NA,dim=c(length(lon),length(lat),12))
for(m in 1:12){
  tmparray = array(NA,dim=c(length(lon),length(lat),nummodels))
  for(i in 1:nummodels){
    tmparray[,,i] = tasmaxvardatalist[[i]][,,m]
  }
  unweightedmeanarray[,,m]=apply(tmparray,c(1,2),mean,na.rm=TRUE)
  message("Finished unweighted ensemble mean calculation for month: ",m)
}

###
# Get weighted and unweighted RMSEs vs. all 5 obs data

unweightedRMSE = NULL

for(i in 1:length(obsidx)){
  
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  
  obsdat = as.vector(tasmaxvardatalist[[obsidx[i]]])
  moddat = as.vector(unweightedmeanarray)
  RMSEval = RMSE(exp=moddat,obs=obsdat)
  frametmp = data.frame(obs,RMSEval)
  frametmp$obs = as.character(frametmp$obs)
  unweightedRMSE = rbind(unweightedRMSE,frametmp)
  unweightedRMSE$obs = as.character(unweightedRMSE$obs)
}

weightedRMSE = NULL
for(i in 1:length(obsidx)){
  obsdat = as.vector(tasmaxvardatalist[[obsidx[i]]])
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  frametmp2 = NULL
  for(c in 1:length(Duvec)){
    Du = Duvec[c]
    moddat = as.vector(weightedmeanarray[,,,c])
    RMSEval = RMSE(exp=moddat,obs=obsdat)
    frametmp = data.frame(obs,Du,RMSEval)
    frametmp2 = rbind(frametmp2,frametmp)
  }
  frametmp2$RMSErel = frametmp2$RMSEval/unweightedRMSE$RMSEval[i]
  frametmp2$obs = as.character(frametmp2$obs)
  weightedRMSE = rbind(weightedRMSE,frametmp2)
  weightedRMSE$obs = as.character(weightedRMSE$obs)
}

tasmaxweightedRMSE = weightedRMSE

###############
# Switching to tasmin

var = "tasmin"

EDQMfiles = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*EDQM*historical*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",var,"_*histclimo.nc",sep=""),intern=TRUE)
MACALfiles = system(paste("ls /home/woot0002/MACA_LIVNEH/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)
MACAMfiles = system(paste("ls /home/woot0002/MACA_METDATA/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)
BCCAfiles = system(paste("ls /home/woot0002/BCCA/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)

DAYMETfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*daymet*.nc",sep=""),intern=TRUE)
LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*livneh*.nc",sep=""),intern=TRUE)
PRISMfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*prism*.nc",sep=""),intern=TRUE)
METDATAfile = system(paste("ls /data2/3to5/I35/METDATA/",var,"_histclimo_mon.nc",sep=""),intern=TRUE)
MAURERfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_maurer_histclimo_mon.nc",sep=""),intern=TRUE)

###
# create full filelist + metadata table

#EDQM
filelist1 = do.call("rbind",strsplit(EDQMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,5],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$DS = "EDQM"
filelist2$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filelist2$training = rep(c("Daymet","Livneh","PRISM"),3)
filelist2$exp = filelist2[,5]
EDQMdat = filelist2[,c(10,12,9,11)]

#LOCA
filelist1 = do.call("rbind",strsplit(LOCAfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCAdat = filelist2[,c(2,3,4,6)]
names(LOCAdat) = names(EDQMdat)

#MACA Livneh
filelist1 = do.call("rbind",strsplit(MACALfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
filelist2$DS = "MACA"
MACALdat = filelist2[,c(3,4,7,6)]
names(MACALdat) = names(EDQMdat)

#MACA METDATA
filelist1 = do.call("rbind",strsplit(MACAMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "METDATA"
filelist2$DS = "MACA"
MACAMdat = filelist2[,c(3,4,8,7)]
names(MACAMdat) = names(EDQMdat)

#BCCA
filelist1 = do.call("rbind",strsplit(BCCAfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "MAURER"
filelist2$DS = "BCCA"
BCCAdat = filelist2[,c(5,7,10,9)]
names(BCCAdat) = names(EDQMdat)

#All metadata
GCM = rep(NA,5)
exp = rep(NA,5)
DS = rep(NA,5)
training = c("DAYMET","LIVNEH","PRISM","METDATA","MAURER")
obsdat = data.frame(GCM,exp,DS,training)

metadat = rbind(EDQMdat,LOCAdat)
metadat = rbind(metadat,MACALdat)
metadat = rbind(metadat,MACAMdat)
metadat = rbind(metadat,BCCAdat)
metadat = rbind(metadat,obsdat)

# all files
allfiles = c(EDQMfiles,LOCAfiles,MACALfiles,MACAMfiles,BCCAfiles,DAYMETfile,LIVNEHfile,PRISMfile,METDATAfile,MAURERfile)

###
# Get data

tasmaxvardatalist = list()
for(i in 1:length(allfiles)){
  nctest = nc_open(allfiles[i])
  tasmaxvardatalist[[i]] = ncvar_get(nctest,nctest$var[[1]]$name)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

###
# Get weighted and unweighted ensemble mean
weightedmeanarray = array(0,dim=c(length(lon),length(lat),12,length(Duvec)))

for(m in 1:12){
  for(c in 1:length(Duvec)){
    for(i in 1:nummodels){
      tmp = tasmaxvardatalist[[i]][,,m]
      weightedmeanarray[,,m,c]=weightedmeanarray[,,m,c]+(tmp*WA[i,c])
    }
  }
  message("Finished weighted ensemble mean calculation for month: ",m)
}

unweightedmeanarray = array(NA,dim=c(length(lon),length(lat),12))
for(m in 1:12){
  tmparray = array(NA,dim=c(length(lon),length(lat),nummodels))
  for(i in 1:nummodels){
    tmparray[,,i] = tasmaxvardatalist[[i]][,,m]
  }
  unweightedmeanarray[,,m]=apply(tmparray,c(1,2),mean,na.rm=TRUE)
  message("Finished unweighted ensemble mean calculation for month: ",m)
}

###
# Get weighted and unweighted RMSEs vs. all 5 obs data

unweightedRMSE = NULL

for(i in 1:length(obsidx)){
  
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  
  obsdat = as.vector(tasmaxvardatalist[[obsidx[i]]])
  moddat = as.vector(unweightedmeanarray)
  RMSEval = RMSE(exp=moddat,obs=obsdat)
  frametmp = data.frame(obs,RMSEval)
  frametmp$obs = as.character(frametmp$obs)
  unweightedRMSE = rbind(unweightedRMSE,frametmp)
  unweightedRMSE$obs = as.character(unweightedRMSE$obs)
}

weightedRMSE = NULL
for(i in 1:length(obsidx)){
  obsdat = as.vector(tasmaxvardatalist[[obsidx[i]]])
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  frametmp2 = NULL
  for(c in 1:length(Duvec)){
    Du = Duvec[c]
    moddat = as.vector(weightedmeanarray[,,,c])
    RMSEval = RMSE(exp=moddat,obs=obsdat)
    frametmp = data.frame(obs,Du,RMSEval)
    frametmp2 = rbind(frametmp2,frametmp)
  }
  frametmp2$RMSErel = frametmp2$RMSEval/unweightedRMSE$RMSEval[i]
  frametmp2$obs = as.character(frametmp2$obs)
  weightedRMSE = rbind(weightedRMSE,frametmp2)
  weightedRMSE$obs = as.character(weightedRMSE$obs)
}

tasminweightedRMSE = weightedRMSE

############
# for pr

var = "pr"

EDQMfiles = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*QDM*historical*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",var,"_*histclimo.nc",sep=""),intern=TRUE)
MACALfiles = system(paste("ls /home/woot0002/MACA_LIVNEH/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)
MACAMfiles = system(paste("ls /home/woot0002/MACA_METDATA/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)
BCCAfiles = system(paste("ls /home/woot0002/BCCA/regrid/*",var,"_*histclimo.nc",sep=""),intern=TRUE)

DAYMETfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*daymet*.nc",sep=""),intern=TRUE)
LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*livneh*.nc",sep=""),intern=TRUE)
PRISMfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_day*prism*.nc",sep=""),intern=TRUE)
METDATAfile = system(paste("ls /data2/3to5/I35/METDATA/",var,"_histclimo_mon.nc",sep=""),intern=TRUE)
MAURERfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_maurer_histclimo_mon.nc",sep=""),intern=TRUE)

###
# create full filelist + metadata table

#EDQM
filelist1 = do.call("rbind",strsplit(EDQMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,5],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$DS = "EDQM"
filelist2$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filelist2$training = rep(c("Daymet","Livneh","PRISM"),3)
filelist2$exp = filelist2[,5]
EDQMdat = filelist2[,c(10,12,9,11)]

#LOCA
filelist1 = do.call("rbind",strsplit(LOCAfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCAdat = filelist2[,c(2,3,4,6)]
names(LOCAdat) = names(EDQMdat)

#MACA Livneh
filelist1 = do.call("rbind",strsplit(MACALfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
filelist2$DS = "MACA"
MACALdat = filelist2[,c(3,4,7,6)]
names(MACALdat) = names(EDQMdat)

#MACA METDATA
filelist1 = do.call("rbind",strsplit(MACAMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "METDATA"
filelist2$DS = "MACA"
MACAMdat = filelist2[,c(3,4,8,7)]
names(MACAMdat) = names(EDQMdat)

#BCCA
filelist1 = do.call("rbind",strsplit(BCCAfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "MAURER"
filelist2$DS = "BCCA"
BCCAdat = filelist2[,c(5,7,10,9)]
names(BCCAdat) = names(EDQMdat)

#All metadata
GCM = rep(NA,5)
exp = rep(NA,5)
DS = rep(NA,5)
training = c("DAYMET","LIVNEH","PRISM","METDATA","MAURER")
obsdat = data.frame(GCM,exp,DS,training)

metadat = rbind(EDQMdat,LOCAdat)
metadat = rbind(metadat,MACALdat)
metadat = rbind(metadat,MACAMdat)
metadat = rbind(metadat,BCCAdat)
metadat = rbind(metadat,obsdat)

# all files
allfiles = c(EDQMfiles,LOCAfiles,MACALfiles,MACAMfiles,BCCAfiles,DAYMETfile,LIVNEHfile,PRISMfile,METDATAfile,MAURERfile)

###
# Get data

tasmaxvardatalist = list()
for(i in 1:length(allfiles)){
  nctest = nc_open(allfiles[i])
  tasmaxvardatalist[[i]] = ncvar_get(nctest,nctest$var[[1]]$name)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

###
# Get weighted and unweighted ensemble mean
weightedmeanarray = array(0,dim=c(length(lon),length(lat),12,length(Duvec)))

for(m in 1:12){
  for(c in 1:length(Duvec)){
    for(i in 1:nummodels){
      tmp = tasmaxvardatalist[[i]][,,m]
      weightedmeanarray[,,m,c]=weightedmeanarray[,,m,c]+(tmp*WA[i,c])
    }
  }
  message("Finished weighted ensemble mean calculation for month: ",m)
}

unweightedmeanarray = array(NA,dim=c(length(lon),length(lat),12))
for(m in 1:12){
  tmparray = array(NA,dim=c(length(lon),length(lat),nummodels))
  for(i in 1:nummodels){
    tmparray[,,i] = tasmaxvardatalist[[i]][,,m]
  }
  unweightedmeanarray[,,m]=apply(tmparray,c(1,2),mean,na.rm=TRUE)
  message("Finished unweighted ensemble mean calculation for month: ",m)
}

###
# Get weighted and unweighted RMSEs vs. all 5 obs data

unweightedRMSE = NULL

for(i in 1:length(obsidx)){
  
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  
  obsdat = as.vector(tasmaxvardatalist[[obsidx[i]]])
  moddat = as.vector(unweightedmeanarray)
  RMSEval = RMSE(exp=moddat,obs=obsdat)
  frametmp = data.frame(obs,RMSEval)
  frametmp$obs = as.character(frametmp$obs)
  unweightedRMSE = rbind(unweightedRMSE,frametmp)
  unweightedRMSE$obs = as.character(unweightedRMSE$obs)
}

weightedRMSE = NULL
for(i in 1:length(obsidx)){
  obsdat = as.vector(tasmaxvardatalist[[obsidx[i]]])
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  frametmp2 = NULL
  for(c in 1:length(Duvec)){
    Du = Duvec[c]
    moddat = as.vector(weightedmeanarray[,,,c])
    RMSEval = RMSE(exp=moddat,obs=obsdat)
    frametmp = data.frame(obs,Du,RMSEval)
    frametmp2 = rbind(frametmp2,frametmp)
  }
  frametmp2$RMSErel = frametmp2$RMSEval/unweightedRMSE$RMSEval[i]
  frametmp2$obs = as.character(frametmp2$obs)
  weightedRMSE = rbind(weightedRMSE,frametmp2)
  weightedRMSE$obs = as.character(weightedRMSE$obs)
}

prweightedRMSE = weightedRMSE

###########
# plotting

for(i in 1:5){
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  
  cols=c("black","red","blue","purple","green")
  idx = which(tasmaxweightedRMSE$obs==obs)
  tmpframe = tasmaxweightedRMSE[idx,]
  if(i==1){
    plot(RMSErel~Du,data=tmpframe,type="l",lwd=2,col="black",ylim=range(tasmaxweightedRMSE$RMSErel),main="Tasmax RMSE scores",ylab="RMSE score of weighted mean relative to unweighted mean")
  } else {
    lines(RMSErel~Du,data=tmpframe,col=cols[i],lwd=2)
  }
  
}

legend("bottomright",c("DAYMET","LIVNEH","PRISM","METDATA","MAURER"),col=cols,lwd=2)


for(i in 1:5){
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  
  cols=c("black","red","blue","purple","green")
  idx = which(tasminweightedRMSE$obs==obs)
  tmpframe = tasminweightedRMSE[idx,]
  if(i==1){
    plot(RMSErel~Du,data=tmpframe,type="l",lwd=2,col="black",main="Tasmin RMSE scores",ylim=range(tasminweightedRMSE$RMSErel),ylab="RMSE score of weighted mean relative to unweighted mean")
  } else {
    lines(RMSErel~Du,data=tmpframe,col=cols[i],lwd=2)
  }
  
}

legend("topright",c("DAYMET","LIVNEH","PRISM","METDATA","MAURER"),col=cols,lwd=2)


for(i in 1:5){
  if(i==1) obs ="DAYMET"
  if(i==2) obs ="LIVNEH"
  if(i==3) obs ="PRISM"
  if(i==4) obs ="METDATA"
  if(i==5) obs="MAURER"
  
  cols=c("black","red","blue","purple","green")
  idx = which(prweightedRMSE$obs==obs)
  tmpframe = prweightedRMSE[idx,]
  if(i==1){
    plot(RMSErel~Du,data=tmpframe,type="l",lwd=2,col="black",main="Pr RMSE scores",ylim=range(prweightedRMSE$RMSErel),ylab="RMSE score of weighted mean relative to unweighted mean")
  } else {
    lines(RMSErel~Du,data=tmpframe,col=cols[i],lwd=2)
  }
  
}

legend("bottomright",c("DAYMET","LIVNEH","PRISM","METDATA","MAURER"),col=cols,lwd=2)

###

tasmaxRMSEagg = tasminRMSEagg = prRMSEagg = NULL

for(c in 1:length(Duvec)){
  
  Du = Duvec[c]
  
  tmp = subset(tasmaxweightedRMSE,Du==Duvec[c])
  tmaxtmp = apply(tmp[,3:4],2,mean,na.rm=TRUE)
  RMSEval = tmaxtmp[1]
  RMSErel = tmaxtmp[2]
  tmpframe = data.frame(Du,RMSEval,RMSErel)
  tasmaxRMSEagg = rbind(tasmaxRMSEagg,tmpframe)
  
  tmp = subset(tasminweightedRMSE,Du==Duvec[c])
  tmaxtmp = apply(tmp[,3:4],2,mean,na.rm=TRUE)
  RMSEval = tmaxtmp[1]
  RMSErel = tmaxtmp[2]
  tmpframe = data.frame(Du,RMSEval,RMSErel)
  tasminRMSEagg = rbind(tasminRMSEagg,tmpframe)
  
  tmp = subset(prweightedRMSE,Du==Duvec[c])
  tmaxtmp = apply(tmp[,3:4],2,mean,na.rm=TRUE)
  RMSEval = tmaxtmp[1]
  RMSErel = tmaxtmp[2]
  tmpframe = data.frame(Du,RMSEval,RMSErel)
  prRMSEagg = rbind(prRMSEagg,tmpframe)
  
}

###

RMSEagg = rbind(tasmaxRMSEagg,tasminRMSEagg)
RMSEagg = rbind(RMSEagg,prRMSEagg)
RMSEagg2 = aggregate(RMSEagg$RMSErel,by=list(Du=RMSEagg$Du),mean,na.rm=TRUE)

plot(x~Du,data=RMSEagg2,type="l",lwd=2,col="black",main="All RMSE scores by Du radius",ylim=range(c(prRMSEagg$RMSErel,tasmaxRMSEagg$RMSErel,tasminRMSEagg$RMSErel)),ylab="RMSE score of weighted mean relative to unweighted mean")
lines(RMSErel~Du,data=tasmaxRMSEagg,col="red",lwd=1)
lines(RMSErel~Du,data=tasminRMSEagg,col="blue",lwd=1)
lines(RMSErel~Du,data=prRMSEagg,col="green",lwd=1)
abline(v=0.6,lty=3)
legend("bottomright",c("Multivariate Mean","tasmax","tasmin","pr"),col=c("black","red","blue","green"),lwd=c(3,1,1,1))

###

rownames(Wumat) = simnames[1:nummodels]
products = c("EDQM_Daymet","EDQM_Livneh","EDQM_PRISM","LOCA_Livneh","MACA_Livneh","MACA_METDATA","BCCA_MAURER")

avgWumat = matrix(NA,nrow=length(products),ncol=ncol(Wumat))

for(i in 1:length(products)){
  idx = which(paste(metadat$DS,metadat$training,sep="_")==products[i])
  avgWumat[i,]=apply(Wumat[idx,],2,mean,na.rm=TRUE)
}

cols = c("red","blue","green","purple","pink","orange","black")
for(i in 1:length(products)){
  
  if(i==1){
    plot(avgWumat[i,]~Duvec,col=cols[i],lwd=2,type="l",xlim=c(0,1),ylab="Model Independence Weight",xlab="Du")
  } else {
    lines(avgWumat[i,]~Duvec,col=cols[i],lwd=2,type="l")
  }
  
}
legend("topright",legend=products,col=cols,lwd=2)
abline(v=0.4,lty=3)
abline(h=1/length(products),lty=3)
