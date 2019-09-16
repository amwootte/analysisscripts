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

Dqvec = Dcvec = seq(0.05,2,by=0.01)
Du = 0.48

##
# similarity weight calculation
S_top = -(normRMSEmat/Du)^2
S = exp(S_top)

Ru = c()
for(i in 1:nummodels){
  Ru[i]=1+sum(S[1:nummodels,i],na.rm=TRUE)
}

Wu = Ru^(-1)

##
# quality weight calculation
Wqmat = matrix(NA,nrow=nummodels,ncol=length(Dqvec))
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
  
  Q_top = -(normRMSEmat[useidx,i]/Dqvec)^2
  Wqmat[i,] = exp(Q_top)
}

##
# consitency weight calculation
Wcmat = matrix(NA,nrow=nummodels,ncol=length(Dcvec))
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
  
  C_top = -(mean(normRMSEmat[useidx,i],na.rm=TRUE)/Dcvec)^2
  Wcmat[i,] = exp(C_top)
}

########
# Calculate final weight

W = WA = array(NA,dim=c(nummodels,length(Dqvec),length(Dcvec)))

for(i in 1:length(Dqvec)){
  for(j in 1:length(Dcvec)){
  W[,i,j] = Wqmat[,i]*Wcmat[,j]*Wu
  WA[,i,j] = W[,i,j]/sum(W[,i,j])
  }
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
weightedmeanarray = array(0,dim=c(length(lon),length(lat),12,length(Dqvec),length(Dcvec)))

for(i in 1:nummodels){
  for(m in 1:12){
    
    tmp = tasmaxvardatalist[[i]][,,m]
    
    for(q in 1:length(Dqvec)){
      for(c in 1:length(Dcvec)){
        
        
        
        
      }
    }
  }
}






##
# set radii
Du = 0.48 # original 0.48 - ???? does this need to be 0.48 * distance between best performing model and obs.
Dq = 0.8
Dc = 0.8



W2 = Wc*Wu
WA2 = W2/sum(W2)

#########
# 

plot(c_dist~t_dist,pch=19,xlim=range(normRMSEmat,na.rm=TRUE),ylim=range(normRMSEmat,na.rm=TRUE))
abline(coef=c(0,1),lty=2)

plot(Wc~Wq,pch=19,xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1),lty=2)

plot(Wc~Wu,pch=19,xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1),lty=2)

simmetdat = metadat[1:nummodels,]
simmetdat$product = paste(simmetdat$DS,simmetdat$training,sep="_")
simmetdat$Wc = Wc
simmetdat$Wu = Wu
simmetdat$Wq = Wq

simmetdat$col = NA

products=unique(simmetdat$product)
cols = c("red","blue","green","purple","pink","orange","black")

for(p in 1:length(products)){
  simmetdat$col[which(simmetdat$product==products[p])] = cols[p]
  if(p==1){
    plot(Wc~Wq,data=simmetdat[which(simmetdat$product==products[p]),],xlim=c(0,1),ylim=c(0,1),col=cols[p])
  } else {
    points(Wc~Wq,data=simmetdat[which(simmetdat$product==products[p]),],col=cols[p]) 
  }
}
abline(coef=c(0,1),lty=2)
legend("topright",legend=products,col=cols,pch=1)

library(ggplot2) 

#Here's the 1st plot 
dg<-qplot(Wq,Wc,colour=Wu,shape=product,data=simmetdat) 
dg + scale_colour_gradient2(low="red", high="blue",mid = "gray96", 
                            midpoint = 0.5,limits=c(0,1),breaks=seq(0,1,by=0.1),guide="legend") +
  xlim(0,1)+ylim(0,1)


meandat = aggregate(simmetdat[,6:8],by=list(product=simmetdat$product),mean,na.rm=TRUE)
