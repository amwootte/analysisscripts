#####################
#
# RMSE calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = varin = "rx1day"

if(var=="tmin"){
  varin="tasmin"
}
if(var=="tmax"){
  varin="tasmax"
}
if(var=="tmax95"){
  varin="tasmax"
}
if(var=="tmin32"){
  varin="tasmin"
}
if(var=="r1mm"){
  varin="pr"
}
if(var=="pr50"){
  varin="pr"
}
if(var=="rx1day"){
  varin="pr"
}
if(var=="rx5day"){
  varin="pr"
}

if(varin=="pr"){
  EDQMfiles = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*QDM*historical*.nc",sep=""),intern=TRUE)
} else {
  EDQMfiles = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*EDQM*historical*.nc",sep=""),intern=TRUE)
}
PARMfiles = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*PARM*historical*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*histclimo.nc",sep=""),intern=TRUE)
MACALfiles = system(paste("ls /home/woot0002/MACA_LIVNEH/regrid/*",varin,"_*histclimo.nc",sep=""),intern=TRUE)
MACAMfiles = system(paste("ls /home/woot0002/MACA_METDATA/regrid/*",varin,"_*histclimo.nc",sep=""),intern=TRUE)
BCCAfiles = system(paste("ls /home/woot0002/BCCA/regrid/*",varin,"_*histclimo.nc",sep=""),intern=TRUE)

DAYMETfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*daymet*.nc",sep=""),intern=TRUE)
LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)
PRISMfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*prism*.nc",sep=""),intern=TRUE)
if(var!="tmax" & var!="tmin"){
  METDATAfile = system(paste("ls /data2/3to5/I35/METDATA/",var,"_histclimo_mon.nc",sep=""),intern=TRUE)
  MAURERfile = system(paste("ls /home/woot0002/monthlyclimo/",var,"_maurer_histclimo_mon.nc",sep=""),intern=TRUE)
} else {
  METDATAfile = system(paste("ls /data2/3to5/I35/METDATA/",varin,"_histclimo_mon.nc",sep=""),intern=TRUE)
  MAURERfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_maurer_histclimo_mon.nc",sep=""),intern=TRUE)
}
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

#PARM
filelist1 = do.call("rbind",strsplit(PARMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,5],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$DS = "PARM"
filelist2$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filelist2$training = rep(c("Daymet","Livneh","PRISM"),3)
filelist2$exp = filelist2[,5]
PARMdat = filelist2[,c(10,12,9,11)]

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

metadat = rbind(EDQMdat,PARMdat)
metadat = rbind(metadat,LOCAdat)
metadat = rbind(metadat,MACALdat)
metadat = rbind(metadat,MACAMdat)
metadat = rbind(metadat,BCCAdat)
metadat = rbind(metadat,obsdat)

# all files
allfiles = c(EDQMfiles,PARMfiles,LOCAfiles,MACALfiles,MACAMfiles,BCCAfiles,DAYMETfile,LIVNEHfile,PRISMfile,METDATAfile,MAURERfile)

####
# Extra Functions
weightfunc = function(datafile,latweightmat){
  results=array(data=NA,dim=c(nrow(datafile),ncol(datafile),12))
  for(i in 1:12){
    results[,,i]=datafile[,,i]*latweightmat
  }
  return(results)
}

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

####
# Gather data

ncvarname = paste(var,"climo",sep="")

vardatalist = list()
for(i in 1:length(allfiles)){
  nctest = nc_open(allfiles[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  vardatalist[[i]] = ncvar_get(nctest,nctest$var[[idx]]$name)
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(vardatalist,mean,na.rm=TRUE)

####
# Apply the area weighting
vardatalatweight = list()
latweights = cos((lat*pi)/180)
latweightmat = matrix(rep(latweights,each=length(lon)),nrow=length(lon),ncol=length(latweights))
latweightvector = as.vector(latweightmat)

for(i in 1:length(allfiles)){
  vardatalatweight[[i]]=weightfunc(vardatalist[[i]],latweightmat = latweightmat)
}
sapply(vardatalatweight,mean,na.rm=TRUE)

####
# RMSE matrix calculation

RMSEmat = matrix(NA,nrow=length(allfiles),ncol=(length(allfiles)))

for(i in 1:length(allfiles)){
  for(j in 1:length(allfiles)){
    if(i!=j){
      RMSEmat[i,j] = RMSE(exp=as.vector(vardatalatweight[[i]]),obs=as.vector(vardatalatweight[[j]]))
    }
  }
  message("RMSE calculated for ",i," / ",length(allfiles))
}

####
# normalize RMSE matrix

obsdatvec = c()
for(i in 131:135){
  obsdatvec = c(obsdatvec,vardatalatweight[[i]])
}

#normRMSEmat = RMSEmat/sd(obsdatvec,na.rm=TRUE) # version 1
normRMSEmat = RMSEmat/mean(RMSEmat,na.rm=TRUE) # version 2

####
# save output in Rdata files

save(list=c("RMSEmat","normRMSEmat","metadat","allfiles"),file=paste("/home/woot0002/RMSEfiles/",var,"_RMSEmats_v2.Rdata",sep=""))

####
# heatmaps

library(ggplot2)
library(reshape2)
get_upper_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

simnames = c()
for(i in 1:nrow(metadat)){
if(i<=(nrow(metadat)-5)){
  simnames[i]= paste(metadat[i,1],metadat[i,2],metadat[i,3],metadat[i,4],sep="_")
} else {
  simnames[i]= metadat[i,4]
}
}

rownames(normRMSEmat) = colnames(normRMSEmat) = simnames


del.RMSE = get_upper_tri(normRMSEmat)

#del.earlymat = ifelse(del.earlymat==0,NA,del.earlymat)

melted_cormat1 <- melt(del.RMSE, na.rm = TRUE)

topend = ceiling(max(normRMSEmat,na.rm=TRUE))

### Heatmap
pdf(paste(var,"_NRMSEheatmap_v2.pdf",sep=""),height=15,width=15)
ggplot(data = melted_cormat1, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.25), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 5, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 5, hjust = 1))+
  coord_fixed()
dev.off()




