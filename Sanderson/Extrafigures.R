################
#
# Extra analysis for Downscaling independence analysis.
#
#

library(ncdf4)
library(maps)
library(fields)
library(sp)
library(mailR)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

var = varin = "pr"

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

GCMfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*histclimo.nc",sep=""),intern=TRUE)

GCMprojfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)
LOCAprojfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*projclimo.nc",sep=""),intern=TRUE)

LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)

load("/home/woot0002/DS_ind/GCMlist.Rdata")

GCM_hfiles = GCM_pfiles = LOCA_hfiles = LOCA_pfiles = c()

for(i in 1:length(GCMlist)){
  
  GCM_hfiles[i] = GCMfiles[grep(paste(GCMlist[i],"_",sep=""),GCMfiles)]
  GCM_pfiles[i] = GCMprojfiles[grep(paste(GCMlist[i],"_",sep=""),GCMprojfiles)]
  
  LOCA_hfiles[i] = LOCAfiles[grep(paste(GCMlist[i],"_",sep=""),LOCAfiles)]
  LOCA_pfiles[i] = LOCAprojfiles[grep(paste(GCMlist[i],"_",sep=""),LOCAprojfiles)]
  
}

###
# create full filelist + metadata table - historical

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_hfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMhdat = filelist2[,c(2,3,4,6)]
names(GCMhdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_hfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCAhdat = filelist2[,c(2,3,4,6)]
names(LOCAhdat) = names(GCMhdat)

#All metadata
GCM = rep(NA,1)
exp = rep(NA,1)
DS = rep(NA,1)
training = "LIVNEH"
obsdat = data.frame(GCM,exp,DS,training)

GCMhdat = rbind(GCMhdat,obsdat)
LOCAhdat= rbind(LOCAhdat,obsdat)

# all files
GCMgroup = c(GCM_hfiles,LIVNEHfile)
LOCAgroup = c(LOCA_hfiles,LIVNEHfile)

###
# create full filelist + metadata table - projected

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_pfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMpdat = filelist2[,c(2,3,4,6)]
names(GCMpdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_pfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCApdat = filelist2[,c(2,3,4,6)]
names(LOCApdat) = names(GCMpdat)

# all files
GCMpgroup = GCM_pfiles
LOCApgroup = LOCA_pfiles

#########
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

######
# Gather data

ncvarname = paste(var,"climo",sep="")

### GCM hist + Livneh
GCMhvardatalist = list()
for(i in 1:length(GCMgroup)){
  nctest = nc_open(GCMgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  GCMhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(GCMhvardatalist,mean,na.rm=TRUE)

### GCM projected change
GCMpvardatalist = list()
for(i in 1:length(GCMpgroup)){
  nctest = nc_open(GCMpgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  GCMpvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)-GCMhvardatalist[[i]]
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(GCMpvardatalist,mean,na.rm=TRUE)


### LOCA historical + Livneh

LOCAhvardatalist = list()
for(i in 1:length(LOCAgroup)){
  nctest = nc_open(LOCAgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  LOCAhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCAhvardatalist,mean,na.rm=TRUE)

### LOCA projected change

LOCApvardatalist = list()
for(i in 1:length(LOCApgroup)){
  nctest = nc_open(LOCApgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
  LOCApvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE) - LOCAhvardatalist[[i]]
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCApvardatalist,mean,na.rm=TRUE)

######
# plotting

pdf("/home/woot0002/DS_ind/Sanderson_hist_pr.pdf",onefile=TRUE,height=5,width=15)
for(i in 1:length(GCMlist)){
  
  GCMdat = GCMhvardatalist[[i]]
  LOCAdat = LOCAhvardatalist[[i]]
  Livnehdat = GCMhvardatalist[[27]]
  
  GCMdat = ifelse(is.na(LOCAdat)==FALSE,GCMdat,NA)
  Livnehdat = ifelse(is.na(LOCAdat)==FALSE,Livnehdat,NA)
  
  datavec = c(GCMdat,LOCAdat,Livnehdat)
  datacolorbar=colorramp(datavec,colorchoice="whitetogreen",Blimit=30,use_fixed_scale = FALSE,fixed_scale=c(0,1400))
  
  par(mfrow=c(1,3))
  
  testsfc1 = list(x=lon-360,y=lat,z=GCMdat)
  surface(testsfc1,type="I",main=paste("Raw: ",GCMlist[i],sep=""),zlim=datacolorbar[[1]],col=datacolorbar[[3]],breaks=datacolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(GCMdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(GCMdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(GCMdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
  testsfc1 = list(x=lon-360,y=lat,z=LOCAdat)
  surface(testsfc1,type="I",main=paste("LOCA: ",GCMlist[i],sep=""),zlim=datacolorbar[[1]],col=datacolorbar[[3]],breaks=datacolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(LOCAdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(LOCAdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(LOCAdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
  testsfc1 = list(x=lon-360,y=lat,z=Livnehdat)
  surface(testsfc1,type="I",main="Livneh v. 1.2",zlim=datacolorbar[[1]],col=datacolorbar[[3]],breaks=datacolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(Livnehdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(Livnehdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(Livnehdat,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
}

dev.off()

########
# plot projected change, most / least independent

load("/home/woot0002/DS_ind/Sanderson_EnsembleWeights_v4.Rdata")

idxorder = order(GCMhdat$Wuc,decreasing=TRUE)

GCMhdat2 = GCMhdat[idxorder,]
LOCAhdat2 = LOCAhdat[idxorder,]

GCMmost = GCMpvardatalist[[idxorder[1]]]
LOCAmost = LOCApvardatalist[[idxorder[1]]]

GCMleast = GCMpvardatalist[[idxorder[length(idxorder)]]]
LOCAleast = LOCApvardatalist[[idxorder[length(idxorder)]]]

GCMmost = ifelse(is.na(LOCAmost)==FALSE,GCMmost,NA)
GCMleast = ifelse(is.na(LOCAleast)==FALSE,GCMleast,NA)

datavec = c(GCMmost,LOCAmost,GCMleast,LOCAleast)
datacolorbar=colorramp(datavec,colorchoice="browntogreen",Blimit=30,use_fixed_scale = FALSE,fixed_scale=c(0,1400))

pdf("/home/woot0002/DS_ind/Sanderson_change_pr.pdf",onefile=TRUE,height=10,width=10)

par(mfrow=c(2,2))

testsfc1 = list(x=lon-360,y=lat,z=GCMmost)
surface(testsfc1,type="I",main=paste("Raw: ",GCMlist[idxorder[1]],sep=""),zlim=datacolorbar[[1]],col=datacolorbar[[3]],breaks=datacolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
text(-108.85,28.45,labels=paste("MAX = ",round(max(GCMmost,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,27.45,labels=paste("MEAN = ",round(mean(GCMmost,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,26.45,labels=paste("MIN = ",round(min(GCMmost,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lon-360,y=lat,z=GCMleast)
surface(testsfc1,type="I",main=paste("Raw: ",GCMlist[idxorder[length(idxorder)]],sep=""),zlim=datacolorbar[[1]],col=datacolorbar[[3]],breaks=datacolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
text(-108.85,28.45,labels=paste("MAX = ",round(max(GCMleast,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,27.45,labels=paste("MEAN = ",round(mean(GCMleast,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,26.45,labels=paste("MIN = ",round(min(GCMleast,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lon-360,y=lat,z=LOCAmost)
surface(testsfc1,type="I",main=paste("LOCA: ",GCMlist[idxorder[1]],sep=""),zlim=datacolorbar[[1]],col=datacolorbar[[3]],breaks=datacolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
text(-108.85,28.45,labels=paste("MAX = ",round(max(LOCAmost,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,27.45,labels=paste("MEAN = ",round(mean(LOCAmost,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,26.45,labels=paste("MIN = ",round(min(LOCAmost,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lon-360,y=lat,z=LOCAleast)
surface(testsfc1,type="I",main=paste("LOCA: ",GCMlist[idxorder[length(idxorder)]],sep=""),zlim=datacolorbar[[1]],col=datacolorbar[[3]],breaks=datacolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
text(-108.85,28.45,labels=paste("MAX = ",round(max(LOCAleast,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,27.45,labels=paste("MEAN = ",round(mean(LOCAleast,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-108.85,26.45,labels=paste("MIN = ",round(min(LOCAleast,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

dev.off()

##############
# means calcs

GCMhmat = sapply(GCMhvardatalist,mean,na.rm=TRUE)
LOCAhmat =sapply(LOCAhvardatalist,mean,na.rm=TRUE)
GCMcmat =sapply(GCMpvardatalist,mean,na.rm=TRUE)
LOCAcmat = sapply(LOCApvardatalist,mean,na.rm=TRUE)

for(i in 1:length(GCMlist)){
  for(m in 1:12){
    GCMhmat[i,m] = mean(GCMhvardatalist[[i]][,,m],na.rm=TRUE)
    LOCAhmat[i,m] = mean(LOCAhvardatalist[[i]][,,m],na.rm=TRUE)
    GCMcmat[i,m] = mean(GCMpvardatalist[[i]][,,m],na.rm=TRUE)
    LOCAcmat[i,m] = mean(LOCApvardatalist[[i]][,,m],na.rm=TRUE)
  }
}










