#####################
#
# RMSE calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = varin = "pr50"

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

GCMfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*histclimo.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*histclimo.nc",sep=""),intern=TRUE)

LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)

###
# create full filelist + metadata table

#GCMs
filelist1 = do.call("rbind",strsplit(GCMfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
GCMdat = filelist2[,c(2,3,4,6)]
names(GCMdat) = c("GCM","exp","DS","training")


#LOCA
filelist1 = do.call("rbind",strsplit(LOCAfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCAdat = filelist2[,c(2,3,4,6)]
names(LOCAdat) = names(GCMdat)

LOCAfiles = LOCAfiles[which(LOCAdat$GCM %in% GCMdat$GCM)]
LOCAdat = LOCAdat[which(LOCAdat$GCM %in% GCMdat$GCM),]

#All metadata
GCM = rep(NA,1)
exp = rep(NA,1)
DS = rep(NA,1)
training = "LIVNEH"
obsdat = data.frame(GCM,exp,DS,training)

GCMdat = rbind(GCMdat,obsdat)
LOCAdat= rbind(LOCAdat,obsdat)

# all files
GCMgroup = c(GCMfiles,LIVNEHfile)
LOCAgroup = c(LOCAfiles,LIVNEHfile)

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

GCMvardatalist = list()
for(i in 1:length(GCMgroup)){
  nctest = nc_open(GCMgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  GCMvardatalist[[i]] = ncvar_get(nctest,nctest$var[[idx]]$name)
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(GCMvardatalist,mean,na.rm=TRUE)

## 

LOCAvardatalist = list()
for(i in 1:length(LOCAgroup)){
  nctest = nc_open(LOCAgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==var)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  LOCAvardatalist[[i]] = ncvar_get(nctest,nctest$var[[idx]]$name)
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCAvardatalist,mean,na.rm=TRUE)

####
# Apply the area weighting
GCMvardatalatweight = LOCAvardatalatweight = list()
latweights = cos((lat*pi)/180)
latweightmat = matrix(rep(latweights,each=length(lon)),nrow=length(lon),ncol=length(latweights))
latweightvector = as.vector(latweightmat)

for(i in 1:length(GCMgroup)){
  GCMvardatalatweight[[i]]=weightfunc(GCMvardatalist[[i]],latweightmat = latweightmat)
  LOCAvardatalatweight[[i]]=weightfunc(LOCAvardatalist[[i]],latweightmat = latweightmat)
}

####
# RMSE matrix calculation

GCMRMSEmat = LOCARMSEmat = matrix(NA,nrow=length(GCMgroup),ncol=(length(GCMgroup)))

for(i in 1:length(GCMgroup)){
  for(j in 1:length(GCMgroup)){
    if(i!=j){
      GCMRMSEmat[i,j] = RMSE(exp=as.vector(GCMvardatalatweight[[i]]),obs=as.vector(GCMvardatalatweight[[j]]))
      LOCARMSEmat[i,j] = RMSE(exp=as.vector(LOCAvardatalatweight[[i]]),obs=as.vector(LOCAvardatalatweight[[j]]))
    }
  }
  message("RMSE calculated for ",i," / ",length(GCMgroup))
}

####
# normalize RMSE matrix

normGCMRMSEmat = GCMRMSEmat/mean(GCMvardatalatweight[[27]],na.rm=TRUE)
normLOCARMSEmat = LOCARMSEmat/mean(LOCAvardatalatweight[[27]],na.rm=TRUE)

####
# save output in Rdata files

save(list=c("GCMRMSEmat","normGCMRMSEmat","GCMdat","GCMgroup"),file=paste("/home/woot0002/RMSEfiles/",var,"_RMSEmats_GCM.Rdata",sep=""))
save(list=c("LOCARMSEmat","normLOCARMSEmat","LOCAdat","LOCAgroup"),file=paste("/home/woot0002/RMSEfiles/",var,"_RMSEmats_LOCA.Rdata",sep=""))

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

simnamesGCM = simnamesLOCA = c()
for(i in 1:nrow(GCMdat)){
if(i<=(nrow(GCMdat)-1)){
  simnamesGCM[i]= paste(GCMdat[i,1],GCMdat[i,2],GCMdat[i,3],GCMdat[i,4],sep="_")
  simnamesLOCA[i]= paste(LOCAdat[i,1],LOCAdat[i,2],LOCAdat[i,3],LOCAdat[i,4],sep="_")
} else {
  simnamesGCM[i]= GCMdat[i,4]
  simnamesLOCA[i]= LOCAdat[i,4]
}
}

rownames(normGCMRMSEmat) = colnames(normGCMRMSEmat) = simnamesGCM
rownames(normLOCARMSEmat) = colnames(normLOCARMSEmat) = simnamesLOCA

del.RMSEGCM = get_upper_tri(normGCMRMSEmat)
del.RMSELOCA = get_upper_tri(normLOCARMSEmat)

#del.earlymat = ifelse(del.earlymat==0,NA,del.earlymat)

melted_cormatGCM <- melt(del.RMSEGCM, na.rm = TRUE)
melted_cormatLOCA <- melt(del.RMSELOCA, na.rm = TRUE)

topend = ceiling(max(c(normGCMRMSEmat,normLOCARMSEmat),na.rm=TRUE))

### Heatmap
pdf(paste(var,"_NRMSEheatmap_GCM.pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormatGCM, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = 1, limit = c(0,topend), breaks=seq(0,topend,by=0.25), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

pdf(paste(var,"_NRMSEheatmap_LOCA.pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormatLOCA, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = 1, limit = c(0,topend), breaks=seq(0,topend,by=0.25), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()


