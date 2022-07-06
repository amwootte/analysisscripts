#####################
#
# RMSE calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = varin = "pr"
type="ann"
statemask = "new mexico"

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
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)

GCMprojfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)
LOCAprojfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)

LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)

load("/home/woot0002/DS_ind/manuscript1/GCMlist.Rdata")

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

#####
# get domain mask

if(statemask!="full"){
  test = nc_open(paste("/home/woot0002/DS_ind/",statemask,"_mask.nc",sep=""))
  regionmask = ncvar_get(test,"mask")
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
}

#########
####
# Extra Functions
weightfunc = function(datafile,latweightmat){
  
  if(is.na(dim(datafile)[3])==FALSE){
    Tsize = dim(datafile)[3]
  } else {
    Tsize=1
  }
  
  results=array(data=NA,dim=c(nrow(datafile),ncol(datafile),Tsize))
  for(i in 1:Tsize){
    if(Tsize>1){
      results[,,i]=datafile[,,i]*latweightmat
    } else {
      results[,,i]=datafile*latweightmat
    }
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
  if(type!="ann"){
    GCMhvardatalist[[i]] = tmp
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      GCMhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(statemask=="full"){
        GCMhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist[[i]],NA)
      } else {
        GCMhvardatalist[[i]] = ifelse(regionmask==1,GCMhvardatalist[[i]],NA)
      }
    } 
    if(var=="tmax" | var=="tmin"){
      GCMhvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(statemask=="full"){
        GCMhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist[[i]],NA)
      } else {
        GCMhvardatalist[[i]] = ifelse(regionmask==1,GCMhvardatalist[[i]],NA)
      }
    } 
  }
  
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
  if(type!="ann"){
    GCMpvardatalist[[i]] = tmp-GCMhvardatalist[[i]]
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      GCMpvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(statemask=="full"){
        GCMpvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist[[i]],NA)-GCMhvardatalist[[i]]
      } else {
        GCMpvardatalist[[i]] = ifelse(regionmask==1,GCMpvardatalist[[i]],NA)-GCMhvardatalist[[i]]
      }
    } 
    if(var=="tmax" | var=="tmin"){
      GCMpvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(statemask=="full"){
        GCMpvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist[[i]],NA)-GCMhvardatalist[[i]]
      } else {
        GCMpvardatalist[[i]] = ifelse(regionmask==1,GCMpvardatalist[[i]],NA)-GCMhvardatalist[[i]]
      }
    } 
  }
  
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
  if(type!="ann"){
    GCMhvardatalist[[i]] = tmp
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      LOCAhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(statemask=="full"){
        LOCAhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist[[i]],NA)
      } else {
        LOCAhvardatalist[[i]] = ifelse(regionmask==1,LOCAhvardatalist[[i]],NA)
      }
    } 
    if(var=="tmax" | var=="tmin"){
      LOCAhvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(statemask=="full"){
        LOCAhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist[[i]],NA)
      } else {
        LOCAhvardatalist[[i]] = ifelse(regionmask==1,LOCAhvardatalist[[i]],NA)
      }
    } 
  }

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
  if(type!="ann"){
    LOCApvardatalist[[i]] = tmp-LOCAhvardatalist[[i]]
  } else {
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      LOCApvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(statemask=="full"){
        LOCApvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist[[i]],NA)-LOCAhvardatalist[[i]]
      } else {
        LOCApvardatalist[[i]] = ifelse(regionmask==1,LOCApvardatalist[[i]],NA)-LOCAhvardatalist[[i]]
      }
    } 
    if(var=="tmax" | var=="tmin"){
      LOCApvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(statemask=="full"){
        LOCApvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist[[i]],NA)-LOCAhvardatalist[[i]]
      } else {
        LOCApvardatalist[[i]] = ifelse(regionmask==1,LOCApvardatalist[[i]],NA)-LOCAhvardatalist[[i]]
      }
    } 
  }
  
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCApvardatalist,mean,na.rm=TRUE)


####
# Apply the area weighting
GCMhvardatalatweight = LOCAhvardatalatweight = list()
latweights = cos((lat*pi)/180)
latweightmat = matrix(rep(latweights,each=length(lon)),nrow=length(lon),ncol=length(latweights))
latweightvector = as.vector(latweightmat)

for(i in 1:length(GCMgroup)){
  GCMhvardatalatweight[[i]]=weightfunc(GCMhvardatalist[[i]],latweightmat = latweightmat)
  LOCAhvardatalatweight[[i]]=weightfunc(LOCAhvardatalist[[i]],latweightmat = latweightmat)
}

###

GCMpvardatalatweight = LOCApvardatalatweight = list()
latweights = cos((lat*pi)/180)
latweightmat = matrix(rep(latweights,each=length(lon)),nrow=length(lon),ncol=length(latweights))
latweightvector = as.vector(latweightmat)

for(i in 1:length(GCMpgroup)){
  GCMpvardatalatweight[[i]]=weightfunc(GCMpvardatalist[[i]],latweightmat = latweightmat)
  LOCApvardatalatweight[[i]]=weightfunc(LOCApvardatalist[[i]],latweightmat = latweightmat)
}


####
# RMSE matrix calculation

GCMhRMSEmat = LOCAhRMSEmat = matrix(NA,nrow=length(GCMgroup),ncol=(length(GCMgroup)))

for(i in 1:length(GCMgroup)){
  for(j in 1:length(GCMgroup)){
    if(i!=j){
      GCMhRMSEmat[i,j] = RMSE(exp=as.vector(GCMhvardatalatweight[[i]]),obs=as.vector(GCMhvardatalatweight[[j]]))
      LOCAhRMSEmat[i,j] = RMSE(exp=as.vector(LOCAhvardatalatweight[[i]]),obs=as.vector(LOCAhvardatalatweight[[j]]))
    }
  }
  message("RMSE calculated for ",i," / ",length(GCMgroup))
}

###

GCMpRMSEmat = LOCApRMSEmat = matrix(NA,nrow=length(GCMpgroup),ncol=(length(GCMpgroup)))

for(i in 1:length(GCMpgroup)){
  for(j in 1:length(GCMpgroup)){
    if(i!=j){
      GCMpRMSEmat[i,j] = RMSE(exp=as.vector(GCMpvardatalatweight[[i]]),obs=as.vector(GCMpvardatalatweight[[j]]))
      LOCApRMSEmat[i,j] = RMSE(exp=as.vector(LOCApvardatalatweight[[i]]),obs=as.vector(LOCApvardatalatweight[[j]]))
    }
  }
  message("RMSE calculated for ",i," / ",length(GCMpgroup))
}

####
# normalize hist RMSE matrix

normGCMhRMSEmat = GCMhRMSEmat/sd(GCMhvardatalatweight[[27]],na.rm=TRUE)
normLOCAhRMSEmat = LOCAhRMSEmat/sd(LOCAhvardatalatweight[[27]],na.rm=TRUE)

####
# normalize RMSE change matrix

normGCMpRMSEmat = GCMpRMSEmat/abs(mean(sapply(GCMpvardatalist,sd,na.rm=TRUE)))
normLOCApRMSEmat = LOCApRMSEmat/abs(mean(sapply(LOCApvardatalist,sd,na.rm=TRUE)))

normGCMcRMSEmat = normGCMhRMSEmat
normGCMcRMSEmat[1:26,1:26] = normGCMpRMSEmat

normLOCAcRMSEmat = normLOCAhRMSEmat
normLOCAcRMSEmat[1:26,1:26] = normLOCApRMSEmat

####
# save output in Rdata files

save(list=c("normGCMhRMSEmat","GCMhdat","GCMgroup","normGCMcRMSEmat","GCMpdat","GCMpgroup"),file=paste("/home/woot0002/RMSEfiles/",var,"_RMSEmats_combo_GCM_v4_",type,"_",statemask,".Rdata",sep=""))
save(list=c("normLOCAhRMSEmat","LOCAhdat","LOCAgroup","normLOCAcRMSEmat","LOCApdat","LOCApgroup"),file=paste("/home/woot0002/RMSEfiles/",var,"_RMSEmats_combo_LOCA_v4_",type,"_",statemask,".Rdata",sep=""))

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
for(i in 1:nrow(GCMhdat)){
if(i<=(nrow(GCMhdat)-1)){
  simnamesGCM[i]= paste(GCMhdat[i,1],GCMhdat[i,2],GCMhdat[i,3],GCMhdat[i,4],sep="_")
  simnamesLOCA[i]= paste(LOCAhdat[i,1],LOCAhdat[i,2],LOCAhdat[i,3],LOCAhdat[i,4],sep="_")
} else {
  simnamesGCM[i]= as.character(GCMhdat[i,4])
  simnamesLOCA[i]= as.character(LOCAhdat[i,4])
}
}

rownames(normGCMhRMSEmat) = colnames(normGCMhRMSEmat) = simnamesGCM
rownames(normLOCAhRMSEmat) = colnames(normLOCAhRMSEmat) = simnamesLOCA

del.RMSEGCM = get_upper_tri(normGCMhRMSEmat)
del.RMSELOCA = get_upper_tri(normLOCAhRMSEmat)

#del.earlymat = ifelse(del.earlymat==0,NA,del.earlymat)

melted_cormatGCM <- melt(del.RMSEGCM, na.rm = TRUE)
melted_cormatLOCA <- melt(del.RMSELOCA, na.rm = TRUE)

topend = ceiling(max(c(normGCMhRMSEmat,normLOCAhRMSEmat),na.rm=TRUE))

### Heatmap
pdf(paste("/home/woot0002/DS_ind/",var,"_NRMSEheatmap_GCM_hist_v4_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormatGCM, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

pdf(paste("/home/woot0002/DS_ind/",var,"_NRMSEheatmap_LOCA_hist_v4_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormatLOCA, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

###
# projected change heat maps

rownames(normGCMcRMSEmat) = colnames(normGCMcRMSEmat) = simnamesGCM
rownames(normLOCAcRMSEmat) = colnames(normLOCAcRMSEmat) = simnamesLOCA

del.RMSEGCM = get_upper_tri(normGCMcRMSEmat)
del.RMSELOCA = get_upper_tri(normLOCAcRMSEmat)

#del.earlymat = ifelse(del.earlymat==0,NA,del.earlymat)

melted_cormatGCM <- melt(del.RMSEGCM, na.rm = TRUE)
melted_cormatLOCA <- melt(del.RMSELOCA, na.rm = TRUE)

topend = ceiling(max(c(normGCMcRMSEmat,normLOCAcRMSEmat),na.rm=TRUE))

### Heatmap
pdf(paste("/home/woot0002/DS_ind/",var,"_NRMSEheatmap_GCM_change_v4_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormatGCM, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

pdf(paste("/home/woot0002/DS_ind/",var,"_NRMSEheatmap_LOCA_change_v4_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormatLOCA, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()




