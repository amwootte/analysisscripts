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

########
# Focusing RMSE and Bias dot plots - one variable at a time

varin = "pr"
var = "pr"

weightingfiles = system(paste("ls /home/woot0002/DS_ind/WeightedMeansVars_",var,"_*_WU*_SA*wBMA.Rdata",sep=""),intern=TRUE)
BMApostfiles = system(paste("ls /home/woot0002/DS_ind/BMAPost_WeightedMeansVars_",var,"_*_WU*_SA*.Rdata",sep=""),intern=TRUE)

########
# grid basics for later

GCMfiles_pr = system("ls /home/woot0002/GCMs/regrid/pr_*histclimo*.nc",intern=TRUE)
nctest = nc_open(GCMfiles_pr[1])
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")
nc_close(nctest)

########
# all weighting scheme means plus best BMA

setwd("/home/woot0002/DS_ind/")

meansdat_all = NULL 

for(i in 1:length(weightingfiles)){
  load(weightingfiles[i])
  if(var=="tmax"){
    meansdat_pr = meansdat_tmax
  }
  filelist = strsplit(weightingfiles[i],"_",fixed=TRUE)
  WUdomain = substr(filelist[[1]][5],3,nchar(filelist[[1]][5]))
  SAdomain = substr(filelist[[1]][6],3,nchar(filelist[[1]][6]))
  meansdat_pr$wg = filelist[[1]][4]
  meansdat_pr$var = filelist[[1]][3]
  meansdat_pr$WUdomain = WUdomain
  meansdat_pr$SAdomain = SAdomain
  meansdat_all = rbind(meansdat_all,meansdat_pr)
}

######

meansdat_all$region <- factor(meansdat_all$region,levels = c('full','new mexico','louisiana'),ordered = TRUE)
#meansdat_all$wg <- factor(meansdat_all$wg,levels = c('tmax','pr','mv'),ordered = TRUE)
meansdat_all$wg <- factor(meansdat_all$wg,levels = c('tmax','pr'),ordered = TRUE)
meansdat_all$WUdomain <- factor(meansdat_all$WUdomain,levels = c('full','new mexico','louisiana'),ordered = TRUE)
meansdat_all$SAdomain <- factor(meansdat_all$SAdomain,levels = c('full','new mexico','louisiana'),ordered = TRUE)

########
# BMA posteriors

meansdatBMA_all = NULL 

for(i in 1:length(BMApostfiles)){
  load(BMApostfiles[i])
  if(var=="tmax"){
    meansdatBMA_pr = meansdatBMA_tmax
  }
  filelist = strsplit(BMApostfiles[i],"_",fixed=TRUE)
  WUdomain = substr(filelist[[1]][6],3,nchar(filelist[[1]][6]))
  SAdomain = substr(filelist[[1]][7],3,(nchar(filelist[[1]][7])-6))
  meansdatBMA_pr$wg = filelist[[1]][5]
  meansdatBMA_pr$var = filelist[[1]][4]
  meansdatBMA_pr$WUdomain = WUdomain
  meansdatBMA_pr$SAdomain = SAdomain
  meansdatBMA_all = rbind(meansdatBMA_all,meansdatBMA_pr)
}

######

meansdatBMA_all$region <- factor(meansdatBMA_all$region,levels = c('full','new mexico','louisiana'),ordered = TRUE)
#meansdat_all$wg <- factor(meansdat_all$wg,levels = c('tmax','pr','mv'),ordered = TRUE)
meansdatBMA_all$wg <- factor(meansdatBMA_all$wg,levels = c('tmax','pr'),ordered = TRUE)
meansdatBMA_all$WUdomain <- factor(meansdatBMA_all$WUdomain,levels = c('full','new mexico','louisiana'),ordered = TRUE)
meansdatBMA_all$SAdomain <- factor(meansdatBMA_all$SAdomain,levels = c('full','new mexico','louisiana'),ordered = TRUE)

#######
load("/home/woot0002/DS_ind/manuscript1/GCMlist.Rdata")

GCMfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)

GCMprojfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)
LOCAprojfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)

LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)

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
    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      GCMhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      GCMhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      GCMhvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      GCMhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist[[i]],NA)
    } 
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon");
  nc_close(nctest)
}

sapply(GCMhvardatalist,mean,na.rm=TRUE)

### GCM projected change
GCMpvardatalist = list()
for(i in 1:length(GCMpgroup)){
  nctest = nc_open(GCMpgroup[i])
  idx = which(names(nctest$var)==ncvarname)
  if(length(idx)==0){
    idx = which(names(nctest$var)==varin)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)

    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      GCMpvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      GCMpvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      GCMpvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      GCMpvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist[[i]],NA)
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
    idx = which(names(nctest$var)==varin)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)

    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      LOCAhvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      LOCAhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      LOCAhvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      LOCAhvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist[[i]],NA)
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
    idx = which(names(nctest$var)==varin)
    if(length(idx)==0){
      idx = which(names(nctest$var)==varin)
    }
  }
  
  tmp = ncvar_get(nctest,nctest$var[[idx]]$name)

    if(var=="pr" | var=="pr50" | var=="tmax95" | var=="tmin32"){
      LOCApvardatalist[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      LOCApvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist[[i]],NA)
    } 
    if(var=="tmax" | var=="tmin"){
      LOCApvardatalist[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      LOCApvardatalist[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist[[i]],NA)
    } 
  
  #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
  if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
  nc_close(nctest)
}

sapply(LOCApvardatalist,mean,na.rm=TRUE)

GCMchange = LOCAchange = GCMproj = LOCAproj = GCMhist = LOCAhist = array(NA,dim=c(length(lon),ncol=length(lat),26))
OBS = LOCAhvardatalist[[27]]
for(i in 1:26){
  GCMchange[,,i] = GCMpvardatalist[[i]]-GCMhvardatalist[[i]]
  LOCAchange[,,i] = LOCApvardatalist[[i]]-LOCAhvardatalist[[i]]
  GCMproj[,,i] = GCMpvardatalist[[i]]
  LOCAproj[,,i] = LOCApvardatalist[[i]]
  GCMhist[,,i] = GCMhvardatalist[[i]]
  LOCAhist[,,i] = LOCAhvardatalist[[i]]
}

######
# state mask

statelist = c("louisiana","new mexico")
regionmasklist = list()
for(s in 1:length(statelist)){
  test = nc_open(paste("/home/woot0002/DS_ind/",statelist[s],"_mask.nc",sep=""))
  regionmasklist[[s]] = ncvar_get(test,"mask")
  nc_close(test)
}

memsdat = NULL

for(m in 1:26){
  
  GCMhtmp = GCMhvardatalist[[m]]
  LOCAhtmp = LOCAhvardatalist[[m]]
  GCMptmp = GCMpvardatalist[[m]]
  LOCAptmp = LOCApvardatalist[[m]]
  
  for(s in 1:3){
    if(s<3){
      regionidx = which(regionmasklist[[s]]==1)
      histmeans_GCM = mean(GCMhtmp[regionidx],na.rm=TRUE)
      histmeans_LOCA = mean(LOCAhtmp[regionidx],na.rm=TRUE)
      
      changemeans_GCM = mean(GCMptmp[regionidx]-GCMhtmp[regionidx],na.rm=TRUE)
      changemeans_LOCA = mean(LOCAptmp[regionidx]-LOCAhtmp[regionidx],na.rm=TRUE)
      
      obs = rep(mean(OBS[regionidx],na.rm=TRUE),2)
      region = rep(statelist[s],2)
      
      GCMbias = GCMhtmp[regionidx]-OBS[regionidx]
      LOCAbias = LOCAhtmp[regionidx]-OBS[regionidx]
      
      RMSE_GCM = sqrt(mean((GCMbias)^2,na.rm=TRUE))
      RMSE_LOCA = sqrt(mean((LOCAbias)^2,na.rm=TRUE))
      
    } else {
      histmeans_GCM = mean(GCMhtmp,na.rm=TRUE)
      histmeans_LOCA = mean(LOCAhtmp,na.rm=TRUE)
      
      changemeans_GCM = mean(GCMptmp-GCMhtmp,na.rm=TRUE)
      changemeans_LOCA = mean(LOCAptmp-LOCAhtmp,na.rm=TRUE)
      
      obs = rep(mean(OBS,na.rm=TRUE),2)
      region = rep("full",2)
      
      GCMbias = GCMhtmp-OBS
      LOCAbias = LOCAhtmp-OBS
      
      RMSE_GCM = sqrt(mean((GCMbias)^2,na.rm=TRUE))
      RMSE_LOCA = sqrt(mean((LOCAbias)^2,na.rm=TRUE))
    }
    
    histmeans = c(histmeans_GCM,histmeans_LOCA)
    changemeans = c(changemeans_GCM,changemeans_LOCA)
    bias = histmeans-obs
    rmse = c(RMSE_GCM,RMSE_LOCA)
    GCM = rep(GCMlist[m],2)
    DS = c("CMIP5","LOCA")
    memsframe = data.frame(GCM,region,DS,histmeans,obs,changemeans,bias,rmse)
    memsdat = rbind(memsdat,memsframe)
  }
  message("GCM calcs finished for GCM ",m," / 26")
}

memsdat$region <- factor(memsdat$region,levels = c('full','new mexico','louisiana'),ordered = TRUE)

save(list=c("memsdat","meansdat_all"),file=paste("/home/woot0002/DS_ind/",var,"_datavals_projchange.Rdata",sep=""))

ggplot(memsdat, aes(x=region, y=histmeans, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=memsdat, mapping=aes(x=region,y=obs),size=2.5,shape=23,fill="blue") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("Hist Ens. Members")+xlab("")+ylab("mm")

ggplot(memsdat, aes(x=region, y=bias, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+geom_hline(yintercept=0,linetype="dashed") + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("Bias of Hist Ens. Members")+xlab("")+ylab("Bias (deg C)")

ggplot(memsdat, aes(x=region, y=rmse, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+geom_hline(yintercept=0,linetype="dashed") + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("RMSE of Hist Ens. Members")+xlab("")+ylab("RMSE deg C")

ggplot(memsdat, aes(x=region, y=changemeans, fill=DS)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+geom_hline(yintercept=0,linetype="dashed") + scale_fill_manual(values=c( "cyan", "yellow")) +
  ggtitle("Ens. Members Projected Change")+xlab("")+ylab("Change (mm)")


ggplot(meansdat_all, aes(x=region, y=bias))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+geom_hline(yintercept=0,linetype="dashed")+ggtitle("Bias of Ens. Means against Livneh")+xlab("Region to which the weighting is applied")+ylab("Bias (deg C)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))
ggplot(meansdat_all, aes(x=region, y=rmse))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+geom_hline(yintercept=0,linetype="dashed")+ggtitle("RMSE of Ens. Means against Livneh")+xlab("Region to which the weighting is applied")+ylab("RMSE (deg C)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))
ggplot(meansdat_all, aes(x=region, y=changemeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Mean Projected Changes - RCP8.5 2070-2099")+xlab("Region to which the weighting is applied")+ylab("Change (deg C)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))
ggplot(meansdat_all, aes(x=region, y=changemeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) + ggtitle("Ens. Mean Projected Changes - RCP8.5 2070-2099")+xlab("Region to which the weighting is applied")+ylab("Change (mm)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))

############
# plotting all

meansdat_CMIP5 = subset(meansdat_all,DS=="CMIP5")
meansdat_LOCA = subset(meansdat_all,DS=="LOCA")

meansdatBMA_CMIP5 = subset(meansdatBMA_all,DS=="CMIP5")
meansdatBMA_LOCA = subset(meansdatBMA_all,DS=="LOCA")

#+ylim(0,1600)

ggplot(meansdatBMA_CMIP5, aes(x=region, y=histmeans)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_CMIP5, mapping=aes(x=region,y=histmeans,colour = factor(group),shape=factor(DS)),size=2.5) + geom_point(data=meansdat_CMIP5, mapping=aes(x=region,y=obs),size=2.5,shape=23,fill="blue") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean")+xlab("")+ylab("Mean (deg C)") +facet_grid(cols=vars(WUdomain),rows=vars(wg))

ggplot(meansdatBMA_LOCA, aes(x=region, y=histmeans)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_LOCA, mapping=aes(x=region,y=histmeans,colour = factor(group),shape=factor(DS)),size=2.5) + geom_point(data=meansdat_LOCA, mapping=aes(x=region,y=obs),size=2.5,shape=23,fill="blue") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean")+xlab("")+ylab("Mean (deg C)") +facet_grid(cols=vars(WUdomain),rows=vars(wg))


ggplot(meansdatBMA_CMIP5, aes(x=region, y=bias)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_CMIP5, mapping=aes(x=region,y=bias,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean Bias")+xlab("")+ylab("Bias (deg C)") +facet_grid(cols=vars(WUdomain),rows=vars(wg))+geom_hline(yintercept=0,linetype="dashed")

ggplot(meansdatBMA_LOCA, aes(x=region, y=bias)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_LOCA, mapping=aes(x=region,y=bias,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean Bias")+xlab("")+ylab("Bias (deg C)") +facet_grid(cols=vars(WUdomain),rows=vars(wg))+geom_hline(yintercept=0,linetype="dashed")

ggplot(meansdatBMA_CMIP5, aes(x=region, y=rmse)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_CMIP5, mapping=aes(x=region,y=rmse,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean RMSE")+xlab("")+ylab("RMSE (mm)") +facet_grid(cols=vars(WUdomain),rows=vars(wg))+geom_hline(yintercept=0,linetype="dashed")
ggsave(paste("/home/woot0002/DS_ind/RMSE_CMIP5_",varin,".pdf",sep=""),device = "pdf",width=8,height=7)

ggplot(meansdatBMA_LOCA, aes(x=region, y=rmse)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_LOCA, mapping=aes(x=region,y=rmse,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Hist Ens. Mean RMSE")+xlab("")+ylab("RMSE (mm)") +facet_grid(cols=vars(WUdomain),rows=vars(wg))+geom_hline(yintercept=0,linetype="dashed")
ggsave(paste("/home/woot0002/DS_ind/RMSE_LOCA_",varin,".pdf",sep=""),device = "pdf",width=8,height=7)


ggplot(meansdatBMA_CMIP5, aes(x=region, y=changemeans)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_CMIP5, mapping=aes(x=region,y=changemeans,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Ens. Mean Projected Change")+xlab("")+ylab("Change (deg C)")+ylim(0,7)+facet_grid(cols=vars(WUdomain),rows=vars(wg))+geom_hline(yintercept=0,linetype="dashed")

#+ylim(0,6.5) 

ggplot(meansdatBMA_LOCA, aes(x=region, y=changemeans)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + geom_point(data=meansdat_LOCA, mapping=aes(x=region,y=changemeans,colour = factor(group),shape=factor(DS)),size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Ens. Mean Projected Change")+xlab("")+ylab("Change (deg C)")+ylim(0,7)+facet_grid(cols=vars(WUdomain),rows=vars(wg))+geom_hline(yintercept=0,linetype="dashed")


#
##########





ggplot(meansdat_CMIP5, aes(x=region, y=bias))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+geom_hline(yintercept=0,linetype="dashed")+ggtitle("Bias of Ens. Means against Livneh")+xlab("Region to which the weighting is applied")+ylab("Bias (mm)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))
ggplot(meansdat_CMIP5, aes(x=region, y=rmse))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+geom_hline(yintercept=0,linetype="dashed")+ggtitle("RMSE of Ens. Means against Livneh")+xlab("Region to which the weighting is applied")+ylab("RMSE (mm)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))
ggplot(meansdat_CMIP5, aes(x=region, y=changemeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Mean Projected Changes - RCP8.5 2070-2099")+xlab("Region to which the weighting is applied")+ylab("Change (mm)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))

ggplot(meansdat_LOCA, aes(x=region, y=bias))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+ggtitle("Bias of Ens. Means against Livneh")+xlab("Region to which the weighting is applied")+ylab("Bias (mm)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))
ggplot(meansdat_LOCA, aes(x=region, y=rmse))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+ggtitle("RMSE of Ens. Means against Livneh")+xlab("Region to which the weighting is applied")+ylab("RMSE (mm)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))
ggplot(meansdat_LOCA, aes(x=region, y=changemeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) + ggtitle("Ens. Mean Projected Changes - RCP8.5 2070-2099")+xlab("Region to which the weighting is applied")+ylab("Change (mm)")+facet_grid(cols=vars(WUdomain),rows=vars(wg))

