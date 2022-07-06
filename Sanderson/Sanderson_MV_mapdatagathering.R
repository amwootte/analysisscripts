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

var = varin = "pr"
type="ann"

WUs=SAs=c("full","louisiana","new mexico")

GCM_bias_tmax = GCM_bias_pr = LOCA_bias_tmax = LOCA_bias_pr = GCM_change_tmax = GCM_change_pr = LOCA_change_tmax = LOCA_change_pr = array(NA,dim=c(190,140,5,3,3))

for(w in 1:3){
  for(s in 1:3){
    weightingused = WUs[w]
    stateapplied = SAs[s]

    if(weightingused=="full"){
      load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,".Rdata",sep=""))
      BMAweights_GCM = read.table(paste("best_BMA_combo_",var,".txt",sep=""))
      BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,".txt",sep=""))
    } else {
      load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,"_",weightingused,".Rdata",sep=""))
      BMAweights_GCM = read.table(paste("best_BMA_combo_",var,"_",weightingused,".txt",sep=""))
      BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
    }

    GCMhdat$BMA = t(BMAweights_GCM)[,1]
    LOCAhdat$BMA = t(BMAweights_LOCA)[,1]

#####
# get domain mask

    if(stateapplied!="full"){
      test = nc_open(paste("/home/woot0002/DS_ind/",stateapplied,"_mask.nc",sep=""))
      regionmask = ncvar_get(test,"mask")
      lon = ncvar_get(test,"lon")
      lat = ncvar_get(test,"lat")
      nc_close(test)
    }

####

    GCMweights= GCMhdat
    LOCAweights = LOCAhdat

    # precip files
    GCMfiles_pr = system("ls /home/woot0002/GCMs/regrid/pr_*histclimo*.nc",intern=TRUE)
    LOCAfiles_pr = system("ls /home/woot0002/LOCA/regrid/pr_*histclimo*.nc",intern=TRUE)
    GCMprojfiles_pr = system("ls /home/woot0002/GCMs/regrid/pr_*projclimo*.nc",intern=TRUE)
    LOCAprojfiles_pr = system("ls /home/woot0002/LOCA/regrid/pr_*projclimo*.nc",intern=TRUE)
    LIVNEHfile_pr = system("ls /home/woot0002/monthlyclimo/pr_day*livneh*.nc",intern=TRUE)

    # tasmax files
    GCMfiles_tmax = system("ls /home/woot0002/GCMs/regrid/tasmax_*histclimo*.nc",intern=TRUE)
    LOCAfiles_tmax = system("ls /home/woot0002/LOCA/regrid/tasmax_*histclimo*.nc",intern=TRUE)
    GCMprojfiles_tmax = system("ls /home/woot0002/GCMs/regrid/tasmax_*projclimo*.nc",intern=TRUE)
    LOCAprojfiles_tmax = system("ls /home/woot0002/LOCA/regrid/tasmax_*projclimo*.nc",intern=TRUE)
    LIVNEHfile_tmax = system("ls /home/woot0002/monthlyclimo/tasmax_day*livneh*.nc",intern=TRUE)

    # subset files down
    load("/home/woot0002/DS_ind/manuscript1/GCMlist.Rdata")

    GCM_hfiles_pr = GCM_pfiles_pr = LOCA_hfiles_pr = LOCA_pfiles_pr = c()
    GCM_hfiles_tmax = GCM_pfiles_tmax = LOCA_hfiles_tmax = LOCA_pfiles_tmax = c()

    for(i in 1:length(GCMlist)){
      #pr 
      GCM_hfiles_pr[i] = GCMfiles_pr[grep(paste(GCMlist[i],"_",sep=""),GCMfiles_pr)]
      GCM_pfiles_pr[i] = GCMprojfiles_pr[grep(paste(GCMlist[i],"_",sep=""),GCMprojfiles_pr)]
      LOCA_hfiles_pr[i] = LOCAfiles_pr[grep(paste(GCMlist[i],"_",sep=""),LOCAfiles_pr)]
      LOCA_pfiles_pr[i] = LOCAprojfiles_pr[grep(paste(GCMlist[i],"_",sep=""),LOCAprojfiles_pr)]
      #tmax
      GCM_hfiles_tmax[i] = GCMfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),GCMfiles_tmax)]
      GCM_pfiles_tmax[i] = GCMprojfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),GCMprojfiles_tmax)]
      LOCA_hfiles_tmax[i] = LOCAfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),LOCAfiles_tmax)]
      LOCA_pfiles_tmax[i] = LOCAprojfiles_tmax[grep(paste(GCMlist[i],"_",sep=""),LOCAprojfiles_tmax)]
    }

    ###
    # create full filelist + metadata table - historical

    #GCMs
    filelist1 = do.call("rbind",strsplit(GCM_hfiles_pr,"/",fixed=TRUE))
    filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
    filelist2 = as.data.frame(filelist2)
    filelist2$training = "NA"
    GCMhdat = filelist2[,c(2,3,4,6)]
    names(GCMhdat) = c("GCM","exp","DS","training")

    #LOCA
    filelist1 = do.call("rbind",strsplit(LOCA_hfiles_pr,"/",fixed=TRUE))
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
    GCMgroup_pr = c(GCM_hfiles_pr,LIVNEHfile_pr)
    LOCAgroup_pr = c(LOCA_hfiles_pr,LIVNEHfile_pr)

    GCMgroup_tmax = c(GCM_hfiles_tmax,LIVNEHfile_tmax)
    LOCAgroup_tmax = c(LOCA_hfiles_tmax,LIVNEHfile_tmax)

    ###
    # create full filelist + metadata table - projected

    #GCMs
    filelist1 = do.call("rbind",strsplit(GCM_pfiles_pr,"/",fixed=TRUE))
    filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
    filelist2 = as.data.frame(filelist2)
    filelist2$training = "NA"
    GCMpdat = filelist2[,c(2,3,4,6)]
    names(GCMpdat) = c("GCM","exp","DS","training")

    #LOCA
    filelist1 = do.call("rbind",strsplit(LOCA_pfiles_pr,"/",fixed=TRUE))
    filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
    filelist2 = as.data.frame(filelist2)
    filelist2$training = "Livneh"
    LOCApdat = filelist2[,c(2,3,4,6)]
    names(LOCApdat) = names(GCMpdat)

    # all files
    GCMpgroup_pr = GCM_pfiles_pr
    LOCApgroup_pr = LOCA_pfiles_pr
    GCMpgroup_tmax = GCM_pfiles_tmax
    LOCApgroup_tmax = LOCA_pfiles_tmax

    ######
    # Gather data

    ncvarname = "prclimo"
    ### GCM hist + Livneh - pr
    GCMhvardatalist_pr = list()
    for(i in 1:length(GCMgroup_pr)){
      nctest = nc_open(GCMgroup_pr[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      if(w==1 & s==1 & i==1){
        lon = ncvar_get(nctest,"lon")
        lat = ncvar_get(nctest,"lat")
      }
      GCMhvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(stateapplied=="full"){
        GCMhvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist_pr[[i]],NA)
      } else {
        GCMhvardatalist_pr[[i]] = ifelse(regionmask==1,GCMhvardatalist_pr[[i]],NA)
      }
      GCMhvardatalist_pr[[i]] = ifelse(GCMhvardatalist_pr[[i]]==0,NA,GCMhvardatalist_pr[[i]])
      #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon");
      nc_close(nctest)
    }

    sapply(GCMhvardatalist_pr,mean,na.rm=TRUE)

    ### GCM projected change - pr
    GCMpvardatalist_pr = list()
    for(i in 1:length(GCMpgroup_pr)){
      nctest = nc_open(GCMpgroup_pr[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      GCMpvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(stateapplied=="full"){
        GCMpvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist_pr[[i]],NA)
      } else {
        GCMpvardatalist_pr[[i]] = ifelse(regionmask==1,GCMpvardatalist_pr[[i]],NA)
      }
      #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
      nc_close(nctest)
    }

    sapply(GCMpvardatalist_pr,mean,na.rm=TRUE)

    ### LOCA historical + Livneh - pr

    LOCAhvardatalist_pr = list()
    for(i in 1:length(LOCAgroup_pr)){
      nctest = nc_open(LOCAgroup_pr[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      LOCAhvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(stateapplied=="full"){
        LOCAhvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist_pr[[i]],NA)
      } else{
        LOCAhvardatalist_pr[[i]] = ifelse(regionmask==1,LOCAhvardatalist_pr[[i]],NA)
      }
      LOCAhvardatalist_pr[[i]] = ifelse(LOCAhvardatalist_pr[[i]]==0,NA,LOCAhvardatalist_pr[[i]])
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
      nc_close(nctest)
    }

    sapply(LOCAhvardatalist_pr,mean,na.rm=TRUE)

    ### LOCA projected change - pr

    LOCApvardatalist_pr = list()
    for(i in 1:length(LOCApgroup_pr)){
      nctest = nc_open(LOCApgroup_pr[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      LOCApvardatalist_pr[[i]] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      if(stateapplied=="full"){
        LOCApvardatalist_pr[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist_pr[[i]],NA)
      } else{
        LOCApvardatalist_pr[[i]] = ifelse(regionmask==1,LOCApvardatalist_pr[[i]],NA)
      }
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
      nc_close(nctest)
    }

    sapply(LOCApvardatalist_pr,mean,na.rm=TRUE)

    ######
    # Gather Data 2

    ncvarname = "tmaxclimo"
    ### GCM hist + Livneh - tmax
    GCMhvardatalist_tmax = list()
    for(i in 1:length(GCMgroup_tmax)){
      nctest = nc_open(GCMgroup_tmax[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      GCMhvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(stateapplied=="full"){
        GCMhvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMhvardatalist_tmax[[i]],NA)
      } else{
        GCMhvardatalist_tmax[[i]] = ifelse(regionmask==1,GCMhvardatalist_tmax[[i]],NA)
      }
      #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon");
      nc_close(nctest)
    }

    sapply(GCMhvardatalist_tmax,mean,na.rm=TRUE)

    ### GCM projected change - tmax
    GCMpvardatalist_tmax = list()
    for(i in 1:length(GCMpgroup_tmax)){
      nctest = nc_open(GCMpgroup_tmax[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      GCMpvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(stateapplied=="full"){
        GCMpvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,GCMpvardatalist_tmax[[i]],NA)
      } else{
        GCMpvardatalist_tmax[[i]] = ifelse(regionmask==1,GCMpvardatalist_tmax[[i]],NA)
      }
      #vardatalist[[i]] = ncvar_get(nctest,ncvarname)
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
      nc_close(nctest)
    }

    sapply(GCMpvardatalist_tmax,mean,na.rm=TRUE)

    ### LOCA historical + Livneh - tmax

    LOCAhvardatalist_tmax = list()
    for(i in 1:length(LOCAgroup_tmax)){
      nctest = nc_open(LOCAgroup_tmax[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      LOCAhvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(stateapplied=="full"){
        LOCAhvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCAhvardatalist_tmax[[i]],NA)
      } else{
        LOCAhvardatalist_tmax[[i]] = ifelse(regionmask==1,LOCAhvardatalist_tmax[[i]],NA)
      }
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
      nc_close(nctest)
    }

    sapply(LOCAhvardatalist_tmax,mean,na.rm=TRUE)

    ### LOCA projected change - tmax

    LOCApvardatalist_tmax = list()
    for(i in 1:length(LOCApgroup_tmax)){
      nctest = nc_open(LOCApgroup_tmax[i])
      idx = which(names(nctest$var)==ncvarname)
      tmp = ncvar_get(nctest,nctest$var[[idx]]$name)
      LOCApvardatalist_tmax[[i]] = apply(tmp,c(1,2),mean,na.rm=TRUE)
      if(stateapplied=="full"){
        LOCApvardatalist_tmax[[i]] = ifelse(is.na(tmp[,,1])==FALSE,LOCApvardatalist_tmax[[i]],NA)
      } else{
        LOCApvardatalist_tmax[[i]] = ifelse(regionmask==1,LOCApvardatalist_tmax[[i]],NA)
      }
      if(i==1) lat = ncvar_get(nctest,"lat"); lon=ncvar_get(nctest,"lon")
      nc_close(nctest)
    }

    sapply(LOCApvardatalist_tmax,mean,na.rm=TRUE)


    #######
    # projected changes - _pr

    GCMchange_pr = LOCAchange_pr = GCMproj_pr = LOCAproj_pr = GCMhist_pr = LOCAhist_pr = array(NA,dim=c(length(lon),ncol=length(lat),26))
    OBS_pr = LOCAhvardatalist_pr[[27]]
    if(var=="pr"){
      if(min(OBS_pr,na.rm=TRUE)==0){
      OBS_pr = ifelse(OBS_pr==0,NA,OBS_pr)
      }
    }

    for(i in 1:26){
      GCMchange_pr[,,i] = GCMpvardatalist_pr[[i]]-GCMhvardatalist_pr[[i]]
      LOCAchange_pr[,,i] = LOCApvardatalist_pr[[i]]-LOCAhvardatalist_pr[[i]]
      GCMproj_pr[,,i] = GCMpvardatalist_pr[[i]]
      LOCAproj_pr[,,i] = LOCApvardatalist_pr[[i]]
      GCMhist_pr[,,i] = GCMhvardatalist_pr[[i]]
      LOCAhist_pr[,,i] = LOCAhvardatalist_pr[[i]]
    }

    #######
    # projected changes - _tmax

    GCMchange_tmax = LOCAchange_tmax = GCMproj_tmax = LOCAproj_tmax = GCMhist_tmax = LOCAhist_tmax = array(NA,dim=c(length(lon),ncol=length(lat),26))
    OBS_tmax = LOCAhvardatalist_tmax[[27]]
    for(i in 1:26){
      GCMchange_tmax[,,i] = GCMpvardatalist_tmax[[i]]-GCMhvardatalist_tmax[[i]]
      LOCAchange_tmax[,,i] = LOCApvardatalist_tmax[[i]]-LOCAhvardatalist_tmax[[i]]
      GCMproj_tmax[,,i] = GCMpvardatalist_tmax[[i]]
      LOCAproj_tmax[,,i] = LOCApvardatalist_tmax[[i]]
      GCMhist_tmax[,,i] = GCMhvardatalist_tmax[[i]]
      LOCAhist_tmax[,,i] = LOCAhvardatalist_tmax[[i]]
    }


    ######
    # prep weights

    GCMweights$Wh = (GCMweights$Wuh*GCMweights$Wqh)/sum(GCMweights$Wuh*GCMweights$Wqh)
    GCMweights$Wc = (GCMweights$Wuc*GCMweights$Wqc)/sum(GCMweights$Wuc*GCMweights$Wqc)
    GCMweights$Ws = GCMweights$Wqh/sum(GCMweights$Wqh)

    LOCAweights$Wh = (LOCAweights$Wuh*LOCAweights$Wqh)/sum(LOCAweights$Wuh*LOCAweights$Wqh)
    LOCAweights$Wc = (LOCAweights$Wuc*LOCAweights$Wqc)/sum(LOCAweights$Wuc*LOCAweights$Wqc)
    LOCAweights$Ws = LOCAweights$Wqh/sum(LOCAweights$Wqh)

    ######
    # Calculate historical means (weighted and unweighted) - pr

    GCMunweightedmean_hist_pr = apply(GCMhist_pr,c(1,2),mean,na.rm=TRUE)
    LOCAunweightedmean_hist_pr = apply(LOCAhist_pr,c(1,2),mean,na.rm=TRUE)

    GCMskillmean_hist_pr = GCMSIhmean_hist_pr = GCMSIcmean_hist_pr = GCMBMAmean_hist_pr = GCMunweightedmean_hist_pr
    LOCAskillmean_hist_pr = LOCASIhmean_hist_pr = LOCASIcmean_hist_pr = LOCABMAmean_hist_pr = LOCAunweightedmean_hist_pr

    for(i in 1:26){
  
      ## skill mean
      tmpG = GCMhist_pr[,,i]*GCMweights$Ws[i]
      tmpL = LOCAhist_pr[,,i]*LOCAweights$Ws[i]
      if(i==1){
        GCMskillmean_hist_pr = tmpG
        LOCAskillmean_hist_pr = tmpL
      } else {
        GCMskillmean_hist_pr = GCMskillmean_hist_pr+tmpG
        LOCAskillmean_hist_pr = LOCAskillmean_hist_pr+tmpL
      }
  
      ## skill+ind hist only
      tmpG = GCMhist_pr[,,i]*GCMweights$Wh[i]
      tmpL = LOCAhist_pr[,,i]*LOCAweights$Wh[i]
      if(i==1){
        GCMSIhmean_hist_pr = tmpG
        LOCASIhmean_hist_pr = tmpL
      } else {
        GCMSIhmean_hist_pr = GCMSIhmean_hist_pr+tmpG
        LOCASIhmean_hist_pr = LOCASIhmean_hist_pr+tmpL
      }
  
      ## skill+ind hist and change
      tmpG = GCMhist_pr[,,i]*GCMweights$Wc[i]
      tmpL = LOCAhist_pr[,,i]*LOCAweights$Wc[i]
      if(i==1){
        GCMSIcmean_hist_pr = tmpG
        LOCASIcmean_hist_pr = tmpL
      } else {
        GCMSIcmean_hist_pr = GCMSIcmean_hist_pr+tmpG
        LOCASIcmean_hist_pr = LOCASIcmean_hist_pr+tmpL
      }
  
      ## BMA hist and change
      tmpG = GCMhist_pr[,,i]*GCMweights$BMA[i]
      tmpL = LOCAhist_pr[,,i]*LOCAweights$BMA[i]
      if(i==1){
        GCMBMAmean_hist_pr = tmpG
        LOCABMAmean_hist_pr = tmpL
      } else {
        GCMBMAmean_hist_pr = GCMBMAmean_hist_pr+tmpG
        LOCABMAmean_hist_pr = LOCABMAmean_hist_pr+tmpL
      }
  
    }

    LOCA_bias_pr[,,1,w,s] = LOCAunweightedmean_hist_pr-OBS_pr
    LOCA_bias_pr[,,2,w,s] = LOCAskillmean_hist_pr-OBS_pr
    LOCA_bias_pr[,,3,w,s] = LOCASIhmean_hist_pr-OBS_pr
    LOCA_bias_pr[,,4,w,s] = LOCASIcmean_hist_pr-OBS_pr
    LOCA_bias_pr[,,5,w,s] = LOCABMAmean_hist_pr-OBS_pr

    GCM_bias_pr[,,1,w,s] = GCMunweightedmean_hist_pr-OBS_pr
    GCM_bias_pr[,,2,w,s] = GCMskillmean_hist_pr-OBS_pr
    GCM_bias_pr[,,3,w,s] = GCMSIhmean_hist_pr-OBS_pr
    GCM_bias_pr[,,4,w,s] = GCMSIcmean_hist_pr-OBS_pr
    GCM_bias_pr[,,5,w,s] = GCMBMAmean_hist_pr-OBS_pr

    ######
    # Calculate historical means (weighted and unweighted) - tmax

    GCMunweightedmean_hist_tmax = apply(GCMhist_tmax,c(1,2),mean,na.rm=TRUE)
    LOCAunweightedmean_hist_tmax = apply(LOCAhist_tmax,c(1,2),mean,na.rm=TRUE)

    GCMskillmean_hist_tmax = GCMSIhmean_hist_tmax = GCMSIcmean_hist_tmax = GCMBMAmean_hist_tmax = GCMunweightedmean_hist_tmax
    LOCAskillmean_hist_tmax = LOCASIhmean_hist_tmax = LOCASIcmean_hist_tmax = LOCABMAmean_hist_tmax = LOCAunweightedmean_hist_tmax

    for(i in 1:26){
  
      ## skill mean
      tmpG = GCMhist_tmax[,,i]*GCMweights$Ws[i]
      tmpL = LOCAhist_tmax[,,i]*LOCAweights$Ws[i]
      if(i==1){
        GCMskillmean_hist_tmax = tmpG
        LOCAskillmean_hist_tmax = tmpL
      } else {
        GCMskillmean_hist_tmax = GCMskillmean_hist_tmax+tmpG
        LOCAskillmean_hist_tmax = LOCAskillmean_hist_tmax+tmpL
      }
  
      ## skill+ind hist only
      tmpG = GCMhist_tmax[,,i]*GCMweights$Wh[i]
      tmpL = LOCAhist_tmax[,,i]*LOCAweights$Wh[i]
      if(i==1){
        GCMSIhmean_hist_tmax = tmpG
        LOCASIhmean_hist_tmax = tmpL
      } else {
        GCMSIhmean_hist_tmax = GCMSIhmean_hist_tmax+tmpG
        LOCASIhmean_hist_tmax = LOCASIhmean_hist_tmax+tmpL
      }
  
      ## skill+ind hist and change
      tmpG = GCMhist_tmax[,,i]*GCMweights$Wc[i]
      tmpL = LOCAhist_tmax[,,i]*LOCAweights$Wc[i]
      if(i==1){
        GCMSIcmean_hist_tmax = tmpG
        LOCASIcmean_hist_tmax = tmpL
      } else {
        GCMSIcmean_hist_tmax = GCMSIcmean_hist_tmax+tmpG
        LOCASIcmean_hist_tmax = LOCASIcmean_hist_tmax+tmpL
      }
  
      ## BMA hist and change
      tmpG = GCMhist_tmax[,,i]*GCMweights$BMA[i]
      tmpL = LOCAhist_tmax[,,i]*LOCAweights$BMA[i]
      if(i==1){
        GCMBMAmean_hist_tmax = tmpG
        LOCABMAmean_hist_tmax = tmpL
      } else {
        GCMBMAmean_hist_tmax = GCMBMAmean_hist_tmax+tmpG
        LOCABMAmean_hist_tmax = LOCABMAmean_hist_tmax+tmpL
      }
  
    }


    LOCA_bias_tmax[,,1,w,s] = LOCAunweightedmean_hist_tmax-OBS_tmax
    LOCA_bias_tmax[,,2,w,s] = LOCAskillmean_hist_tmax-OBS_tmax
    LOCA_bias_tmax[,,3,w,s] = LOCASIhmean_hist_tmax-OBS_tmax
    LOCA_bias_tmax[,,4,w,s] = LOCASIcmean_hist_tmax-OBS_tmax
    LOCA_bias_tmax[,,5,w,s] = LOCABMAmean_hist_tmax-OBS_tmax

    GCM_bias_tmax[,,1,w,s] = GCMunweightedmean_hist_tmax-OBS_tmax
    GCM_bias_tmax[,,2,w,s] = GCMskillmean_hist_tmax-OBS_tmax
    GCM_bias_tmax[,,3,w,s] = GCMSIhmean_hist_tmax-OBS_tmax
    GCM_bias_tmax[,,4,w,s] = GCMSIcmean_hist_tmax-OBS_tmax
    GCM_bias_tmax[,,5,w,s] = GCMBMAmean_hist_tmax-OBS_tmax

    ######
    # Calculate change means (weighted and unweighted) - pr only

    GCMunweightedmean_change_pr = apply(GCMchange_pr,c(1,2),mean,na.rm=TRUE)
    LOCAunweightedmean_change_pr = apply(LOCAchange_pr,c(1,2),mean,na.rm=TRUE)

    GCMskillmean_change_pr = GCMSIhmean_change_pr = GCMSIcmean_change_pr = GCMBMAmean_change_pr = GCMunweightedmean_change_pr
    LOCAskillmean_change_pr = LOCASIhmean_change_pr = LOCASIcmean_change_pr = LOCABMAmean_change_pr = LOCAunweightedmean_change_pr

    for(i in 1:26){
  
      ## skill mean
      tmpG = GCMchange_pr[,,i]*GCMweights$Ws[i]
      tmpL = LOCAchange_pr[,,i]*LOCAweights$Ws[i]
      if(i==1){
        GCMskillmean_change_pr = tmpG
        LOCAskillmean_change_pr = tmpL
      } else {
        GCMskillmean_change_pr = GCMskillmean_change_pr+tmpG
        LOCAskillmean_change_pr = LOCAskillmean_change_pr+tmpL
      }
  
      ## skill+ind hist only
      tmpG = GCMchange_pr[,,i]*GCMweights$Wh[i]
      tmpL = LOCAchange_pr[,,i]*LOCAweights$Wh[i]
      if(i==1){
        GCMSIhmean_change_pr = tmpG
        LOCASIhmean_change_pr = tmpL
      } else {
        GCMSIhmean_change_pr = GCMSIhmean_change_pr+tmpG
        LOCASIhmean_change_pr = LOCASIhmean_change_pr+tmpL
      }
  
      ## skill+ind hist and change
      tmpG = GCMchange_pr[,,i]*GCMweights$Wc[i]
      tmpL = LOCAchange_pr[,,i]*LOCAweights$Wc[i]
      if(i==1){
        GCMSIcmean_change_pr = tmpG
        LOCASIcmean_change_pr = tmpL
      } else {
        GCMSIcmean_change_pr = GCMSIcmean_change_pr+tmpG
        LOCASIcmean_change_pr = LOCASIcmean_change_pr+tmpL
      }
  
      ## BMA hist and change
      tmpG = GCMchange_pr[,,i]*GCMweights$BMA[i]
      tmpL = LOCAchange_pr[,,i]*LOCAweights$BMA[i]
      if(i==1){
        GCMBMAmean_change_pr = tmpG
        LOCABMAmean_change_pr = tmpL
      } else {
        GCMBMAmean_change_pr = GCMBMAmean_change_pr+tmpG
        LOCABMAmean_change_pr = LOCABMAmean_change_pr+tmpL
      }
  
    }

    LOCA_change_pr[,,1,w,s] = LOCAunweightedmean_change_pr
    LOCA_change_pr[,,2,w,s] = LOCAskillmean_change_pr
    LOCA_change_pr[,,3,w,s] = LOCASIhmean_change_pr
    LOCA_change_pr[,,4,w,s] = LOCASIcmean_change_pr
    LOCA_change_pr[,,5,w,s] = LOCABMAmean_change_pr

    GCM_change_pr[,,1,w,s] = GCMunweightedmean_change_pr
    GCM_change_pr[,,2,w,s] = GCMskillmean_change_pr
    GCM_change_pr[,,3,w,s] = GCMSIhmean_change_pr
    GCM_change_pr[,,4,w,s] = GCMSIcmean_change_pr
    GCM_change_pr[,,5,w,s] = GCMBMAmean_change_pr


    ######
    # Calculate change means (weighted and unweighted) - tmax only

    GCMunweightedmean_change_tmax = apply(GCMchange_tmax,c(1,2),mean,na.rm=TRUE)
    LOCAunweightedmean_change_tmax = apply(LOCAchange_tmax,c(1,2),mean,na.rm=TRUE)

    GCMskillmean_change_tmax = GCMSIhmean_change_tmax = GCMSIcmean_change_tmax = GCMBMAmean_change_tmax = GCMunweightedmean_change_tmax
    LOCAskillmean_change_tmax = LOCASIhmean_change_tmax = LOCASIcmean_change_tmax = LOCABMAmean_change_tmax = LOCAunweightedmean_change_tmax

    for(i in 1:26){
  
      ## skill mean
      tmpG = GCMchange_tmax[,,i]*GCMweights$Ws[i]
      tmpL = LOCAchange_tmax[,,i]*LOCAweights$Ws[i]
      if(i==1){
        GCMskillmean_change_tmax = tmpG
        LOCAskillmean_change_tmax = tmpL
      } else {
        GCMskillmean_change_tmax = GCMskillmean_change_tmax+tmpG
        LOCAskillmean_change_tmax = LOCAskillmean_change_tmax+tmpL
      }
  
      ## skill+ind hist only
      tmpG = GCMchange_tmax[,,i]*GCMweights$Wh[i]
      tmpL = LOCAchange_tmax[,,i]*LOCAweights$Wh[i]
      if(i==1){
        GCMSIhmean_change_tmax = tmpG
        LOCASIhmean_change_tmax = tmpL
      } else {
        GCMSIhmean_change_tmax = GCMSIhmean_change_tmax+tmpG
        LOCASIhmean_change_tmax = LOCASIhmean_change_tmax+tmpL
      }
  
      ## skill+ind hist and change
      tmpG = GCMchange_tmax[,,i]*GCMweights$Wc[i]
      tmpL = LOCAchange_tmax[,,i]*LOCAweights$Wc[i]
      if(i==1){
        GCMSIcmean_change_tmax = tmpG
        LOCASIcmean_change_tmax = tmpL
      } else {
        GCMSIcmean_change_tmax = GCMSIcmean_change_tmax+tmpG
        LOCASIcmean_change_tmax = LOCASIcmean_change_tmax+tmpL
      }
  
      ## BMA hist and change
      tmpG = GCMchange_tmax[,,i]*GCMweights$BMA[i]
      tmpL = LOCAchange_tmax[,,i]*LOCAweights$BMA[i]
      if(i==1){
        GCMBMAmean_change_tmax = tmpG
        LOCABMAmean_change_tmax = tmpL
      } else {
        GCMBMAmean_change_tmax = GCMBMAmean_change_tmax+tmpG
        LOCABMAmean_change_tmax = LOCABMAmean_change_tmax+tmpL
      }
    }


    LOCA_change_tmax[,,1,w,s] = LOCAunweightedmean_change_tmax
    LOCA_change_tmax[,,2,w,s] = LOCAskillmean_change_tmax
    LOCA_change_tmax[,,3,w,s] = LOCASIhmean_change_tmax
    LOCA_change_tmax[,,4,w,s] = LOCASIcmean_change_tmax
    LOCA_change_tmax[,,5,w,s] = LOCABMAmean_change_tmax

    GCM_change_tmax[,,1,w,s] = GCMunweightedmean_change_tmax
    GCM_change_tmax[,,2,w,s] = GCMskillmean_change_tmax
    GCM_change_tmax[,,3,w,s] = GCMSIhmean_change_tmax
    GCM_change_tmax[,,4,w,s] = GCMSIcmean_change_tmax
    GCM_change_tmax[,,5,w,s] = GCMBMAmean_change_tmax

    message("Finished Calcs for ",w," / 3 weights and ",s," / 3 applications")
  }
}
lon = lon-360

save(list=c("LOCA_change_tmax","LOCA_change_pr","GCM_change_tmax","GCM_change_pr","LOCA_bias_tmax","LOCA_bias_pr","GCM_bias_tmax","GCM_bias_pr","lon","lat","GCMweights","LOCAweights","WUs","SAs","var"),file=paste("/home/woot0002/DS_ind/valuesarrays_WU",var,".Rdata",sep=""))

