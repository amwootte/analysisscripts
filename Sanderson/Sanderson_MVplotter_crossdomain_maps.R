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
# Focusing Bias and change maps

WUs = c("full","louisiana","new mexico")
SAs = c("full","louisiana","new mexico")
varbased = c("pr","tmax","mv")
weightscheme = c("unweighted","skill","SI-h","SI-c")

GCM_bias_tmax = GCM_change_tmax = GCM_bias_pr = GCM_change_pr = array(NA,dim=c(190,140,length(varbased),length(WUs),length(SAs),length(weightscheme)))
LOCA_bias_tmax = LOCA_change_tmax = LOCA_bias_pr = LOCA_change_pr = array(NA,dim=c(190,140,length(varbased),length(WUs),length(SAs),length(weightscheme)))

#####

for(i in 1:length(varbased)){
  for(j in 1:length(WUs)){
    for(k in 1:length(SAs)){
      load(paste("/home/woot0002/DS_ind/Biasesbasedon_",varbased[i],"_WU",WUs[j],"_SA",WUs[k],"_weighting.Rdata",sep=""))
      load(paste("/home/woot0002/DS_ind/Changesbasedon_",varbased[i],"_WU",WUs[j],"_SA",WUs[k],"_weighting.Rdata",sep=""))
      for(l in 1:length(weightscheme)){
        if(l==1){
          GCM_bias_tmax[,,i,j,k,l]=GCMunweightedmean_bias_tmax
          LOCA_bias_tmax[,,i,j,k,l]=LOCAunweightedmean_bias_tmax
          GCM_bias_pr[,,i,j,k,l]=GCMunweightedmean_bias_pr
          LOCA_bias_pr[,,i,j,k,l]=LOCAunweightedmean_bias_pr
          GCM_change_tmax[,,i,j,k,l]=GCMunweightedmean_change_tmax
          LOCA_change_tmax[,,i,j,k,l]=LOCAunweightedmean_change_tmax
          GCM_change_pr[,,i,j,k,l]=GCMunweightedmean_change_pr
          LOCA_change_pr[,,i,j,k,l]=LOCAunweightedmean_change_pr
        }
        
        if(l==2){
          GCM_bias_tmax[,,i,j,k,l]=GCMskillmean_bias_tmax
          LOCA_bias_tmax[,,i,j,k,l]=LOCAskillmean_bias_tmax
          GCM_bias_pr[,,i,j,k,l]=GCMskillmean_bias_pr
          LOCA_bias_pr[,,i,j,k,l]=LOCAskillmean_bias_pr
          GCM_change_tmax[,,i,j,k,l]=GCMskillmean_change_tmax
          LOCA_change_tmax[,,i,j,k,l]=LOCAskillmean_change_tmax
          GCM_change_pr[,,i,j,k,l]=GCMskillmean_change_pr
          LOCA_change_pr[,,i,j,k,l]=LOCAskillmean_change_pr
        }
        
        if(l==3){
          GCM_bias_tmax[,,i,j,k,l]=GCMSIhmean_bias_tmax
          LOCA_bias_tmax[,,i,j,k,l]=LOCASIhmean_bias_tmax
          GCM_bias_pr[,,i,j,k,l]=GCMSIhmean_bias_pr
          LOCA_bias_pr[,,i,j,k,l]=LOCASIhmean_bias_pr
          GCM_change_tmax[,,i,j,k,l]=GCMSIhmean_change_tmax
          LOCA_change_tmax[,,i,j,k,l]=LOCASIhmean_change_tmax
          GCM_change_pr[,,i,j,k,l]=GCMSIhmean_change_pr
          LOCA_change_pr[,,i,j,k,l]=LOCASIhmean_change_pr
        }
        
        if(l==4){
          GCM_bias_tmax[,,i,j,k,l]=GCMSIcmean_bias_tmax
          LOCA_bias_tmax[,,i,j,k,l]=LOCASIcmean_bias_tmax
          GCM_bias_pr[,,i,j,k,l]=GCMSIcmean_bias_pr
          LOCA_bias_pr[,,i,j,k,l]=LOCASIcmean_bias_pr
          GCM_change_tmax[,,i,j,k,l]=GCMSIcmean_change_tmax
          LOCA_change_tmax[,,i,j,k,l]=LOCASIcmean_change_tmax
          GCM_change_pr[,,i,j,k,l]=GCMSIcmean_change_pr
          LOCA_change_pr[,,i,j,k,l]=LOCASIcmean_change_pr
        }
        
      } 
    }
  }
}

GCMfiles_pr = system("ls /home/woot0002/GCMs/regrid/pr_*histclimo*.nc",intern=TRUE)
nctest = nc_open(GCMfiles_pr[1])
lat = ncvar_get(nctest,"lat")
lon=ncvar_get(nctest,"lon")
nc_close(nctest)

# weight based on full domain
i=3 # weight based on pr
k=3 # weight applied to full domain

if(i==1){
  varbase = "pr"
}
if(i==2){
  varbase = "tmax"
}
if(i==3){
  varbase = "mv"
}

if(k==1){
  xlimrange = range(lon)
  ylimrange = range(lat)
  appto = "full"
}

if(k==2){
  xlimrange = c(-94.5,-89.5)
  ylimrange = c(28,34)
  appto = "louisiana"
}

if(k==3){
  xlimrange = c(-109.5,-102.5)
  ylimrange = c(30,39)
  appto = "new mexico"
}

biascolorbar_pr = colorramp(GCM_bias_pr,colorchoice="bluetored",Blimit=30,use_fixed_scale = TRUE, fixed_scale = c(-1200,1200))
biascolorbar_tmax = colorramp(GCM_bias_tmax,colorchoice="bluetored",Blimit=30,use_fixed_scale = FALSE, fixed_scale = c(-10,10))

changecolorbar_pr = colorramp(GCM_change_pr,colorchoice="browntogreen",Blimit=30,use_fixed_scale = TRUE,fixed_scale = c(-150,150))
changecolorbar_tmax = colorramp(GCM_change_tmax,colorchoice="bluetored",Blimit=30,use_fixed_scale = FALSE,fixed_scale = c(-6,6))

pdf(paste("/home/woot0002/DS_ind/CMIP5maps_prbias_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)

par(mfrow=c(3,4))

for(j in 1:3){ # weighting based on j domain
for(h in 1:4){ # weighting scheme used
    testsfc = list(x=lon,y=lat,z=GCM_bias_pr[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("GCM bias pr \n",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=biascolorbar_pr[[1]],col=biascolorbar_pr[[3]],breaks=biascolorbar_pr[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
}
}

dev.off()

pdf(paste("/home/woot0002/DS_ind/CMIP5maps_tmaxbias_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)

par(mfrow=c(3,4))
for(j in 1:3){
  for(h in 1:4){
    testsfc = list(x=lon,y=lat,z=GCM_bias_tmax[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("GCM bias tmax ",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=biascolorbar_tmax[[1]],col=biascolorbar_tmax[[3]],breaks=biascolorbar_tmax[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
  }
}
dev.off()

pdf(paste("/home/woot0002/DS_ind/CMIP5maps_prchange_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)

par(mfrow=c(3,4))
for(j in 1:3){
  for(h in 1:4){
    testsfc = list(x=lon,y=lat,z=GCM_change_pr[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("GCM change pr ",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=changecolorbar_pr[[1]],col=changecolorbar_pr[[3]],breaks=changecolorbar_pr[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
  }
}

dev.off()

pdf(paste("/home/woot0002/DS_ind/CMIP5maps_tmaxchange_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)

par(mfrow=c(3,4))
for(j in 1:3){
  for(h in 1:4){
    testsfc = list(x=lon,y=lat,z=GCM_change_tmax[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("GCM change tmax ",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=changecolorbar_tmax[[1]],col=changecolorbar_tmax[[3]],breaks=changecolorbar_tmax[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
  }
}

dev.off()





biascolorbar_pr = colorramp(LOCA_bias_pr,colorchoice="bluetored",Blimit=30,use_fixed_scale = TRUE, fixed_scale = c(-200,200))
biascolorbar_tmax = colorramp(LOCA_bias_tmax,colorchoice="bluetored",Blimit=30,use_fixed_scale = FALSE, fixed_scale = c(-10,10))

changecolorbar_pr = colorramp(LOCA_change_pr,colorchoice="browntogreen",Blimit=30,use_fixed_scale = TRUE,fixed_scale = c(-150,150))
changecolorbar_tmax = colorramp(LOCA_change_tmax,colorchoice="bluetored",Blimit=20,use_fixed_scale = FALSE,fixed_scale = c(-6,6))

pdf(paste("/home/woot0002/DS_ind/LOCAmaps_prbias_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)
par(mfrow=c(3,4))
for(j in 1:3){
  for(h in 1:4){
    testsfc = list(x=lon,y=lat,z=LOCA_bias_pr[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("LOCA bias pr ",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=biascolorbar_pr[[1]],col=biascolorbar_pr[[3]],breaks=biascolorbar_pr[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
  }
}
dev.off()

pdf(paste("/home/woot0002/DS_ind/LOCAmaps_tmaxbias_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)
par(mfrow=c(3,4))
for(j in 1:3){
  for(h in 1:4){
    testsfc = list(x=lon,y=lat,z=LOCA_bias_tmax[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("LOCA bias tmax ",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=biascolorbar_tmax[[1]],col=biascolorbar_tmax[[3]],breaks=biascolorbar_tmax[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
  }
}
dev.off()

pdf(paste("/home/woot0002/DS_ind/LOCAmaps_prchange_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)
par(mfrow=c(3,4))
for(j in 1:3){
  for(h in 1:4){
    testsfc = list(x=lon,y=lat,z=LOCA_change_pr[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("LOCA change pr ",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=changecolorbar_pr[[1]],col=changecolorbar_pr[[3]],breaks=changecolorbar_pr[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
  }
}
dev.off()

pdf(paste("/home/woot0002/DS_ind/LOCAmaps_tmaxchange_basedon-",varbase,"_appto-",appto,".pdf",sep=""),width=20,height=15)
par(mfrow=c(3,4))
for(j in 1:3){
  for(h in 1:4){
    testsfc = list(x=lon,y=lat,z=LOCA_change_tmax[,,i,j,k,h])
    surface(testsfc,type="I",main=paste("LOCA change tmax ",weightscheme[h]," based on ",WUs[j],sep=""),ylim=ylimrange,xlim=xlimrange,zlim=changecolorbar_tmax[[1]],col=changecolorbar_tmax[[3]],breaks=changecolorbar_tmax[[2]],xlab="Longitude",ylab="Latitude")
    map("state",add=TRUE)
  }
}
dev.off()

