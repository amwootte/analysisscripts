###################
#
# R stippling test

library(ncdf4)
library(maps)
library(fields)
library(sp)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

# set basic arguments
varname = "rx1day" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2070,2099)
colorchoice = "browntoblue"
colorlimits = c(-20,20)
BINS = 30
USEFS = FALSE
typepart = "thres"
percbound = 20
ensthres = 18
threshold = 20

type=1
# type=1 - individual point
# type=2 - based on pdf

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

############

varin = varname
if(varname=="tmax85" | varname=="tmax90" | varname=="tmax95" | varname=="tmax100" | varname=="gsl") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="prcptot" | varname=="prmean" | varname=="prmed" | varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd" | varname=="R99" | varname=="R95" | varname=="R90"  | varname=="R99freq" | varname=="R95freq" | varname=="R90freq") varin="pr"

########
# Find all file names

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

projnotes = paste(projnotes,collapse=",")
histnotes = paste(histnotes,collapse=",")

histlist = paste(histfilelist,collapse=",")
projlist = paste(projfilelist,collapse=",")

###########

test = nc_open(paste("/data2/3to5/I35/all_mems/",varname,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_ann.nc",sep=""))
memhist = ncvar_get(test,"histmean")
memproj = ncvar_get(test,"projmeandiff")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

test = nc_open(paste("/data2/3to5/I35/ens_means/",varname,"_ensmean_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_ann.nc",sep=""))
meanproj85 = ncvar_get(test,"projmeandiff_rcp85")
meanproj45 = ncvar_get(test,"projmeandiff_rcp45")
meanproj26 = ncvar_get(test,"projmeandiff_rcp26")
nc_close(test)

rcp26idx = which(projfilebreakdown$scen=="rcp26")
rcp45idx = which(projfilebreakdown$scen=="rcp45")
rcp85idx = which(projfilebreakdown$scen=="rcp85")

  # RCP2.6 close to ens mean
  if(type==1){
    ensmean = apply(memproj[,,rcp26idx],c(1,2),mean,na.rm=TRUE)
    if(typepart=="chgper"){
    chgper = memproj[,,rcp26idx]
    for(i in 1:27){
      top = memproj[,,rcp26idx[i]]-ensmean
      bottom = abs(ensmean)
      chgper[,,i] = (top/bottom)*100
    }
    
    pb = c(-percbound,percbound)
    tmppos = ifelse(chgper>pb[2],1,0)
    tmpmid = ifelse(chgper<=pb[2] & chgper>= pb[1],1,0)
    tmpneg = ifelse(chgper< pb[1],1,0)
    }
    
    if(typepart=="sdbound"){
      sdmat = matrix(NA,nrow=length(lon),ncol=length(lat))
      for(R in 1:length(lon)){
        for(C in 1:length(lat)){
          if(all(is.na(memproj[R,C,rcp26idx])==FALSE)==TRUE){
            sdmat[R,C] = sd(memproj[R,C,rcp26idx],na.rm=TRUE)
          }
        }
      }
      
      upperbd = ensmean+1*sdmat
      lowerbd = ensmean-1*sdmat
      
      tmpmid = memproj[,,rcp26idx]
      
      for(i in 1:27){
        tmpmid[,,i] = ifelse(memproj[,,rcp26idx[i]]<=upperbd & memproj[,,rcp26idx[i]]>= lowerbd,1,0) 
      }
    }  
    
    if(typepart=="thres"){
      
      upperbd = ensmean+threshold
      lowerbd = ensmean-threshold
      
      tmpmid = memproj[,,rcp26idx]
      
      for(i in 1:27){
        tmpmid[,,i] = ifelse(memproj[,,rcp26idx[i]]<=upperbd & memproj[,,rcp26idx[i]]>= lowerbd,1,0) 
      }
    }  
    
    
    agreemid26 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
  }
  
  if(type==2){
    ensmean=meanproj26
    
    if(typepart=="chgper"){
    quants = quantile(memproj[,,rcp26idx],probs=seq(0.01,0.99,by=0.01),na.rm=TRUE)
    quantdiff = abs(quants-mean(ensmean,na.rm=TRUE))
    minidx=which(quantdiff==min(quantdiff))
    lwbd = minidx-percbound
    if(lwbd<1) lwbd=1
    upbd = minidx+percbound
    if(upbd>99) upbd=99
    
    quantrange = c(quantdiff[lwbd],quantdiff[upbd])
    
    upperbd = ensmean+quantrange[2]
    lowerbd = ensmean-quantrange[1]
    
    tmpmid = memproj[,,rcp26idx]
    for(i in 1:27){
      tmpmid[,,i] = ifelse(memproj[,,rcp26idx[i]]<=upperbd & memproj[,,rcp26idx[i]]>=lowerbd,1,0)  
    }
    }
    
    if(typepart=="sdbound"){
      sdval = sd(memproj[,,rcp26idx],na.rm=TRUE)
      upperbd = ensmean+1*sdval
      lowerbd = ensmean-1*sdval
      
      tmpmid = memproj[,,rcp26idx]
      for(i in 1:27){
        tmpmid[,,i] = ifelse(memproj[,,rcp26idx[i]]<=upperbd & memproj[,,rcp26idx[i]]>=lowerbd,1,0)  
      }
    }
    
    agreemid26 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
  }

  # RCP4.5 close to ens mean
if(type==1){  
ensmean = apply(memproj[,,rcp45idx],c(1,2),mean,na.rm=TRUE)
  
if(typepart=="chgper"){
chgper = memproj[,,rcp45idx]
  for(i in 1:27){
    top = memproj[,,rcp45idx[i]]-ensmean
    bottom = abs(ensmean)
    chgper[,,i] = (top/bottom)*100
  }
  
  tmppos = ifelse(chgper>pb[2],1,0)
  tmpmid = ifelse(chgper<=pb[2] & chgper>= pb[1],1,0)
  tmpneg = ifelse(chgper< pb[1],1,0)
}
  
  if(typepart=="sdbound"){
    sdmat = matrix(NA,nrow=length(lon),ncol=length(lat))
    for(R in 1:length(lon)){
      for(C in 1:length(lat)){
        if(all(is.na(memproj[R,C,rcp45idx])==FALSE)==TRUE){
          sdmat[R,C] = sd(memproj[R,C,rcp45idx],na.rm=TRUE)
        }
      }
    }
    
    upperbd = ensmean+1*sdmat
    lowerbd = ensmean-1*sdmat
    
    tmpmid = memproj[,,rcp45idx]
    
    for(i in 1:27){
      tmpmid[,,i] = ifelse(memproj[,,rcp45idx[i]]<=upperbd & memproj[,,rcp45idx[i]]>= lowerbd,1,0) 
    }
  }

if(typepart=="thres"){
  
  upperbd = ensmean+threshold
  lowerbd = ensmean-threshold
  
  tmpmid = memproj[,,rcp45idx]
  
  for(i in 1:27){
    tmpmid[,,i] = ifelse(memproj[,,rcp45idx[i]]<=upperbd & memproj[,,rcp45idx[i]]>= lowerbd,1,0) 
  }
}  

  
  agreemid45 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
}

  if(type==2){
    ensmean=meanproj45
    if(typepart=="chgper"){
    quants = quantile(memproj[,,rcp45idx],probs=seq(0.01,0.99,by=0.01),na.rm=TRUE)
    quantdiff = abs(quants-mean(ensmean,na.rm=TRUE))
    minidx=which(quantdiff==min(quantdiff))
    lwbd = minidx-percbound
    if(lwbd<1) lwbd=1
    upbd = minidx+percbound
    if(upbd>99) upbd=99
    
    quantrange = c(quantdiff[lwbd],quantdiff[upbd])
    
    upperbd = ensmean+quantrange[2]
    lowerbd = ensmean-quantrange[1]
    
    tmpmid = memproj[,,rcp45idx]
    for(i in 1:27){
      tmpmid[,,i] = ifelse(memproj[,,rcp45idx[i]]<=upperbd & memproj[,,rcp45idx[i]]>=lowerbd,1,0)  
    }
    }
    
    if(typepart=="sdbound"){
      sdval = sd(memproj[,,rcp45idx],na.rm=TRUE)
      upperbd = ensmean+1*sdval
      lowerbd = ensmean-1*sdval
      
      tmpmid = memproj[,,rcp45idx]
      for(i in 1:27){
        tmpmid[,,i] = ifelse(memproj[,,rcp45idx[i]]<=upperbd & memproj[,,rcp45idx[i]]>=lowerbd,1,0)  
      }
    }
    agreemid45 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
  }

  # RCP8.5 close to ens mean
if(type==1){  
ensmean = apply(memproj[,,rcp85idx],c(1,2),mean,na.rm=TRUE)
if(typepart=="chgper"){  
chgper = memproj[,,rcp85idx]
  for(i in 1:27){
    top = memproj[,,rcp85idx[i]]-ensmean
    bottom = abs(ensmean)
    chgper[,,i] = (top/bottom)*100
  }
  
  tmppos = ifelse(chgper>pb[2],1,0)
  tmpmid = ifelse(chgper<=pb[2] & chgper>= pb[1],1,0)
  tmpneg = ifelse(chgper< pb[1],1,0)
}
  
  if(typepart=="sdbound"){
    sdmat = matrix(NA,nrow=length(lon),ncol=length(lat))
    for(R in 1:length(lon)){
      for(C in 1:length(lat)){
        if(all(is.na(memproj[R,C,rcp85idx])==FALSE)==TRUE){
          sdmat[R,C] = sd(memproj[R,C,rcp85idx],na.rm=TRUE)
        }
      }
    }
    
    upperbd = ensmean+1*sdmat
    lowerbd = ensmean-1*sdmat
    
    tmpmid = memproj[,,rcp85idx]
    
    for(i in 1:27){
      tmpmid[,,i] = ifelse(memproj[,,rcp85idx[i]]<=upperbd & memproj[,,rcp85idx[i]]>= lowerbd,1,0) 
    }
  }

if(typepart=="thres"){
  
  upperbd = ensmean+threshold
  lowerbd = ensmean-threshold
  
  tmpmid = memproj[,,rcp85idx]
  
  for(i in 1:27){
    tmpmid[,,i] = ifelse(memproj[,,rcp85idx[i]]<=upperbd & memproj[,,rcp85idx[i]]>= lowerbd,1,0) 
  }
}  

  agreemid85 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
}

if(type==2){
  ensmean=meanproj85
  if(typepart=="chgper"){
  quants = quantile(memproj[,,rcp85idx],probs=seq(0.01,0.99,by=0.01),na.rm=TRUE)
  quantdiff = abs(quants-mean(ensmean,na.rm=TRUE))
  minidx=which(quantdiff==min(quantdiff))
  lwbd = minidx-percbound
  if(lwbd<1) lwbd=1
  upbd = minidx+percbound
  if(upbd>99) upbd=99
  
  quantrange = c(quantdiff[lwbd],quantdiff[upbd])
  
  upperbd = ensmean+quantrange[2]
  lowerbd = ensmean-quantrange[1]
  
  tmpmid = memproj[,,rcp85idx]
  for(i in 1:27){
    tmpmid[,,i] = ifelse(memproj[,,rcp85idx[i]]<=upperbd & memproj[,,rcp85idx[i]]>=lowerbd,1,0)  
  }
  }
  
  if(typepart=="sdbound"){
    sdval = sd(memproj[,,rcp85idx],na.rm=TRUE)
    upperbd = ensmean+1*sdval
    lowerbd = ensmean-1*sdval
    
    tmpmid = memproj[,,rcp85idx]
    for(i in 1:27){
      tmpmid[,,i] = ifelse(memproj[,,rcp85idx[i]]<=upperbd & memproj[,,rcp85idx[i]]>=lowerbd,1,0)  
    }
  }
  
  agreemid85 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
}
  
agree85v2 = ifelse(agreemid85>=ensthres,1,0)
agree45v2 = ifelse(agreemid45>=ensthres,1,0)
agree26v2 = ifelse(agreemid26>=ensthres,1,0)

######
# plot maps

testsfc = list(x=lon-360,y=lat,z=agree85v2)
surface(testsfc,type="I")
map("state",add=TRUE)

###

agree85idx = which(agree85v2==1,arr.ind=TRUE)
agree45idx = which(agree45v2==1,arr.ind=TRUE)
agree26idx = which(agree26v2==1,arr.ind=TRUE)

diffcolorbar = colorramp(c(meanproj85,meanproj45,meanproj26),colorchoice=colorchoice,Blimit=BINS,use_fixed_scale = USEFS,fixed_scale=colorlimits)

if(typepart == "chgper"){
filename = paste("/home/woot0002/stippletest_",varname,"_Perc",percbound,"_Ensthres",ensthres,"_V",type,".pdf",sep="")
}

if(typepart == "sdbound"){
  filename = paste("/home/woot0002/stippletest_",varname,"_1SD_Ensthres",ensthres,"_V",type,".pdf",sep="")
}

if(typepart == "thres"){
  filename = paste("/home/woot0002/stippletest_",varname,"_thres",threshold,"_Ensthres",ensthres,"_V",type,".pdf",sep="")
}


pdf(filename,height=6,width=6,onefile=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=meanproj26)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP2.6",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree26idx[which(1:nrow(agree26idx)%%2==0),1]]-360, lat[agree26idx[which(1:nrow(agree26idx)%%2==0),2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4) #alpha=0.4

testsfc1 = list(x=lon-360,y=lat,z=meanproj45)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree45idx[which(1:nrow(agree45idx)%%2==0),1]]-360, lat[agree45idx[which(1:nrow(agree45idx)%%2==0),2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

testsfc1 = list(x=lon-360,y=lat,z=meanproj85)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree85idx[which(1:nrow(agree85idx)%%2==0),1]]-360, lat[agree85idx[which(1:nrow(agree85idx)%%2==0),2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

dev.off()


