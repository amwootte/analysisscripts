
library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("/data2/3to5/I35/scripts/colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maptools)
library(ggplot2)
library(zoo)
library(lars)

streamflowID = "boisdearc"
loclon = -98.90301639
loclat = 31.19016833
futureperiod = c(2006,2099)

#####
# C-PrEP data

tasmaxfile_hist = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_*historical*.nc",intern=TRUE))

tasminfile_hist = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_*historical*.nc",intern=TRUE))

tasmaxfile_proj = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_day*rcp*.nc",intern=TRUE))

tasminfile_proj = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_day*rcp*.nc",intern=TRUE))


filebreakdown = do.call(rbind,strsplit(tasmaxfile_proj,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(tasmaxfile_hist,"_",fixed=TRUE))
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

#####
# determine which grid cell to pull

nctest = nc_open(tasmaxfile_hist[1])
lon = ncvar_get(nctest,"lon")-360
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

###
# create model grid

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

###
# get cells to use

pointarea = distfunc(loclon,loclat,modelgrid)

###
# pull tmax and tmin data

histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
climatehist = NULL
for(i in 1:length(tasmaxfile_hist)){
  
  ###
  # Get TMAX, and TMIN data
  
  nctest = nc_open(tasmaxfile_hist[i])
  TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  nctest = nc_open(tasminfile_hist[i])
  TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  if(length(TMAX)<length(histdates) | length(TMIN)<length(histdates)){
    idxout = which(substr(histdates,6,10)=="02-29")
    datesin = histdates[-idxout]

    if(length(TMAX)==length(histdates)){
      TMAX = TMAX[-idxout]
    }
    if(length(TMIN)==length(histdates)){
      TMIN = TMIN[-idxout]
    }
  } else {
    datesin = histdates
  }
  
  #####
  # unit conversion
  
  TMAX = TMAX-273.15
  TMIN = TMIN-273.15
  DATE = datesin
  modinfo = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i])
  tmpdat = data.frame(modinfo,DATE,TMAX,TMIN)
  climatehist = rbind(climatehist,tmpdat)
  message("Finished gathering results for file ",i," / ",length(tasmaxfile_hist))
}

###
# pull tmax and tmin data

histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
climatehist = NULL
for(i in 1:length(tasmaxfile_hist)){
  
  ###
  # Get TMAX, and TMIN data
  
  nctest = nc_open(tasmaxfile_hist[i])
  TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  nctest = nc_open(tasminfile_hist[i])
  TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  if(length(TMAX)<length(histdates) | length(TMIN)<length(histdates)){
    idxout = which(substr(histdates,6,10)=="02-29")
    datesin = histdates[-idxout]
    
    if(length(TMAX)==length(histdates)){
      TMAX = TMAX[-idxout]
    }
    if(length(TMIN)==length(histdates)){
      TMIN = TMIN[-idxout]
    }
  } else {
    datesin = histdates
  }
  
  #####
  # unit conversion
  
  TMAX = TMAX-273.15
  TMIN = TMIN-273.15
  DATE = datesin
  modinfo = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")
  tmpdat = data.frame(modinfo,DATE,TMAX,TMIN)
  climatehist = rbind(climatehist,tmpdat)
  message("Finished gathering results for file ",i," / ",length(tasmaxfile_hist))
}

###

###
# pull tmax and tmin data - projected

projdates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
climateproj = NULL
for(i in 1:length(tasmaxfile_proj)){
  
  ###
  # Get TMAX, and TMIN data
  
  nctest = nc_open(tasmaxfile_proj[i])
  TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  nctest = nc_open(tasminfile_proj[i])
  TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  nc_close(nctest)
  
  if(length(TMAX)<length(projdates) | length(TMIN)<length(projdates)){
    idxout = which(substr(projdates,6,10)=="02-29")
    datesin = projdates[-idxout]
    
    if(length(TMAX)==length(projdates)){
      TMAX = TMAX[-idxout]
    }
    if(length(TMIN)==length(projdates)){
      TMIN = TMIN[-idxout]
    }
  } else {
    datesin = projdates
  }
  
  #####
  # unit conversion
  
  TMAX = TMAX-273.15
  TMIN = TMIN-273.15
  DATE = datesin
  modinfo = paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_")
  tmpdat = data.frame(modinfo,DATE,TMAX,TMIN)
  climateproj = rbind(climateproj,tmpdat)
  message("Finished gathering results for file ",i," / ",length(tasmaxfile_proj))
}

#####
# create 5day smooth TMAX and TMIN

TMAX5D = TMIN5D = c()
mods = unique(climatehist$modinfo)

for(i in 1:length(mods)){
  TMAX5D = c(TMAX5D,filter(climatehist$TMAX[which(climatehist$modinfo==mods[i])],rep(1/5,5),sides=2))  
  TMIN5D = c(TMIN5D,filter(climatehist$TMIN[which(climatehist$modinfo==mods[i])],rep(1/5,5),sides=2))
}

climatehist$TMAX5D = TMAX5D
climatehist$TMIN5D = TMIN5D

###

TMAX5D = TMIN5D = c()
mods = unique(climateproj$modinfo)

for(i in 1:length(mods)){
  TMAX5D = c(TMAX5D,filter(climateproj$TMAX[which(climateproj$modinfo==mods[i])],rep(1/5,5),sides=2))  
  TMIN5D = c(TMIN5D,filter(climateproj$TMIN[which(climateproj$modinfo==mods[i])],rep(1/5,5),sides=2))
}

climateproj$TMAX5D = TMAX5D
climateproj$TMIN5D = TMIN5D

#####
# convert air to stream temperature - TMAX

load(paste("TMAX5daylmfit_",streamflowID,".Rdata",sep=""))
climatehist$tmaxfit5d = coefficients(lmfit)[2]*climatehist$TMAX5D+coefficients(lmfit)[1]    #predict(lmfit,newdata=climatehist$TMAX5D)
climateproj$tmaxfit5d = coefficients(lmfit)[2]*climateproj$TMAX5D+coefficients(lmfit)[1]

#####
# convert air to stream temperature - TMIN

load(paste("TMIN5daylmfit_",streamflowID,".Rdata",sep=""))
climatehist$tminfit5d = coefficients(lmfit)[2]*climatehist$TMIN5D+coefficients(lmfit)[1]     #predict(lmfit,newdata=climatehist$TMAX5D)
climateproj$tminfit5d = coefficients(lmfit)[2]*climateproj$TMIN5D+coefficients(lmfit)[1]

######
# climatologies - historical

climatehist$tavgfit5d = (climatehist$tmaxfit5d+climatehist$tminfit5d)/2
climatehist$tavgfit5d4d = filter(climatehist$tavgfit5d,rep(1/4,4),sides=2)
histclimo = aggregate(climatehist[,3:ncol(climatehist)],by=list(modinfo=climatehist$modinfo),mean,na.rm=TRUE)

######
# get quantile values and occurence climo.

LT05_Gperc = 0.7864
LT50_Gperc = 0.9999
LT05_Jperc = 0.9891
LT50_Jperc = 0.9999

histquants_LT05_Gperc = aggregate(climatehist$tavgfit5d,by=list(modinfo=climatehist$modinfo),quantile,probs=LT05_Gperc,na.rm=TRUE)
histquants_LT50_Gperc = aggregate(climatehist$tavgfit5d,by=list(modinfo=climatehist$modinfo),quantile,probs=LT50_Gperc,na.rm=TRUE)

histquants_LT05_Lperc = histquants_LT05_Gperc
histquants_LT50_Lperc = histquants_LT05_Gperc

modinfo = as.character(unique(climatehist$modinfo))
for(i in 1:length(modinfo)){
  idxin = which(as.character(climatehist$modinfo)==modinfo[i])
  modidxin = which(as.character(histquants_LT05_Lperc[,1])==modinfo[i])
  histquants_LT05_Lperc[modidxin,2]=quantile(climatehist$tavgfit5d4d[idxin],probs=LT05_Jperc,na.rm=TRUE)
  modidxin = which(as.character(histquants_LT50_Lperc[,1])==modinfo[i])
  histquants_LT50_Lperc[modidxin,2]=quantile(climatehist$tavgfit5d4d[idxin],probs=LT50_Jperc,na.rm=TRUE)
}

years = as.numeric(substr(climatehist$DATE,1,4))

climatehist$season = NA
DJFidx = which(as.numeric(substr(climatehist$DATE,6,7))>=12 | as.numeric(substr(climatehist$DATE,6,7))<=2)
MAMidx = which(as.numeric(substr(climatehist$DATE,6,7))>=3 & as.numeric(substr(climatehist$DATE,6,7))<=5)
JJAidx = which(as.numeric(substr(climatehist$DATE,6,7))>=6 & as.numeric(substr(climatehist$DATE,6,7))<=8)
SONidx = which(as.numeric(substr(climatehist$DATE,6,7))>=9 & as.numeric(substr(climatehist$DATE,6,7))<=11)

climatehist$season[DJFidx]=1
climatehist$season[MAMidx]=2
climatehist$season[JJAidx]=3
climatehist$season[SONidx]=4


LT05gclimo = LT50gclimo = LT05jclimo = LT50jclimo =  c()
LT05g2dspclimo = LT50g2dspclimo = LT05j2dspclimo = LT50j2dspclimo =  c()


for(i in 1:length(modinfo)){
  idxin = which(as.character(climatehist$modinfo)==modinfo[i])
  tmp = climatehist[idxin,]
  
  modidxin = which(as.character(histquants_LT05_Gperc[,1])==modinfo[i])
  tmp$LT05g = ifelse(tmp$tavgfit5d>=histquants_LT05_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Gperc[,1])==modinfo[i])
  tmp$LT50g = ifelse(tmp$tavgfit5d>=histquants_LT50_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT05_Lperc[,1])==modinfo[i])
  tmp$LT05j = ifelse(tmp$tavgfit5d4d>=histquants_LT05_Lperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Lperc[,1])==modinfo[i])
  tmp$LT50j = ifelse(tmp$tavgfit5d4d>=histquants_LT50_Lperc[modidxin,2],1,0)
  
  yearlycount = aggregate(tmp[,12:15],by=list(year=as.numeric(substr(tmp$DATE,1,4))),sum,na.rm=TRUE)
  
  tmp$year = as.numeric(substr(tmp$DATE,1,4))
  years = unique(tmp$year)
  
  LT05g_2dsp = LT50g_2dsp = LT05j_2dsp = LT50j_2dsp = c()
  
  for(y in 1:length(years)){
    
    yearidx = which(tmp$year==years[y])
    tmpin = tmp$LT05g[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT05g_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT50g[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT50g_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT05j[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT05j_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT50j[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT50j_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
  }
  
  LT05g2dspclimo[i] = mean(LT05g_2dsp)
  LT50g2dspclimo[i] = mean(LT50g_2dsp)
  LT05j2dspclimo[i] = mean(LT05j_2dsp)
  LT50j2dspclimo[i] = mean(LT50j_2dsp)
  
  LT05gclimo[i] = mean(yearlycount$LT05g)
  LT50gclimo[i] = mean(yearlycount$LT50g)
  LT05jclimo[i] = mean(yearlycount$LT05j)
  LT50jclimo[i] = mean(yearlycount$LT50j)
  
}

LThistclimos = data.frame(modinfo,LT05gclimo,LT50gclimo,LT05jclimo,LT50jclimo,LT05g2dspclimo,LT50g2dspclimo,LT05j2dspclimo,LT50j2dspclimo)

###
# seasonal climos

seasonclimo = NULL
for(i in 1:length(modinfo)){
  idxin = which(as.character(climatehist$modinfo)==modinfo[i])
  tmp = climatehist[idxin,]
  
  modidxin = which(as.character(histquants_LT05_Gperc[,1])==modinfo[i])
  tmp$LT05g = ifelse(tmp$tavgfit5d>=histquants_LT05_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Gperc[,1])==modinfo[i])
  tmp$LT50g = ifelse(tmp$tavgfit5d>=histquants_LT50_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT05_Lperc[,1])==modinfo[i])
  tmp$LT05j = ifelse(tmp$tavgfit5d4d>=histquants_LT05_Lperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Lperc[,1])==modinfo[i])
  tmp$LT50j = ifelse(tmp$tavgfit5d4d>=histquants_LT50_Lperc[modidxin,2],1,0)
  
  yearlycount = aggregate(tmp[,12:15],by=list(season=paste(as.numeric(substr(tmp$DATE,1,4)),tmp$season,sep="-")),sum,na.rm=TRUE)
  
  seasontmp = aggregate(yearlycount[,2:5],by=list(season=as.numeric(substr(yearlycount[,1],6,6))),mean,na.rm=TRUE)
  seasontmp$modinfo = modinfo[i]
  seasonclimo = rbind(seasonclimo,seasontmp)
}

LThistseasclimos = seasonclimo

######
# climatologies - projected

climateproj_mid = subset(climateproj,as.numeric(substr(DATE,1,4))>=2036 & as.numeric(substr(DATE,1,4))<=2065)
climateproj_end = subset(climateproj,as.numeric(substr(DATE,1,4))>=2070 & as.numeric(substr(DATE,1,4))<=2099)

climateproj_mid$tavgfit5d = (climateproj_mid$tmaxfit5d+climateproj_mid$tminfit5d)/2
climateproj_end$tavgfit5d = (climateproj_end$tmaxfit5d+climateproj_end$tminfit5d)/2

climateproj_mid$tavgfit5d4d = filter(climateproj_mid$tavgfit5d,rep(1/4,4),sides=2)
climateproj_end$tavgfit5d4d = filter(climateproj_end$tavgfit5d,rep(1/4,4),sides=2)

projclimo_mid = aggregate(climateproj_mid[,3:ncol(climateproj_mid)],by=list(modinfo=climateproj_mid$modinfo),mean,na.rm=TRUE)
projclimo_end = aggregate(climateproj_end[,3:ncol(climateproj_end)],by=list(modinfo=climateproj_end$modinfo),mean,na.rm=TRUE)

climateproj_mid$season = NA
DJFidx = which(as.numeric(substr(climateproj_mid$DATE,6,7))>=12 | as.numeric(substr(climateproj_mid$DATE,6,7))<=2)
MAMidx = which(as.numeric(substr(climateproj_mid$DATE,6,7))>=3 & as.numeric(substr(climateproj_mid$DATE,6,7))<=5)
JJAidx = which(as.numeric(substr(climateproj_mid$DATE,6,7))>=6 & as.numeric(substr(climateproj_mid$DATE,6,7))<=8)
SONidx = which(as.numeric(substr(climateproj_mid$DATE,6,7))>=9 & as.numeric(substr(climateproj_mid$DATE,6,7))<=11)

climateproj_mid$season[DJFidx]=1
climateproj_mid$season[MAMidx]=2
climateproj_mid$season[JJAidx]=3
climateproj_mid$season[SONidx]=4

climateproj_end$season = NA
DJFidx = which(as.numeric(substr(climateproj_end$DATE,6,7))>=12 | as.numeric(substr(climateproj_end$DATE,6,7))<=2)
MAMidx = which(as.numeric(substr(climateproj_end$DATE,6,7))>=3 & as.numeric(substr(climateproj_end$DATE,6,7))<=5)
JJAidx = which(as.numeric(substr(climateproj_end$DATE,6,7))>=6 & as.numeric(substr(climateproj_end$DATE,6,7))<=8)
SONidx = which(as.numeric(substr(climateproj_end$DATE,6,7))>=9 & as.numeric(substr(climateproj_end$DATE,6,7))<=11)

climateproj_end$season[DJFidx]=1
climateproj_end$season[MAMidx]=2
climateproj_end$season[JJAidx]=3
climateproj_end$season[SONidx]=4


#######
# threshold climos - mid-century

modinfo = as.character(unique(climateproj_mid$modinfo))
LT05gclimo = LT50gclimo = LT05jclimo = LT50jclimo =  c()
LT05g2dspclimo = LT50g2dspclimo = LT05j2dspclimo = LT50j2dspclimo =  c()

for(i in 1:length(modinfo)){
  idxin = which(as.character(climateproj_mid$modinfo)==modinfo[i])
  tmp = climateproj_mid[idxin,]
  
  modidxin = which(as.character(histquants_LT05_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05g = ifelse(tmp$tavgfit5d>=histquants_LT05_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50g = ifelse(tmp$tavgfit5d>=histquants_LT50_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT05_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05j = ifelse(tmp$tavgfit5d4d>=histquants_LT05_Lperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50j = ifelse(tmp$tavgfit5d4d>=histquants_LT50_Lperc[modidxin,2],1,0)
  
  yearlycount = aggregate(tmp[,12:15],by=list(year=as.numeric(substr(tmp$DATE,1,4))),sum,na.rm=TRUE)
  
  tmp$year = as.numeric(substr(tmp$DATE,1,4))
  years = unique(tmp$year)
  
  LT05g_2dsp = LT50g_2dsp = LT05j_2dsp = LT50j_2dsp = c()
  
  for(y in 1:length(years)){
    
    yearidx = which(tmp$year==years[y])
    tmpin = tmp$LT05g[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT05g_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT50g[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT50g_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT05j[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT05j_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT50j[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT50j_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
  }
  
  LT05g2dspclimo[i] = mean(LT05g_2dsp)
  LT50g2dspclimo[i] = mean(LT50g_2dsp)
  LT05j2dspclimo[i] = mean(LT05j_2dsp)
  LT50j2dspclimo[i] = mean(LT50j_2dsp)
  
  LT05gclimo[i] = mean(yearlycount$LT05g)
  LT50gclimo[i] = mean(yearlycount$LT50g)
  LT05jclimo[i] = mean(yearlycount$LT05j)
  LT50jclimo[i] = mean(yearlycount$LT50j)
  
}

LTprojmidclimos = data.frame(modinfo,LT05gclimo,LT50gclimo,LT05jclimo,LT50jclimo,LT05g2dspclimo,LT50g2dspclimo,LT05j2dspclimo,LT50j2dspclimo)

###
#
modinfo = as.character(unique(climateproj_mid$modinfo))
seasonclimo = NULL
for(i in 1:length(modinfo)){
  idxin = which(as.character(climateproj_mid$modinfo)==modinfo[i])
  tmp =climateproj_mid[idxin,]
  
  modidxin = which(as.character(histquants_LT05_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05g = ifelse(tmp$tavgfit5d>=histquants_LT05_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50g = ifelse(tmp$tavgfit5d>=histquants_LT50_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT05_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05j = ifelse(tmp$tavgfit5d4d>=histquants_LT05_Lperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50j = ifelse(tmp$tavgfit5d4d>=histquants_LT50_Lperc[modidxin,2],1,0)
  
  yearlycount = aggregate(tmp[,12:15],by=list(season=paste(as.numeric(substr(tmp$DATE,1,4)),tmp$season,sep="-")),sum,na.rm=TRUE)
  
  seasontmp = aggregate(yearlycount[,2:5],by=list(season=as.numeric(substr(yearlycount[,1],6,6))),mean,na.rm=TRUE)
  seasontmp$modinfo = modinfo[i]
  seasonclimo = rbind(seasonclimo,seasontmp)
}

LTprojmidseasclimos = seasonclimo

#######
# end-century

modinfo = as.character(unique(climateproj_end$modinfo))
LT05gclimo = LT50gclimo = LT05jclimo = LT50jclimo =  c()
LT05g2dspclimo = LT50g2dspclimo = LT05j2dspclimo = LT50j2dspclimo =  c()

for(i in 1:length(modinfo)){
  idxin = which(as.character(climateproj_end$modinfo)==modinfo[i])
  tmp = climateproj_end[idxin,]
  
  modidxin = which(as.character(histquants_LT05_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05g = ifelse(tmp$tavgfit5d>=histquants_LT05_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50g = ifelse(tmp$tavgfit5d>=histquants_LT50_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT05_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05j = ifelse(tmp$tavgfit5d4d>=histquants_LT05_Lperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50j = ifelse(tmp$tavgfit5d4d>=histquants_LT50_Lperc[modidxin,2],1,0)
  
  yearlycount = aggregate(tmp[,12:15],by=list(year=as.numeric(substr(tmp$DATE,1,4))),sum,na.rm=TRUE)
  
  tmp$year = as.numeric(substr(tmp$DATE,1,4))
  years = unique(tmp$year)
  
  LT05g_2dsp = LT50g_2dsp = LT05j_2dsp = LT50j_2dsp = c()
  
  for(y in 1:length(years)){
    
    yearidx = which(tmp$year==years[y])
    tmpin = tmp$LT05g[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT05g_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT50g[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT50g_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT05j[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT05j_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
    tmpin = tmp$LT50j[yearidx]
    tmpin[which(is.na(tmpin)==TRUE)]=0
    LT50j_2dsp[y]=spell_length_calc(tmpin,premasked=TRUE,spell_len=2,outtype="count",thold=1)
  }
  
  LT05g2dspclimo[i] = mean(LT05g_2dsp)
  LT50g2dspclimo[i] = mean(LT50g_2dsp)
  LT05j2dspclimo[i] = mean(LT05j_2dsp)
  LT50j2dspclimo[i] = mean(LT50j_2dsp)
  
  LT05gclimo[i] = mean(yearlycount$LT05g)
  LT50gclimo[i] = mean(yearlycount$LT50g)
  LT05jclimo[i] = mean(yearlycount$LT05j)
  LT50jclimo[i] = mean(yearlycount$LT50j)
  
}

LTprojendclimos = data.frame(modinfo,LT05gclimo,LT50gclimo,LT05jclimo,LT50jclimo,LT05g2dspclimo,LT50g2dspclimo,LT05j2dspclimo,LT50j2dspclimo)

###
modinfo = as.character(unique(climateproj_end$modinfo))
seasonclimo = NULL
for(i in 1:length(modinfo)){
  idxin = which(as.character(climateproj_end$modinfo)==modinfo[i])
  tmp =climateproj_end[idxin,]
  
  modidxin = which(as.character(histquants_LT05_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05g = ifelse(tmp$tavgfit5d>=histquants_LT05_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Gperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50g = ifelse(tmp$tavgfit5d>=histquants_LT50_Gperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT05_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT05j = ifelse(tmp$tavgfit5d4d>=histquants_LT05_Lperc[modidxin,2],1,0)
  
  modidxin = which(as.character(histquants_LT50_Lperc[,1])==substr(modinfo[i],7,nchar(modinfo[i])))
  tmp$LT50j = ifelse(tmp$tavgfit5d4d>=histquants_LT50_Lperc[modidxin,2],1,0)
  
  yearlycount = aggregate(tmp[,12:15],by=list(season=paste(as.numeric(substr(tmp$DATE,1,4)),tmp$season,sep="-")),sum,na.rm=TRUE)
  
  seasontmp = aggregate(yearlycount[,2:5],by=list(season=as.numeric(substr(yearlycount[,1],6,6))),mean,na.rm=TRUE)
  seasontmp$modinfo = modinfo[i]
  seasonclimo = rbind(seasonclimo,seasontmp)
}

LTprojendseasclimos = seasonclimo

#######
# combine climos

histclimo = merge(histclimo,LThistclimos,by="modinfo")
projclimo_mid = merge(projclimo_mid,LTprojmidclimos,by="modinfo")
projclimo_end = merge(projclimo_end,LTprojendclimos,by="modinfo")

#######
# projected changes

projclimoanom_mid = projclimo_mid
projclimoanom_end = projclimo_end

for(i in 1:nrow(projclimo_mid)){
  
  histidx = grep(substr(as.character(projclimo_mid$modinfo[i]),7,nchar(as.character(projclimo_mid$modinfo[i]))),as.character(histclimo$modinfo))
  projclimoanom_mid[i,2:ncol(projclimoanom_mid)] = projclimo_mid[i,2:ncol(projclimo_mid)]-histclimo[histidx,2:ncol(histclimo)]
  projclimoanom_end[i,2:ncol(projclimoanom_end)] = projclimo_end[i,2:ncol(projclimo_end)]-histclimo[histidx,2:ncol(histclimo)]
  
}

###
# threshold seasonal changes

LTprojclimoanom_mid = LTprojmidseasclimos
LTprojclimoanom_end = LTprojendseasclimos

mods = unique(LTprojendseasclimos$modinfo)
season = 1:4

for(i in 1:nrow(LTprojmidseasclimos)){
  
  modidx = which(LTprojmidseasclimos$modinfo == mods[i])
  histidx = which(LThistseasclimos$modinfo==substr(as.character(mods[i]),7,nchar(as.character(mods[i]))))
  LTprojclimoanom_mid[modidx,2:5] = LTprojmidseasclimos[modidx,2:5]-LThistseasclimos[histidx,2:5]
  LTprojclimoanom_end[modidx,2:5] = LTprojendseasclimos[modidx,2:5]-LThistseasclimos[histidx,2:5]

}

LTprojclimoanom_mid$period = "mid-century"
LTprojclimoanom_end$period = "end-of-century"

######
# plotting results

projclimoanom_mid$scen = substr(as.character(projclimoanom_mid$modinfo),1,5)
projclimoanom_end$scen = substr(as.character(projclimoanom_end$modinfo),1,5)

projclimoanom_mid$period = "mid-century"
projclimoanom_end$period = "end-of-century"

projclimoanom = rbind(projclimoanom_mid,projclimoanom_end)

projclimoanom$period = factor(projclimoanom$period,levels=c("mid-century","end-of-century"))

LTseasclimoanom = rbind(LTprojclimoanom_mid,LTprojclimoanom_end)
LTseasclimoanom$scen = substr(as.character(LTseasclimoanom$modinfo),1,5)
LTseasclimoanom$period = factor(LTseasclimoanom$period,levels=c("mid-century","end-of-century"))

projclimoanom$streamflowID = streamflowID
save(list=c("projclimoanom"),file=paste("/home/woot0002/CPREP_streamtemp_",streamflowID,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))

load(file=paste("/home/woot0002/CPREP_streamtemp_",streamflowID,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))

projclimoanom_end = subset(projclimoanom,period=="end-of-century")

library(ggplot2)

ggplot(projclimoanom_end, aes(x=scen, y=tmaxfit5d,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change High Stream Temperatures ",streamflowID,sep=""))+xlab("")+ylab("Degrees C")+ylim(0,3)
ggsave(paste("/home/woot0002/streamtempchange_high_",streamflowID,".pdf",sep=""),device = "pdf",width=5,height=5)

ggplot(projclimoanom_end, aes(x=scen, y=tminfit5d,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Low Stream Temperatures ",streamflowID,sep=""))+xlab("")+ylab("Degrees C")+ylim(0,3)
ggsave(paste("/home/woot0002/streamtempchange_low_",streamflowID,".pdf",sep=""),device = "pdf",width=5,height=5)

ggplot(projclimoanom_end, aes(x=scen, y=tavgfit5d,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Average Stream Temperatures ",streamflowID,sep=""))+xlab("")+ylab("Degrees C")+ylim(0,3)
ggsave(paste("/home/woot0002/streamtempchange_avg_",streamflowID,".pdf",sep=""),device = "pdf",width=5,height=5)

ggplot(projclimoanom_end, aes(x=scen, y=LT05gclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 24h LT05 Glochidia Threshold ",streamflowID,sep=""))+xlab("")+ylab("days/year")+ylim(0,120)
ggsave(paste("/home/woot0002/streamtempchange_LT05g_",streamflowID,".pdf",sep=""),device = "pdf",width=5,height=5)

ggplot(projclimoanom_end, aes(x=scen, y=LT50gclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 24h LT50 Glochidia Threshold ",streamflowID,sep=""))+xlab("")+ylab("days/year")+ylim(0,120)
ggsave(paste("/home/woot0002/streamtempchange_LT50g_",streamflowID,".pdf",sep=""),device = "pdf",width=5,height=5)

ggplot(projclimoanom_end, aes(x=scen, y=LT05jclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 96h LT05 Juvenile Threshold ",streamflowID,sep=""))+xlab("")+ylab("events/year")+ylim(0,120)
ggsave(paste("/home/woot0002/streamtempchange_LT05j_",streamflowID,".pdf",sep=""),device = "pdf",width=5,height=5)

ggplot(projclimoanom_end, aes(x=scen, y=LT50jclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 96h LT50 Juvenile Threshold ",streamflowID,sep=""))+xlab("")+ylab("events/year")+ylim(0,120)
ggsave(paste("/home/woot0002/streamtempchange_LT50j_",streamflowID,".pdf",sep=""),device = "pdf",width=5,height=5)

aggregate(projclimoanom_end[,2:13],by=list(scen=projclimoanom_end$scen),mean,na.rm=TRUE)
aggregate(projclimoanom_end[,2:13],by=list(scen=projclimoanom_end$scen),sd,na.rm=TRUE)

####
# spells

ggplot(projclimoanom, aes(x=scen, y=LT05g2dspclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of >=4 day spells of \n 24h LT05 Glochidia Threshold ",streamflowID,sep=""))+xlab("")+ylab("events")+facet_wrap(vars(period))

ggplot(projclimoanom, aes(x=scen, y=LT50g2dspclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of >=4 day spells of \n 24h LT50 Glochidia Threshold ",streamflowID,sep=""))+xlab("")+ylab("events")+facet_wrap(vars(period))

ggplot(projclimoanom, aes(x=scen, y=LT05j2dspclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of >=4 day spells of \n 96h LT05 Juvenile Threshold ",streamflowID,sep=""))+xlab("")+ylab("events")+facet_wrap(vars(period))

ggplot(projclimoanom, aes(x=scen, y=LT50j2dspclimo,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of >=4 day spells of \n 96h LT50 Juvenile Threshold ",streamflowID,sep=""))+xlab("")+ylab("events")+facet_wrap(vars(period))

####

ggplot(LTseasclimoanom, aes(x=scen, y=LT05g,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 24h LT05 Glochidia Threshold ",streamflowID,sep=""))+xlab("")+ylab("days/year")+facet_grid(cols =vars(season),rows=vars(period))

ggplot(LTseasclimoanom, aes(x=scen, y=LT50g,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 24h LT50 Glochidia Threshold ",streamflowID,sep=""))+xlab("")+ylab("days/year")+facet_grid(cols =vars(season),rows=vars(period))

ggplot(LTseasclimoanom, aes(x=scen, y=LT05j,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 96h LT05 Juvenile Threshold ",streamflowID,sep=""))+xlab("")+ylab("days/year")+facet_grid(cols =vars(season),rows=vars(period))

ggplot(LTseasclimoanom, aes(x=scen, y=LT50j,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change Frequency of 96h LT50 Juvenile Threshold ",streamflowID,sep=""))+xlab("")+ylab("days/year")+facet_grid(cols =vars(season),rows=vars(period))


####

LTprojmidseasclimos$scen = substr(LTprojmidseasclimos$modinfo,1,5)

ggplot(LTprojmidseasclimos, aes(x=season, y=LT05g,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("Projected Frequency of 24h LT05 Glochidia Threshold Charlies")+xlab("")+ylab("days/year")+facet_wrap(vars(season))

LThistseasclimos$season = factor(LThistseasclimos$season,levels=1:4,labels=c("DJF","MAM","JJA","SON"))

ggplot(LThistseasclimos, aes(x=season, y=LT05g)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("Historical 24h LT05 Glochidia Threshold Charlies")+xlab("")+ylab("days/year")




#+facet_grid(cols=vars(season),rows=vars(period))






