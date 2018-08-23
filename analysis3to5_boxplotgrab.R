#####################
# Gather data from 3^5 output for boxplots
# 

setwd("/home/woot0002/3to5/")
library(ncdf4)
library(maps)
library(maptools)

GCMs = c("CCSM4","MIROC5","MPI-ESM-LR")
GCMnum = 1:3
obscode = c("D","L","P")
exps = c("r6i1p1","r1i1p1","r1i1p1")
obs = c("daymet","livneh","prism")

#####
# Gather obs data first

obsprcptot = obsrnnmm = obsrx1day = list()

for(i in 1:length(obs)){
  OHfile = paste("pr_day_",obs[i],"_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc",sep="")
  dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  test = nc_open(OHfile)
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    times = ncvar_get(test,"time")
    
    if(length(times)<length(dates)){
      dates = dates[-which(substr(dates,6,10)=="02-29")]
    }
    years = unique(as.numeric(substr(dates,1,4)))
    
    array1 = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
    for(y in 1:length(years)){
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
      temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
      array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
      array1[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array1[,,y])
      array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
      array2[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array2[,,y])
      temprain = ifelse(temp>1,1,0)
      array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
      array3[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array3[,,y])
      message("Calcs for year ",years[y]," complete")
    }
    
  nc_close(test)
  
  obsprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  obsrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  obsrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for obs ",obs[i]," complete")
}

######
# Grab historical GCM

MHprcptot = MHrnnmm = MHrx1day = list()

for(i in 1:length(GCMs)){
  MHfile = paste("pr_day_",GCMs[i],"_historical_",exps[i],"_SCCSC0p1_19810101-20051231.nc",sep="")
  dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  test = nc_open(MHfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  
  if(length(times)<length(dates)){
    dates = dates[-which(substr(dates,6,10)=="02-29")]
  }
  years = 1981:2005
  
  array1 = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array1[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  
  MHprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  MHrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  MHrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}

######
# Grab future GCM

MFprcptot = MFrnnmm = MFrx1day = list()

for(i in 1:length(GCMs)){
  MFfile = paste("pr_day_",GCMs[i],"_rcp85_",exps[i],"_SCCSC0p1_20060101-20991231.nc",sep="")
  dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
  test = nc_open(MFfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  
  if(length(times)<length(dates)){
    dates = dates[-which(substr(dates,6,10)=="02-29")]
  }
  years = 2070:2099
  
  array1 = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array1[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  
  MFprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  MFrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  MFrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}

######
# Grab historical QDM

QDMHprcptot = QDMHrnnmm = QDMHrx1day = list()
QDMHorder = c()

for(i in 1:length(GCMs)){
  for(j in 1:length(obs)){
  # pr_day_I35prp1-QDM-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc
  # paste("/data2/3to5/I35/pr/EQDM/pr_day_I35prp1-QDM-A",GCMnum[i],"0",obscode[i],"01K00_historical_",exps[i],"_I35Land_19810101-20051231.nc)
  MHfile = paste("/data2/3to5/I35/pr/EDQM/pr_day_I35prp1-QDM-A",GCMnum[i],"0",obscode[j],"01K00_historical_",exps[i],"_I35Land_19810101-20051231.nc",sep="")
  dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  test = nc_open(MHfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  
  if(length(times)<length(dates)){
    dates = dates[-which(substr(dates,6,10)=="02-29")]
  }
  years = unique(as.numeric(substr(dates,1,4)))
  
  array1 = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array1[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  
  outcode = (i-1)*3+j
  QDMHprcptot[[outcode]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  QDMHrnnmm[[outcode]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  QDMHrx1day[[outcode]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  QDMHorder[outcode] = paste(GCMs[i],obs[j],sep="_")
  
  message("Calcs for obs ",obs[j]," complete")
  }
  message("Calcs for GCM ",GCMs[i]," complete")
}

#####
# Grab future QDM

QDMFprcptot = QDMFrnnmm = QDMFrx1day = list()
QDMForder = c()

for(i in 1:length(GCMs)){
  for(j in 1:length(obs)){
    # pr_day_I35prp1-QDM-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc
    # paste("/data2/3to5/I35/pr/EQDM/pr_day_I35prp1-QDM-A",GCMnum[i],"0",obscode[i],"01K00_historical_",exps[i],"_I35Land_19810101-20051231.nc)
    MHfile = paste("/data2/3to5/I35/pr/EDQM/pr_day_I35prp1-QDM-A",GCMnum[i],"8",obscode[j],"01K00_rcp85_",exps[i],"_I35Land_20060101-20991231.nc",sep="")
    dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
    test = nc_open(MHfile)
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    times = ncvar_get(test,"time")
    
    if(length(times)<length(dates)){
      dates = dates[-which(substr(dates,6,10)=="02-29")]
    }
    years = 2070:2099
    
    array1 = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
    for(y in 1:length(years)){
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
      temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
      array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
      array1[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array1[,,y])
      array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
      array2[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array2[,,y])
      temprain = ifelse(temp>1,1,0)
      array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
      array3[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array3[,,y])
      message("Calcs for year ",years[y]," complete")
    }
    
    nc_close(test)
    
    outcode = (i-1)*3+j
    QDMFprcptot[[outcode]] = apply(array1,c(1,2),mean,na.rm=TRUE)
    QDMFrnnmm[[outcode]] = apply(array3,c(1,2),mean,na.rm=TRUE)
    QDMFrx1day[[outcode]] = apply(array2,c(1,2),mean,na.rm=TRUE)
    QDMForder[outcode] = paste(GCMs[i],obs[j],sep="_")
    
    message("Calcs for obs ",obs[j]," complete")
  }
  message("Calcs for GCM ",GCMs[i]," complete")
}

######
# Grab future DeltaSD

DeltaSDFprcptot = DeltaSDFrnnmm = DeltaSDFrx1day = list()
DeltaSDForder = c()

for(i in 1:length(GCMs)){
  for(j in 1:length(obs)){
    # pr_day_I35prp1-QDM-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc
    # paste("/data2/3to5/I35/pr/EQDM/pr_day_I35prp1-QDM-A",GCMnum[i],"0",obscode[i],"01K00_historical_",exps[i],"_I35Land_19810101-20051231.nc)
    MHfile = paste("/data2/3to5/I35/pr/DeltaSD/pr_day_I35prp1-DeltaSD-A",GCMnum[i],"8",obscode[j],"01K00_rcp85_",exps[i],"_I35Land_20060101-20991231.nc",sep="")
    dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
    test = nc_open(MHfile)
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    times = ncvar_get(test,"time")
    
    if(length(times)<length(dates)){
      dates = dates[-which(substr(dates,6,10)=="02-29")]
    }
    years = 2070:2099
    
    array1 = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
    for(y in 1:length(years)){
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
      temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
      array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
      array1[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array1[,,y])
      array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
      array2[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array2[,,y])
      temprain = ifelse(temp>1,1,0)
      array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
      array3[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array3[,,y])
      message("Calcs for year ",years[y]," complete")
    }
    
    nc_close(test)
    
    outcode = (i-1)*3+j
    DeltaSDFprcptot[[outcode]] = apply(array1,c(1,2),mean,na.rm=TRUE)
    DeltaSDFrnnmm[[outcode]] = apply(array3,c(1,2),mean,na.rm=TRUE)
    DeltaSDFrx1day[[outcode]] = apply(array2,c(1,2),mean,na.rm=TRUE)
    DeltaSDForder[outcode] = paste(GCMs[i],obs[j],sep="_")
    
    message("Calcs for obs ",obs[j]," complete")
  }
  message("Calcs for GCM ",GCMs[i]," complete")
}

#######

save(list=c("lon","lat"),file="Domain3to5.Rdata")
save(list=c("obs","obsprcptot","obsrnnmm","obsrx1day"),file="Obs3to5climo.Rdata")
save(list=c("GCMs","MHprcptot","MHrnnmm","MHrx1day"),file="MH3to5climo.Rdata")
save(list=c("GCMs","MFprcptot","MFrnnmm","MFrx1day"),file="MF3to5climo.Rdata")

save(list=c("QDMHorder","QDMHprcptot","QDMHrnnmm","QDMHrx1day"),file="QDMH3to5climo.Rdata")
save(list=c("QDMForder","QDMFprcptot","QDMFrnnmm","QDMFrx1day"),file="QDMF3to5climo.Rdata")

save(list=c("DeltaSDForder","DeltaSDFprcptot","DeltaSDFrnnmm","DeltaSDFrx1day"),file="DeltaSDF3to5climo.Rdata")



