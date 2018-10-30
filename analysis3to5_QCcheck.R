####################
# Basic QC check for downscaled data analysis
#
# (1) boxplots for climatology (multiple period) vs. observations and GCMs - this is the focus here in this analysis script.
# (2) maps for maximum, mean, median, minimum, SD for each member 

###
# Setup libraries, source functions, arguments 

library(ncdf4)
library(sp)
library(fields)
source("/home/woot0002/scripts/analysisfunctions.R")
load("/home/woot0002/3to5/Mask3to5.Rdata")
setwd("/home/woot0002/3to5/")

DS = "DeltaSD"
vargroup = "precip"
if(vargroup=="precip") varname="pr"

####
# Part 1 climatology calculations for prcptot (annual and seasonal), rnnmm (annual only), rx1day (annual only)

GCMs = c("CCSM4","MIROC5","MPI-ESM-LR")
GCMnum = 1:3
obscode = c("D","L","P")
exps = c("r6i1p1","r1i1p1","r1i1p1")
obs = c("daymet","livneh","prism")
scencheck = "rcp45"
futureperiod = c(2006,2039)

#######
# Gather obs data first

obsprcptot = obsDJFprcptot = obsMAMprcptot = obsJJAprcptot = obsSONprcptot = obsrnnmm = obsrx1day = list()

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
  
  array1 = array1DJF = array1MAM = array1JJA = array1SON = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    dates2 = dates[datesidx]
    datesidxDJF = which((as.numeric(substr(dates2,6,7))==12 | as.numeric(substr(dates2,6,7))<=2))
    datesidxMAM = which((as.numeric(substr(dates2,6,7))>=3 & as.numeric(substr(dates2,6,7))<=5))
    datesidxJJA = which((as.numeric(substr(dates2,6,7))>=6 & as.numeric(substr(dates2,6,7))<=8))
    datesidxSON = which((as.numeric(substr(dates2,6,7))>=9 & as.numeric(substr(dates2,6,7))<=11))
    
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(mask==0,NA,array1[,,y])
    array1DJF[,,y] = apply(temp[,,datesidxDJF],c(1,2),sum,na.rm=TRUE)
    array1DJF[,,y] = ifelse(mask==0,NA,array1DJF[,,y])
    array1MAM[,,y] = apply(temp[,,datesidxMAM],c(1,2),sum,na.rm=TRUE)
    array1MAM[,,y] = ifelse(mask==0,NA,array1MAM[,,y])
    array1JJA[,,y] = apply(temp[,,datesidxJJA],c(1,2),sum,na.rm=TRUE)
    array1JJA[,,y] = ifelse(mask==0,NA,array1JJA[,,y])
    array1SON[,,y] = apply(temp[,,datesidxSON],c(1,2),sum,na.rm=TRUE)
    array1SON[,,y] = ifelse(mask==0,NA,array1SON[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(mask==0,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(mask==0,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  obsprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  obsDJFprcptot[[i]] = apply(array1DJF,c(1,2),mean,na.rm=TRUE)
  obsMAMprcptot[[i]] = apply(array1MAM,c(1,2),mean,na.rm=TRUE)
  obsJJAprcptot[[i]] = apply(array1JJA,c(1,2),mean,na.rm=TRUE)
  obsSONprcptot[[i]] = apply(array1SON,c(1,2),mean,na.rm=TRUE)
  
  obsrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  obsrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for obs ",obs[i]," complete")
}

#######
# Gather historical GCM data 

MHprcptot = MHDJFprcptot = MHMAMprcptot = MHJJAprcptot = MHSONprcptot = MHrnnmm = MHrx1day = list()

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
  
  array1 = array1DJF = array1MAM = array1JJA = array1SON = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    dates2 = dates[datesidx]
    datesidxDJF = which((as.numeric(substr(dates2,6,7))==12 | as.numeric(substr(dates2,6,7))<=2))
    datesidxMAM = which((as.numeric(substr(dates2,6,7))>=3 & as.numeric(substr(dates2,6,7))<=5))
    datesidxJJA = which((as.numeric(substr(dates2,6,7))>=6 & as.numeric(substr(dates2,6,7))<=8))
    datesidxSON = which((as.numeric(substr(dates2,6,7))>=9 & as.numeric(substr(dates2,6,7))<=11))
    
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(mask==0,NA,array1[,,y])
    array1DJF[,,y] = apply(temp[,,datesidxDJF],c(1,2),sum,na.rm=TRUE)
    array1DJF[,,y] = ifelse(mask==0,NA,array1DJF[,,y])
    array1MAM[,,y] = apply(temp[,,datesidxMAM],c(1,2),sum,na.rm=TRUE)
    array1MAM[,,y] = ifelse(mask==0,NA,array1MAM[,,y])
    array1JJA[,,y] = apply(temp[,,datesidxJJA],c(1,2),sum,na.rm=TRUE)
    array1JJA[,,y] = ifelse(mask==0,NA,array1JJA[,,y])
    array1SON[,,y] = apply(temp[,,datesidxSON],c(1,2),sum,na.rm=TRUE)
    array1SON[,,y] = ifelse(mask==0,NA,array1SON[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(mask==0,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(mask==0,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  
  MHprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  MHDJFprcptot[[i]] = apply(array1DJF,c(1,2),mean,na.rm=TRUE)
  MHMAMprcptot[[i]] = apply(array1MAM,c(1,2),mean,na.rm=TRUE)
  MHJJAprcptot[[i]] = apply(array1JJA,c(1,2),mean,na.rm=TRUE)
  MHSONprcptot[[i]] = apply(array1SON,c(1,2),mean,na.rm=TRUE)
  MHrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  MHrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}

#######
# Gather future GCM data 

MFprcptot = MFDJFprcptot = MFMAMprcptot = MFJJAprcptot = MFSONprcptot = MFrnnmm = MFrx1day = list()

for(i in 1:length(GCMs)){
  MFfile = paste("pr_day_",GCMs[i],"_",scencheck,"_",exps[i],"_SCCSC0p1_20060101-20991231.nc",sep="")
  dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
  test = nc_open(MFfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  
  if(length(times)<length(dates)){
    dates = dates[-which(substr(dates,6,10)=="02-29")]
  }
  years = futureperiod[1]:futureperiod[2]
  
  array1 = array1DJF = array1MAM = array1JJA = array1SON = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    dates2 = dates[datesidx]
    datesidxDJF = which((as.numeric(substr(dates2,6,7))==12 | as.numeric(substr(dates2,6,7))<=2))
    datesidxMAM = which((as.numeric(substr(dates2,6,7))>=3 & as.numeric(substr(dates2,6,7))<=5))
    datesidxJJA = which((as.numeric(substr(dates2,6,7))>=6 & as.numeric(substr(dates2,6,7))<=8))
    datesidxSON = which((as.numeric(substr(dates2,6,7))>=9 & as.numeric(substr(dates2,6,7))<=11))
    
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(mask==0,NA,array1[,,y])
    array1DJF[,,y] = apply(temp[,,datesidxDJF],c(1,2),sum,na.rm=TRUE)
    array1DJF[,,y] = ifelse(mask==0,NA,array1DJF[,,y])
    array1MAM[,,y] = apply(temp[,,datesidxMAM],c(1,2),sum,na.rm=TRUE)
    array1MAM[,,y] = ifelse(mask==0,NA,array1MAM[,,y])
    array1JJA[,,y] = apply(temp[,,datesidxJJA],c(1,2),sum,na.rm=TRUE)
    array1JJA[,,y] = ifelse(mask==0,NA,array1JJA[,,y])
    array1SON[,,y] = apply(temp[,,datesidxSON],c(1,2),sum,na.rm=TRUE)
    array1SON[,,y] = ifelse(mask==0,NA,array1SON[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(mask==0,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(mask==0,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  
  MFprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  MFDJFprcptot[[i]] = apply(array1DJF,c(1,2),mean,na.rm=TRUE)
  MFMAMprcptot[[i]] = apply(array1MAM,c(1,2),mean,na.rm=TRUE)
  MFJJAprcptot[[i]] = apply(array1JJA,c(1,2),mean,na.rm=TRUE)
  MFSONprcptot[[i]] = apply(array1SON,c(1,2),mean,na.rm=TRUE)
  MFrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  MFrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}

#####
# DS historical grab

histfiles = system(paste("ls /data2/3to5/I35/",varname,"/",DS,"/*00_historical*.nc",sep=""),intern=TRUE)

DSHprcptot = DSHDJFprcptot = DSHMAMprcptot = DSHJJAprcptot = DSHSONprcptot = DSHrnnmm = DSHrx1day = list()

for(i in 1:length(histfiles)){
  dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  test = nc_open(histfiles[i])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"times")
  
  if(length(times)<length(dates)){
    dates = dates[-which(substr(dates,6,10)=="02-29")]
  }
  years = 1981:2005
  
  array1 = array1DJF = array1MAM = array1JJA = array1SON = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    dates2 = dates[datesidx]
    datesidxDJF = which((as.numeric(substr(dates2,6,7))==12 | as.numeric(substr(dates2,6,7))<=2))
    datesidxMAM = which((as.numeric(substr(dates2,6,7))>=3 & as.numeric(substr(dates2,6,7))<=5))
    datesidxJJA = which((as.numeric(substr(dates2,6,7))>=6 & as.numeric(substr(dates2,6,7))<=8))
    datesidxSON = which((as.numeric(substr(dates2,6,7))>=9 & as.numeric(substr(dates2,6,7))<=11))
    
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(mask==0,NA,array1[,,y])
    array1DJF[,,y] = apply(temp[,,datesidxDJF],c(1,2),sum,na.rm=TRUE)
    array1DJF[,,y] = ifelse(mask==0,NA,array1DJF[,,y])
    array1MAM[,,y] = apply(temp[,,datesidxMAM],c(1,2),sum,na.rm=TRUE)
    array1MAM[,,y] = ifelse(mask==0,NA,array1MAM[,,y])
    array1JJA[,,y] = apply(temp[,,datesidxJJA],c(1,2),sum,na.rm=TRUE)
    array1JJA[,,y] = ifelse(mask==0,NA,array1JJA[,,y])
    array1SON[,,y] = apply(temp[,,datesidxSON],c(1,2),sum,na.rm=TRUE)
    array1SON[,,y] = ifelse(mask==0,NA,array1SON[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(mask==0,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(mask==0,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  
  DSHprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  DSHDJFprcptot[[i]] = apply(array1DJF,c(1,2),mean,na.rm=TRUE)
  DSHMAMprcptot[[i]] = apply(array1MAM,c(1,2),mean,na.rm=TRUE)
  DSHJJAprcptot[[i]] = apply(array1JJA,c(1,2),mean,na.rm=TRUE)
  DSHSONprcptot[[i]] = apply(array1SON,c(1,2),mean,na.rm=TRUE)
  DSHrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  DSHrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for histfiles ",histfiles[i]," complete")
}



#####
# DS future grab

futfiles = system(paste("ls /data2/3to5/I35/",varname,"/",DS,"/pr_day*00_",scencheck,"*.nc",sep=""),intern=TRUE)
DSFprcptot = DSFDJFprcptot = DSFMAMprcptot = DSFJJAprcptot = DSFSONprcptot = DSFrnnmm = DSFrx1day = list()

for(i in 1:length(futfiles)){
  dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
  test = nc_open(futfiles[i])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  
  if(length(times)<length(dates)){
    dates = dates[-which(substr(dates,6,10)=="02-29")]
  }
  years = futureperiod[1]:futureperiod[2]
  
  array1 = array1DJF = array1MAM = array1JJA = array1SON = array2 = array3 = array(NA,dim=c(length(lon),length(lat),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    dates2 = dates[datesidx]
    datesidxDJF = which((as.numeric(substr(dates2,6,7))==12 | as.numeric(substr(dates2,6,7))<=2))
    datesidxMAM = which((as.numeric(substr(dates2,6,7))>=3 & as.numeric(substr(dates2,6,7))<=5))
    datesidxJJA = which((as.numeric(substr(dates2,6,7))>=6 & as.numeric(substr(dates2,6,7))<=8))
    datesidxSON = which((as.numeric(substr(dates2,6,7))>=9 & as.numeric(substr(dates2,6,7))<=11))
    
    temp = ncvar_get(test,"pr",start=c(1,1,datesidx[1]),count=c(-1,-1,length(datesidx)))*86400
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(mask==0,NA,array1[,,y])
    array1DJF[,,y] = apply(temp[,,datesidxDJF],c(1,2),sum,na.rm=TRUE)
    array1DJF[,,y] = ifelse(mask==0,NA,array1DJF[,,y])
    array1MAM[,,y] = apply(temp[,,datesidxMAM],c(1,2),sum,na.rm=TRUE)
    array1MAM[,,y] = ifelse(mask==0,NA,array1MAM[,,y])
    array1JJA[,,y] = apply(temp[,,datesidxJJA],c(1,2),sum,na.rm=TRUE)
    array1JJA[,,y] = ifelse(mask==0,NA,array1JJA[,,y])
    array1SON[,,y] = apply(temp[,,datesidxSON],c(1,2),sum,na.rm=TRUE)
    array1SON[,,y] = ifelse(mask==0,NA,array1SON[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(mask==0,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(mask==0,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  nc_close(test)
  
  DSFprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  DSFDJFprcptot[[i]] = apply(array1DJF,c(1,2),mean,na.rm=TRUE)
  DSFMAMprcptot[[i]] = apply(array1MAM,c(1,2),mean,na.rm=TRUE)
  DSFJJAprcptot[[i]] = apply(array1JJA,c(1,2),mean,na.rm=TRUE)
  DSFSONprcptot[[i]] = apply(array1SON,c(1,2),mean,na.rm=TRUE)
  DSFrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  DSFrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for futfiles ",futfiles[i]," complete")
}

######################

#save(list=c("obs","obsprcptot","obsDJFprcptot","obsMAMprcptot","obsJJAprcptot","obsSONprcptot","obsrnnmm","obsrx1day"),file="Obs3to5climo.Rdata")
#save(list=c("GCMs","MHprcptot","MHDJFprcptot","MHMAMprcptot","MHJJAprcptot","MHSONprcptot","MHrnnmm","MHrx1day"),file="MH3to5climo.Rdata")
save(list=c("GCMs","MFprcptot","MFDJFprcptot","MFMAMprcptot","MFJJAprcptot","MFSONprcptot","MFrnnmm","MFrx1day"),file=paste("MF3to5climo",futureperiod[1],"-",futureperiod[2],".Rdata",sep=""))

#save(list=c("histfiles","DSHprcptot","DSHDJFprcptot","DSHMAMprcptot","DSHJJAprcptot","DSHSONprcptot","DSHrnnmm","DSHrx1day"),file=paste(DS,"H3to5climo.Rdata",sep=""))
save(list=c("futfiles","DSFprcptot","DSFDJFprcptot","DSFMAMprcptot","DSFJJAprcptot","DSFSONprcptot","DSFrnnmm","DSFrx1day"),file=paste(DS,"F3to5climo",futureperiod[1],"-",futureperiod[2],".Rdata",sep=""))

########################

prcptot = DJFprcptot = MAMprcptot = JJAprcptot = SONprcptot = rnnmm = rx1day = c()

for(i in 1:3){
  prcptot = c(prcptot,obsprcptot[[i]])
  DJFprcptot = c(DJFprcptot,obsDJFprcptot[[i]])
  MAMprcptot = c(MAMprcptot,obsMAMprcptot[[i]])
  JJAprcptot = c(JJAprcptot,obsJJAprcptot[[i]])
  SONprcptot = c(SONprcptot,obsSONprcptot[[i]])
  rnnmm = c(rnnmm,obsrnnmm[[i]])
  rx1day = c(rx1day,obsrx1day[[i]])
}

models = obs
model = rep(models,each=length(rx1day)/3)
outputframe = data.frame(model,prcptot,DJFprcptot,MAMprcptot,JJAprcptot,SONprcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = outputframe[-NAidx,]

###

prcptot = DJFprcptot = MAMprcptot = JJAprcptot = SONprcptot = rnnmm = rx1day = c()

for(i in 1:3){
  prcptot = c(prcptot,MHprcptot[[i]])
  DJFprcptot = c(DJFprcptot,MHDJFprcptot[[i]])
  MAMprcptot = c(MAMprcptot,MHMAMprcptot[[i]])
  JJAprcptot = c(JJAprcptot,MHJJAprcptot[[i]])
  SONprcptot = c(SONprcptot,MHSONprcptot[[i]])
  rnnmm = c(rnnmm,MHrnnmm[[i]])
  rx1day = c(rx1day,MHrx1day[[i]])
}

models = paste(GCMs,"-H",sep="")
model = rep(models,each=length(rx1day)/3)
outputframe = data.frame(model,prcptot,DJFprcptot,MAMprcptot,JJAprcptot,SONprcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

###

prcptot = DJFprcptot = MAMprcptot = JJAprcptot = SONprcptot = rnnmm = rx1day = c()

for(i in 1:3){
  prcptot = c(prcptot,MFprcptot[[i]])
  DJFprcptot = c(DJFprcptot,MFDJFprcptot[[i]])
  MAMprcptot = c(MAMprcptot,MFMAMprcptot[[i]])
  JJAprcptot = c(JJAprcptot,MFJJAprcptot[[i]])
  SONprcptot = c(SONprcptot,MFSONprcptot[[i]])
  rnnmm = c(rnnmm,MFrnnmm[[i]])
  rx1day = c(rx1day,MFrx1day[[i]])
}

models = paste(GCMs,"-F",sep="")
model = rep(models,each=length(rx1day)/3)
outputframe = data.frame(model,prcptot,DJFprcptot,MAMprcptot,JJAprcptot,SONprcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

####

prcptot = DJFprcptot = MAMprcptot = JJAprcptot = SONprcptot = rnnmm = rx1day = c()

for(i in 1:9){
  prcptot = c(prcptot,DSHprcptot[[i]])
  DJFprcptot = c(DJFprcptot,DSHDJFprcptot[[i]])
  MAMprcptot = c(MAMprcptot,DSHMAMprcptot[[i]])
  JJAprcptot = c(JJAprcptot,DSHJJAprcptot[[i]])
  SONprcptot = c(SONprcptot,DSHSONprcptot[[i]])
  rnnmm = c(rnnmm,DSHrnnmm[[i]])
  rx1day = c(rx1day,DSHrx1day[[i]])
}

models = paste(DS,paste(rep(GCMs,each=3),rep(obs,3),sep="-"),"H",sep="-")
model = rep(models,each=length(rx1day)/9)
outputframe = data.frame(model,prcptot,DJFprcptot,MAMprcptot,JJAprcptot,SONprcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

####

prcptot = DJFprcptot = MAMprcptot = JJAprcptot = SONprcptot = rnnmm = rx1day = c()

for(i in 1:9){
  prcptot = c(prcptot,DSFprcptot[[i]])
  DJFprcptot = c(DJFprcptot,DSFDJFprcptot[[i]])
  MAMprcptot = c(MAMprcptot,DSFMAMprcptot[[i]])
  JJAprcptot = c(JJAprcptot,DSFJJAprcptot[[i]])
  SONprcptot = c(SONprcptot,DSFSONprcptot[[i]])
  rnnmm = c(rnnmm,DSFrnnmm[[i]])
  rx1day = c(rx1day,DSFrx1day[[i]])
}

models = paste(DS,paste(rep(GCMs,each=3),rep(obs,3),sep="-"),"F",sep="-")
model = rep(models,each=length(rx1day)/9)
outputframe = data.frame(model,prcptot,DJFprcptot,MAMprcptot,JJAprcptot,SONprcptot,rnnmm,rx1day)
NAidx = which(is.na(outputframe$prcptot)==TRUE)
alldat = rbind(alldat,outputframe[-NAidx,])

########

alldat$plotorder = NA
alldat$plotorder[which(alldat$model=="daymet")]=1
alldat$plotorder[which(alldat$model=="livneh")]=2
alldat$plotorder[which(alldat$model=="prism")]=3

alldat$plotorder[which(alldat$model=="CCSM4-H")]=4
alldat$plotorder[which(alldat$model=="MIROC5-H")]=5
alldat$plotorder[which(alldat$model=="MPI-ESM-LR-H")]=6

alldat$plotorder[which(alldat$model=="CCSM4-F")]=7
alldat$plotorder[which(alldat$model=="MIROC5-F")]=8
alldat$plotorder[which(alldat$model=="MPI-ESM-LR-F")]=9

alldat$plotorder[which(alldat$model==paste(DS,"-CCSM4-daymet-H",sep=""))]=10
alldat$plotorder[which(alldat$model==paste(DS,"-MIROC5-daymet-H",sep=""))]=11
alldat$plotorder[which(alldat$model==paste(DS,"-MPI-ESM-LR-daymet-H",sep=""))]=12

alldat$plotorder[which(alldat$model==paste(DS,"-CCSM4-daymet-F",sep=""))]=13
alldat$plotorder[which(alldat$model==paste(DS,"-MIROC5-daymet-F",sep=""))]=14
alldat$plotorder[which(alldat$model==paste(DS,"-MPI-ESM-LR-daymet-F",sep=""))]=15

alldat$plotorder[which(alldat$model==paste(DS,"-CCSM4-livneh-H",sep=""))]=16
alldat$plotorder[which(alldat$model==paste(DS,"-MIROC5-livneh-H",sep=""))]=17
alldat$plotorder[which(alldat$model==paste(DS,"-MPI-ESM-LR-livneh-H",sep=""))]=18

alldat$plotorder[which(alldat$model==paste(DS,"-CCSM4-livneh-F",sep=""))]=19
alldat$plotorder[which(alldat$model==paste(DS,"-MIROC5-livneh-F",sep=""))]=20
alldat$plotorder[which(alldat$model==paste(DS,"-MPI-ESM-LR-livneh-F",sep=""))]=21

alldat$plotorder[which(alldat$model==paste(DS,"-CCSM4-prism-H",sep=""))]=22
alldat$plotorder[which(alldat$model==paste(DS,"-MIROC5-prism-H",sep=""))]=23
alldat$plotorder[which(alldat$model==paste(DS,"-MPI-ESM-LR-prism-H",sep=""))]=24

alldat$plotorder[which(alldat$model==paste(DS,"-CCSM4-prism-F",sep=""))]=25
alldat$plotorder[which(alldat$model==paste(DS,"-MIROC5-prism-F",sep=""))]=26
alldat$plotorder[which(alldat$model==paste(DS,"-MPI-ESM-LR-prism-F",sep=""))]=27

####

labelsin = c(obs,paste(GCMs,"H",sep="-"),paste(GCMs,"F",sep="-"),paste(DS,GCMs,obs[1],"H",sep="-"),paste(DS,GCMs,obs[1],"F",sep="-"),paste(DS,GCMs,obs[2],"H",sep="-"),paste(DS,GCMs,obs[2],"F",sep="-"),paste(DS,GCMs,obs[3],"H",sep="-"),paste(DS,GCMs,obs[3],"F",sep="-"))
locations = c(1,2,3,5,6,7,9,10,11,13,14,15,17,18,19,21,22,23,25,26,27,29,30,31,33,34,35)
cols = c("red","purple","magenta",rep(c("brown","cyan","yellow"),8))

pdf(paste("QCraw_",DS,"_",scencheck,"_",futureperiod[1],"-",futureperiod[2],".pdf",sep=""),onefile=TRUE,width=9,height=9)
boxplot(prcptot~plotorder,data=alldat,main="prcptot Obs, GCM, and DS Values across I35 region grid points",at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(DJFprcptot~plotorder,data=alldat,main="DJFprcptot Obs, GCM, and DS Values across I35 region grid points",at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(MAMprcptot~plotorder,data=alldat,main="MAMprcptot Obs, GCM, and DS Values across I35 region grid points",at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(JJAprcptot~plotorder,data=alldat,main="JJAprcptot Obs, GCM, and DS Values across I35 region grid points",at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(SONprcptot~plotorder,data=alldat,main="SONprcptot Obs, GCM, and DS Values across I35 region grid points",at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(rnnmm~plotorder,data=alldat,main="rnnmm Obs, GCM, and DS Values across I35 region grid points",at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(rx1day~plotorder,data=alldat,main="rnnmm Obs, GCM, and DS Values across I35 region grid points",at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)
dev.off()
####

alldat1 = alldat[which(alldat$model=="CCSM4-F"),2:9]-alldat[which(alldat$model=="CCSM4-H"),2:9]
alldat1$model = "CCSM4"
alldat2 = alldat[which(alldat$model=="MIROC5-F"),2:9]-alldat[which(alldat$model=="MIROC5-H"),2:9]
alldat2$model = "MIROC5"
alldat3 = alldat[which(alldat$model=="MPI-ESM-LR-F"),2:9]-alldat[which(alldat$model=="MPI-ESM-LR-H"),2:9]
alldat3$model = "MPI-ESM-LR"
alldatv2 = rbind(alldat1,alldat2)
alldatv2 = rbind(alldatv2,alldat3)

alldat4 = alldat[which(alldat$model==paste(DS,"-CCSM4-daymet-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-CCSM4-daymet-H",sep="")),2:9]
alldat4$model = paste(DS,"-CCSM4-daymet",sep="")
alldat5 = alldat[which(alldat$model==paste(DS,"-MIROC5-daymet-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-MIROC5-daymet-H",sep="")),2:9]
alldat5$model = paste(DS,"-MIROC5-daymet",sep="")
alldat6 = alldat[which(alldat$model==paste(DS,"-MPI-ESM-LR-daymet-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-MPI-ESM-LR-daymet-H",sep="")),2:9]
alldat6$model = paste(DS,"-MPI-ESM-LR-daymet",sep="")
alldatv2 = rbind(alldatv2,alldat4)
alldatv2 = rbind(alldatv2,alldat5)
alldatv2 = rbind(alldatv2,alldat6)

alldat7 = alldat[which(alldat$model==paste(DS,"-CCSM4-livneh-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-CCSM4-livneh-H",sep="")),2:9]
alldat7$model = paste(DS,"-CCSM4-livneh",sep="")
alldat8 = alldat[which(alldat$model==paste(DS,"-MIROC5-livneh-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-MIROC5-livneh-H",sep="")),2:9]
alldat8$model = paste(DS,"-MIROC5-livneh",sep="")
alldat9 = alldat[which(alldat$model==paste(DS,"-MPI-ESM-LR-livneh-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-MPI-ESM-LR-livneh-H",sep="")),2:9]
alldat9$model = paste(DS,"-MPI-ESM-LR-livneh",sep="")
alldatv2 = rbind(alldatv2,alldat7)
alldatv2 = rbind(alldatv2,alldat8)
alldatv2 = rbind(alldatv2,alldat9)

alldat10 = alldat[which(alldat$model==paste(DS,"-CCSM4-prism-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-CCSM4-prism-H",sep="")),2:9]
alldat10$model = paste(DS,"-CCSM4-prism",sep="")
alldat11 = alldat[which(alldat$model==paste(DS,"-MIROC5-prism-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-MIROC5-prism-H",sep="")),2:9]
alldat11$model = paste(DS,"-MIROC5-prism",sep="")
alldat12 = alldat[which(alldat$model==paste(DS,"-MPI-ESM-LR-prism-F",sep="")),2:9]-alldat[which(alldat$model==paste(DS,"-MPI-ESM-LR-prism-H",sep="")),2:9]
alldat12$model = paste(DS,"-MPI-ESM-LR-prism",sep="")
alldatv2 = rbind(alldatv2,alldat10)
alldatv2 = rbind(alldatv2,alldat11)
alldatv2 = rbind(alldatv2,alldat12)

#####

alldatv2$plotorder = NA

alldatv2$plotorder[which(alldatv2$model=="CCSM4")]=1
alldatv2$plotorder[which(alldatv2$model=="MIROC5")]=2
alldatv2$plotorder[which(alldatv2$model=="MPI-ESM-LR")]=3

alldatv2$plotorder[which(alldatv2$model==paste(DS,"-CCSM4-daymet",sep=""))]=4
alldatv2$plotorder[which(alldatv2$model==paste(DS,"-MIROC5-daymet",sep=""))]=5
alldatv2$plotorder[which(alldatv2$model==paste(DS,"-MPI-ESM-LR-daymet",sep=""))]=6

alldatv2$plotorder[which(alldatv2$model==paste(DS,"-CCSM4-livneh",sep=""))]=7
alldatv2$plotorder[which(alldatv2$model==paste(DS,"-MIROC5-livneh",sep=""))]=8
alldatv2$plotorder[which(alldatv2$model==paste(DS,"-MPI-ESM-LR-livneh",sep=""))]=9

alldatv2$plotorder[which(alldatv2$model==paste(DS,"-CCSM4-prism",sep=""))]=10
alldatv2$plotorder[which(alldatv2$model==paste(DS,"-MIROC5-prism",sep=""))]=11
alldatv2$plotorder[which(alldatv2$model==paste(DS,"-MPI-ESM-LR-prism",sep=""))]=12

###

labelsin = c(GCMs,paste(DS,rep(GCMs,3),rep(obs,each=3),sep="-"))
locations = c(1,2,3,5,6,7,9,10,11,13,14,15)
cols = rep(c("brown","cyan","yellow"),4)

pdf(paste("QCchange_",DS,"_",scencheck,"_",futureperiod[1],"-",futureperiod[2],".pdf",sep=""),onefile=TRUE,width=9,height=9)

boxplot(prcptot~plotorder,data=alldatv2,main=paste("prcptot ",futureperiod[1],"-",futureperiod[2]," minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",sep=""),at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(DJFprcptot~plotorder,data=alldatv2,main=paste("DJFprcptot ",futureperiod[1],"-",futureperiod[2]," minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",sep=""),at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(MAMprcptot~plotorder,data=alldatv2,main=paste("MAMprcptot ",futureperiod[1],"-",futureperiod[2]," minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",sep=""),at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(JJAprcptot~plotorder,data=alldatv2,main=paste("JJAprcptot ",futureperiod[1],"-",futureperiod[2]," minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",sep=""),at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(SONprcptot~plotorder,data=alldatv2,main=paste("SONprcptot ",futureperiod[1],"-",futureperiod[2]," minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",sep=""),at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(rnnmm~plotorder,data=alldatv2,main=paste("rnnmm ",futureperiod[1],"-",futureperiod[2]," minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",sep=""),at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)

boxplot(rx1day~plotorder,data=alldatv2,main=paste("rx1day ",futureperiod[1],"-",futureperiod[2]," minus 1981-2005 \n Values across I35 region grid points (Grouped by Obs)",sep=""),at=locations,ylab="mm",xaxt="n",col=cols)
text(x =  locations, y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labelsin, xpd = TRUE,cex=0.5)
dev.off()


