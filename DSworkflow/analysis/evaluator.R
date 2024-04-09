rm(list=ls(all=TRUE))
gc()

source("/ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore02/3to5/I35/scripts/analysisfunctions.R")
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
library(PCICt)
library(abind)
library(SPEI)

#####

setwd("/ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore04/data/DS_proj/EAA")

varname = "pr"
futureperiod = c(2070,2099)
GCMs =c("EC-Earth3","INM-CM4-8","INM-CM5-0","KACE-1-0-G","KIOST-ESM","MPI-ESM1-2-HR")

obsfile =paste("input/daymet_v4_daily_na_",varname,"_19800101-20141231.nc",sep="")

histfiles = system(paste("ls /ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore04/data/DS_proj/EAA/output/",varname,"_*historical_r1i1p1_*.nc",sep=""),intern=TRUE)
rcp45files = system(paste("ls /ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore04/data/DS_proj/EAA/output/",varname,"_*ssp245*.nc",sep=""),intern=TRUE)
rcp85files = system(paste("ls /ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore04/data/DS_proj/EAA/output/",varname,"_*ssp585*.nc",sep=""),intern=TRUE)

histfilesGCM = system(paste("ls /ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore04/data/DS_proj/EAA/input/",varname,"_*historical_r1i1p1_*.nc",sep=""),intern=TRUE)
rcp45filesGCM = system(paste("ls /ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore04/data/DS_proj/EAA/input/",varname,"_*ssp245_r1i1p1*.nc",sep=""),intern=TRUE)
rcp85filesGCM = system(paste("ls /ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore04/data/DS_proj/EAA/input/",varname,"_*ssp585_r1i1p1*.nc",sep=""),intern=TRUE)
 
varout = "pr"
tempout = "annual"
unitsin = "mm"

rawclimoramp = "whitetogreen"
errorclimoramp = "bluetored"
changeclimoramp = "browntogreen"

####
# set relevant thresholds if needed

if(varout == "tmax85"){
  TH = 302.59
}

if(varout == "tmax90"){
  TH = 305.372
}

if(varout == "tmax95"){
  TH = 308.15
}

if(varout == "tmax100"){
  TH = 310.928
}

if(varout == "tmin32"){
  TH = 273.15
}

if(varout == "tmin28"){
  TH = 270.928
}

if(varout == "mdrn"){
  TH = 0.254
}

if(varout == "r1mm"){
  TH = 1
}

if(varout == "pr25"){
  TH = 25.4
}

if(varout == "pr50"){
  TH = 50.8
}

####
# point locations for evaluation

locnames = c("J17","J27","BCRAGD")
loclon = c(-98.4325,-99.7844,-99.081085)
loclat = c(29.47889,29.20889,29.739953)

inputlocations=data.frame(locnames,loclon,loclat)

####
# shapefile load

shapefile = "/home/woot0002/shapefiles/AquiferZones"
test = readShapePoly(shapefile)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
projection(test) <- CRS("+init=epsg:2278")
test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#message("got shapefile loaded")


####
# load obs file, calculate obs climatology

nc = nc_open(obsfile)
lon = ncvar_get(nc,"lon")
lat = ncvar_get(nc,"lat")
histdates = seq(as.Date("1980-01-01"),as.Date("2005-12-31"),by="day")
histdates = histdates[-which(substr(histdates,6,10)=="02-29")]
histdates = as.character(histdates)
years = unique(as.numeric(substr(histdates,1,4)))
  obsdata = ncvar_get(nc,nc$var[[varname]]$name,start=c(1,1,1),count=c(-1,-1,-1))
  if(varname=="pr"){
    obsdata=obsdata*86400
  }
nc_close(nc)

####
# determine grid points for each location

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

###
# get cells to use

inputlocations$R = NA
inputlocations$C = NA

for(n in 1:3){
  pointarea = distfunc(inputlocations$loclon[n],inputlocations$loclat[n],modelgrid)
  inputlocations$R[n] = pointarea$R[1]
  inputlocations$C[n] = pointarea$C[1]
}

if(varout=="DM0" | varout=="DM1" | varout=="DM2" | varout=="DM3" | varout=="DM4"){
  yearmon = unique(substr(histdates,1,7))
  monthlydat = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
  for(m in 1:length(yearmon)){
    idx2 = which(substr(histdates,1,7)==yearmon[m])
    monthlydat[,,m] = apply(obsdata[,,idx2],c(1,2),sum,na.rm=TRUE)
  }
  SPImat = monthlydat
  for(r in 1:length(lon)){
    for(c in 1:length(lat)){
      SPImat[r,c,] = spi(monthlydat[r,c,],1)$fitted
    }
  }
}


####
if(tempout == "annual"){
yearly_obs = array(NA,dim=c(length(lon),length(lat),length(years)))
for(y in 1:length(years)){
  idx1 = which(as.numeric(substr(histdates,1,4))==years[y])
  if(varout=="tasmax" | varout=="tasmin"){
    yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),mean,na.rm=TRUE)
  }
  if(varout=="pr"){
    yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),sum,na.rm=TRUE)
  }
  
  if(varout=="DM0"){
    idx1 = which(as.numeric(substr(yearmon,1,4))==years[y])
    DM0 = ifelse(SPImat[,,idx1]<=-0.5 & SPImat[,,idx1]>= -0.79,1,0)
    yearly_obs[,,y] = apply(DM0,c(1,2),sum,na.rm=TRUE)
  }
  
  if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
    tmp = ifelse(obsdata[,,idx1]>=TH,1,0)
    yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  }
  if(varout=="tmin32" | varout=="tmin28"){
    tmp = ifelse(obsdata[,,idx1]<=TH,1,0)
    yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
  }
  if(varout=="rx1day"){
    yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),max,na.rm=TRUE)
  }
}

climoobs = apply(yearly_obs,c(1,2),mean,na.rm=TRUE)
obsclimo = ifelse(is.na(obsdata[,,1])==FALSE,climoobs,NA)

###
# location grab

#J17 = obsdata[inputlocations$R[1],inputlocations$C[1],]
#J27 = obsdata[inputlocations$R[2],inputlocations$C[2],]
#BCRAGD = obsdata[inputlocations$R[3],inputlocations$C[3],]
#obsdataframe = data.frame(histdates,J17,J27,BCRAGD)
#names(obsdataframe) = c("dates","J17","J27","BCRAGD")
}

#####
# historical simulation climatologies

if(tempout=="annual"){
histclimo = array(NA,dim=c(length(lon),length(lat),length(histfiles)))
}
histlocationlist = list()

for(i in 1:length(histfiles)){
  
  nc = nc_open(histfiles[i])
  timeunits = nc$var[[varname]]$dim[[3]]$units
  calendarinfo = nc$var[[varname]]$dim[[3]]$calendar
  varunits = nc$var[[varname]]$units
  times = ncvar_get(nc,"time")
  startdate1 = do.call("c",strsplit(timeunits," "))
  origindate = as.PCICt(startdate1[3],cal=calendarinfo) # PCICt package helps with weird calendars
  filedates = origindate+times*86400
  years = unique(as.numeric(substr(filedates,1,4)))
  obsdata = ncvar_get(nc,nc$var[[varname]]$name,start=c(1,1,1),count=c(-1,-1,-1)) 
  nc_close(nc)
  if(varname=="pr"){
    obsdata=obsdata*86400
  }
  
  
  if(varout=="DM0" | varout=="DM1" | varout=="DM2" | varout=="DM3" | varout=="DM4"){
    yearmon = unique(substr(filedates,1,7))
    monthlydat = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
    for(m in 1:length(yearmon)){
      idx2 = which(substr(filedates,1,7)==yearmon[m])
      monthlydat[,,m] = apply(obsdata[,,idx2],c(1,2),sum,na.rm=TRUE)
    }
    SPImat = monthlydat
    for(r in 1:length(lon)){
      for(c in 1:length(lat)){
        SPImat[r,c,] = spi(monthlydat[r,c,],1)$fitted
      }
      message("Finished for R ",r," / ",length(lon))
    }
    
  }
  
  
  if(tempout == "annual"){
    yearly_obs = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      idx1 = which(as.numeric(substr(filedates,1,4))==years[y])
      if(varout=="tasmax" | varout=="tasmin"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="pr"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="DM0"){
        idx1 = which(as.numeric(substr(yearmon,1,4))==years[y])
        DM0 = ifelse(SPImat[,,idx1]<=-0.5 & SPImat[,,idx1]>= -0.79,1,0)
        yearly_obs[,,y] = apply(DM0,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
        tmp = ifelse(obsdata[,,idx1]>=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmin32" | varout=="tmin28"){
        tmp = ifelse(obsdata[,,idx1]<=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="rx1day"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),max,na.rm=TRUE)
      }
    }
    
    climoobs = apply(yearly_obs,c(1,2),mean,na.rm=TRUE)
    histclimo[,,i] = ifelse(is.na(obsdata[,,1])==FALSE,climoobs,NA)
    
    ###
    # location grab
    
    #J17 = obsdata[inputlocations$R[1],inputlocations$C[1],]
    #J27 = obsdata[inputlocations$R[2],inputlocations$C[2],]
    #BCRAGD = obsdata[inputlocations$R[3],inputlocations$C[3],]
    #histdataframe = data.frame(as.character(filedates),J17,J27,BCRAGD)
    #names(histdataframe) = c("dates","J17","J27","BCRAGD")
    #histlocationlist[[i]] = histdataframe
  }
  
  message("Finished for ",GCMs[i]," / ",length(histfiles))
}

#####
# historical simulation climatologies -GCM

if(tempout=="annual"){
  histGCMclimo = array(NA,dim=c(length(lon),length(lat),length(histfiles)))
}
histGCMlocationlist = list()

for(i in 1:length(histfilesGCM)){
  
  nc = nc_open(histfilesGCM[i])
  timeunits = nc$var[[varname]]$dim[[3]]$units
  calendarinfo = nc$var[[varname]]$dim[[3]]$calendar
  varunits = nc$var[[varname]]$units
  times = ncvar_get(nc,"time")
  startdate1 = do.call("c",strsplit(timeunits," "))
  origindate = as.PCICt(startdate1[3],cal=calendarinfo) # PCICt package helps with weird calendars
  filedates = origindate+times*86400
  years = unique(as.numeric(substr(filedates,1,4)))
  obsdata = ncvar_get(nc,nc$var[[varname]]$name,start=c(1,1,1),count=c(-1,-1,-1)) 
  nc_close(nc)
  if(varname=="pr"){
    obsdata=obsdata*86400
  }
  
  
  if(varout=="DM0" | varout=="DM1" | varout=="DM2" | varout=="DM3" | varout=="DM4"){
    yearmon = unique(substr(filedates,1,7))
    monthlydat = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
    for(m in 1:length(yearmon)){
      idx2 = which(substr(filedates,1,7)==yearmon[m])
      monthlydat[,,m] = apply(obsdata[,,idx2],c(1,2),sum,na.rm=TRUE)
    }
    SPImat = monthlydat
    for(r in 1:length(lon)){
      for(c in 1:length(lat)){
        SPImat[r,c,] = spi(monthlydat[r,c,],1)$fitted
      }
      message("Finished for R ",r," / ",length(lon))
    }
  }
  
  
  if(tempout == "annual"){
    yearly_obs = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      idx1 = which(as.numeric(substr(filedates,1,4))==years[y])
      if(varout=="tasmax" | varout=="tasmin"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="pr"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="DM0"){
        idx1 = which(as.numeric(substr(yearmon,1,4))==years[y])
        DM0 = ifelse(SPImat[,,idx1]<=-0.5 & SPImat[,,idx1]>= -0.79,1,0)
        yearly_obs[,,y] = apply(DM0,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
        tmp = ifelse(obsdata[,,idx1]>=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmin32" | varout=="tmin28"){
        tmp = ifelse(obsdata[,,idx1]<=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="rx1day"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),max,na.rm=TRUE)
      }
    }
    
    climoobs = apply(yearly_obs,c(1,2),mean,na.rm=TRUE)
    histGCMclimo[,,i] = ifelse(is.na(obsdata[,,1])==FALSE,climoobs,NA)
    
    ###
    # location grab
    
    #J17 = obsdata[inputlocations$R[1],inputlocations$C[1],]
    #J27 = obsdata[inputlocations$R[2],inputlocations$C[2],]
    #BCRAGD = obsdata[inputlocations$R[3],inputlocations$C[3],]
    #histdataframe = data.frame(as.character(filedates),J17,J27,BCRAGD)
    #names(histdataframe) = c("dates","J17","J27","BCRAGD")
    #histGCMlocationlist[[i]] = histdataframe
  }
  
  message("Finished for ",GCMs[i]," / ",length(histfilesGCM))
}

#####
# save historical files and location data
save(list=c("histGCMclimo","histclimo","obsclimo","histlocationlist","histGCMlocationlist"),file=paste("/home/woot0002/EAA/",varout,"_CMIP6_historicalclimos_EAAv2_v03062024.Rdata",sep=""))
#,"obsdataframe"
#####
# rcp45 simulation climatologies

if(tempout=="annual"){
  rcp45climo = array(NA,dim=c(length(lon),length(lat),length(rcp45files)))
}
rcp45locationlist = list()

for(i in 1:length(rcp45files)){
  
  nc = nc_open(rcp45files[i])
  timeunits = nc$var[[varname]]$dim[[3]]$units
  calendarinfo = nc$var[[varname]]$dim[[3]]$calendar
  varunits = nc$var[[varname]]$units
  times = ncvar_get(nc,"time")
  startdate1 = do.call("c",strsplit(timeunits," "))
  origindate = as.PCICt(startdate1[3],cal=calendarinfo) # PCICt package helps with weird calendars
  filedates = origindate+times*86400
  
  idxdates = which(as.numeric(substr(filedates,1,4))>=futureperiod[1] & as.numeric(substr(filedates,1,4))<=futureperiod[2])
  filedates = filedates[idxdates]
  years = unique(as.numeric(substr(filedates,1,4)))
  obsdata = ncvar_get(nc,nc$var[[varname]]$name,start=c(1,1,idxdates[1]),count=c(-1,-1,length(idxdates))) 
  nc_close(nc)
  if(varname=="pr"){
    obsdata=obsdata*86400
  }
  
  if(varout=="DM0" | varout=="DM1" | varout=="DM2" | varout=="DM3" | varout=="DM4"){
    yearmon = unique(substr(filedates,1,7))
    monthlydat = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
    for(m in 1:length(yearmon)){
      idx2 = which(substr(filedates,1,7)==yearmon[m])
      monthlydat[,,m] = apply(obsdata[,,idx2],c(1,2),sum,na.rm=TRUE)
    }
    SPImat = monthlydat
    for(r in 1:length(lon)){
      for(c in 1:length(lat)){
        SPImat[r,c,] = spi(monthlydat[r,c,],1)$fitted
      }
    }
  }
  
  if(tempout == "annual"){
    yearly_obs = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      idx1 = which(as.numeric(substr(filedates,1,4))==years[y])
      if(varout=="tasmax" | varout=="tasmin"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="pr"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="DM0"){
        idx1 = which(as.numeric(substr(yearmon,1,4))==years[y])
        DM0 = ifelse(SPImat[,,idx1]<=-0.5 & SPImat[,,idx1]>= -0.79,1,0)
        yearly_obs[,,y] = apply(DM0,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
        tmp = ifelse(obsdata[,,idx1]>=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmin32" | varout=="tmin28"){
        tmp = ifelse(obsdata[,,idx1]<=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="rx1day"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),max,na.rm=TRUE)
      }
    }
    
    climoobs = apply(yearly_obs,c(1,2),mean,na.rm=TRUE)
    rcp45climo[,,i] = ifelse(is.na(obsdata[,,1])==FALSE,climoobs,NA)
    
    ###
    # location grab
    
    #J17 = obsdata[inputlocations$R[1],inputlocations$C[1],]
    #J27 = obsdata[inputlocations$R[2],inputlocations$C[2],]
    #BCRAGD = obsdata[inputlocations$R[3],inputlocations$C[3],]
    #histdataframe = data.frame(as.character(filedates),J17,J27,BCRAGD)
    #names(histdataframe) = c("dates","J17","J27","BCRAGD")
    #rcp45locationlist[[i]] = histdataframe
  }
  
  message("Finished for ",GCMs[i]," / ",length(rcp45files))
}

#####
# rcp45 simulation climatologies - GCM

if(tempout=="annual"){
  rcp45GCMclimo = array(NA,dim=c(length(lon),length(lat),length(rcp45filesGCM)))
}
rcp45GCMlocationlist = list()

for(i in 1:length(rcp45filesGCM)){
  
  nc = nc_open(rcp45filesGCM[i])
  timeunits = nc$var[[varname]]$dim[[3]]$units
  calendarinfo = nc$var[[varname]]$dim[[3]]$calendar
  varunits = nc$var[[varname]]$units
  times = ncvar_get(nc,"time")
  startdate1 = do.call("c",strsplit(timeunits," "))
  origindate = as.PCICt(startdate1[3],cal=calendarinfo) # PCICt package helps with weird calendars
  filedates = origindate+times*86400
  
  idxdates = which(as.numeric(substr(filedates,1,4))>=futureperiod[1] & as.numeric(substr(filedates,1,4))<=futureperiod[2])
  filedates = filedates[idxdates]
  years = unique(as.numeric(substr(filedates,1,4)))
  obsdata = ncvar_get(nc,nc$var[[varname]]$name,start=c(1,1,idxdates[1]),count=c(-1,-1,length(idxdates))) 
  nc_close(nc)
  if(varname=="pr"){
    obsdata=obsdata*86400
  }
  
  if(varout=="DM0" | varout=="DM1" | varout=="DM2" | varout=="DM3" | varout=="DM4"){
    yearmon = unique(substr(filedates,1,7))
    monthlydat = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
    for(m in 1:length(yearmon)){
      idx2 = which(substr(filedates,1,7)==yearmon[m])
      monthlydat[,,m] = apply(obsdata[,,idx2],c(1,2),sum,na.rm=TRUE)
    }
    SPImat = monthlydat
    for(r in 1:length(lon)){
      for(c in 1:length(lat)){
        SPImat[r,c,] = spi(monthlydat[r,c,],1)$fitted
      }
    }
  }
  
  if(tempout == "annual"){
    yearly_obs = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      idx1 = which(as.numeric(substr(filedates,1,4))==years[y])
      if(varout=="tasmax" | varout=="tasmin"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="pr"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="DM0"){
        idx1 = which(as.numeric(substr(yearmon,1,4))==years[y])
        DM0 = ifelse(SPImat[,,idx1]<=-0.5 & SPImat[,,idx1]>= -0.79,1,0)
        yearly_obs[,,y] = apply(DM0,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
        tmp = ifelse(obsdata[,,idx1]>=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmin32" | varout=="tmin28"){
        tmp = ifelse(obsdata[,,idx1]<=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="rx1day"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),max,na.rm=TRUE)
      }
    }
    
    climoobs = apply(yearly_obs,c(1,2),mean,na.rm=TRUE)
    rcp45GCMclimo[,,i] = ifelse(is.na(obsdata[,,1])==FALSE,climoobs,NA)
    
    ###
    # location grab
    
    #J17 = obsdata[inputlocations$R[1],inputlocations$C[1],]
    #J27 = obsdata[inputlocations$R[2],inputlocations$C[2],]
    #BCRAGD = obsdata[inputlocations$R[3],inputlocations$C[3],]
    #histdataframe = data.frame(as.character(filedates),J17,J27,BCRAGD)
    #names(histdataframe) = c("dates","J17","J27","BCRAGD")
    #rcp45GCMlocationlist[[i]] = histdataframe
  }
  
  message("Finished for ",GCMs[i]," / ",length(rcp45filesGCM))
}

#####
# save files

save(list=c("rcp45GCMclimo","rcp45climo","rcp45locationlist","rcp45GCMlocationlist"),file=paste("/home/woot0002/EAA/",varout,"_CMIP6_ssp245climos_",futureperiod[1],"-",futureperiod[2],"_EAAv2_v03062024.Rdata",sep=""))

#####
# rcp85 simulation climatologies

if(tempout=="annual"){
  rcp85climo = array(NA,dim=c(length(lon),length(lat),length(rcp85files)))
}
rcp85locationlist = list()

for(i in 1:length(rcp85files)){
  
  nc = nc_open(rcp85files[i])
  timeunits = nc$var[[varname]]$dim[[3]]$units
  calendarinfo = nc$var[[varname]]$dim[[3]]$calendar
  varunits = nc$var[[varname]]$units
  times = ncvar_get(nc,"time")
  startdate1 = do.call("c",strsplit(timeunits," "))
  origindate = as.PCICt(startdate1[3],cal=calendarinfo) # PCICt package helps with weird calendars
  filedates = origindate+times*86400
  
  idxdates = which(as.numeric(substr(filedates,1,4))>=futureperiod[1] & as.numeric(substr(filedates,1,4))<=futureperiod[2])
  filedates = filedates[idxdates]
  years = unique(as.numeric(substr(filedates,1,4)))
  obsdata = ncvar_get(nc,nc$var[[varname]]$name,start=c(1,1,idxdates[1]),count=c(-1,-1,length(idxdates))) 
  nc_close(nc)
  if(varname=="pr"){
    obsdata=obsdata*86400
  }
  
  
  if(varout=="DM0" | varout=="DM1" | varout=="DM2" | varout=="DM3" | varout=="DM4"){
    yearmon = unique(substr(filedates,1,7))
    monthlydat = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
    for(m in 1:length(yearmon)){
      idx2 = which(substr(filedates,1,7)==yearmon[m])
      monthlydat[,,m] = apply(obsdata[,,idx2],c(1,2),sum,na.rm=TRUE)
    }
    SPImat = monthlydat
    for(r in 1:length(lon)){
      for(c in 1:length(lat)){
        SPImat[r,c,] = spi(monthlydat[r,c,],1)$fitted
      }
    }
  }
  
  if(tempout == "annual"){
    yearly_obs = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      idx1 = which(as.numeric(substr(filedates,1,4))==years[y])
      if(varout=="tasmax" | varout=="tasmin"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="pr"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="DM0"){
        idx1 = which(as.numeric(substr(yearmon,1,4))==years[y])
        DM0 = ifelse(SPImat[,,idx1]<=-0.5 & SPImat[,,idx1]>= -0.79,1,0)
        yearly_obs[,,y] = apply(DM0,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
        tmp = ifelse(obsdata[,,idx1]>=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmin32" | varout=="tmin28"){
        tmp = ifelse(obsdata[,,idx1]<=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="rx1day"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),max,na.rm=TRUE)
      }
    }
    
    climoobs = apply(yearly_obs,c(1,2),mean,na.rm=TRUE)
    rcp85climo[,,i] = ifelse(is.na(obsdata[,,1])==FALSE,climoobs,NA)
    
    ###
    # location grab
    
    #J17 = obsdata[inputlocations$R[1],inputlocations$C[1],]
    #J27 = obsdata[inputlocations$R[2],inputlocations$C[2],]
    #BCRAGD = obsdata[inputlocations$R[3],inputlocations$C[3],]
    #histdataframe = data.frame(as.character(filedates),J17,J27,BCRAGD)
    #names(histdataframe) = c("dates","J17","J27","BCRAGD")
    #rcp85locationlist[[i]] = histdataframe
  }
  
  message("Finished for ",GCMs[i]," / ",length(rcp85files))
}

#####
# rcp85 simulation climatologies - GCM

if(tempout=="annual"){
  rcp85GCMclimo = array(NA,dim=c(length(lon),length(lat),length(rcp85filesGCM)))
}
rcp85GCMlocationlist = list()

for(i in 1:length(rcp85filesGCM)){
  
  nc = nc_open(rcp85filesGCM[i])
  timeunits = nc$var[[varname]]$dim[[3]]$units
  calendarinfo = nc$var[[varname]]$dim[[3]]$calendar
  varunits = nc$var[[varname]]$units
  times = ncvar_get(nc,"time")
  startdate1 = do.call("c",strsplit(timeunits," "))
  origindate = as.PCICt(startdate1[3],cal=calendarinfo) # PCICt package helps with weird calendars
  filedates = origindate+times*86400
  
  idxdates = which(as.numeric(substr(filedates,1,4))>=futureperiod[1] & as.numeric(substr(filedates,1,4))<=futureperiod[2])
  filedates = filedates[idxdates]
  years = unique(as.numeric(substr(filedates,1,4)))
  obsdata = ncvar_get(nc,nc$var[[varname]]$name,start=c(1,1,idxdates[1]),count=c(-1,-1,length(idxdates))) 
  nc_close(nc)
  if(varname=="pr"){
    obsdata=obsdata*86400
  }
  
  if(varout=="DM0" | varout=="DM1" | varout=="DM2" | varout=="DM3" | varout=="DM4"){
    yearmon = unique(substr(filedates,1,7))
    monthlydat = array(NA,dim=c(length(lon),length(lat),length(yearmon)))
    for(m in 1:length(yearmon)){
      idx2 = which(substr(filedates,1,7)==yearmon[m])
      monthlydat[,,m] = apply(obsdata[,,idx2],c(1,2),sum,na.rm=TRUE)
    }
    SPImat = monthlydat
    for(r in 1:length(lon)){
      for(c in 1:length(lat)){
        SPImat[r,c,] = spi(monthlydat[r,c,],1)$fitted
      }
    }
  }
  
  if(tempout == "annual"){
    yearly_obs = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      idx1 = which(as.numeric(substr(filedates,1,4))==years[y])
      if(varout=="tasmax" | varout=="tasmin"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="pr"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="DM0"){
        idx1 = which(as.numeric(substr(yearmon,1,4))==years[y])
        DM0 = ifelse(SPImat[,,idx1]<=-0.5 & SPImat[,,idx1]>= -0.79,1,0)
        yearly_obs[,,y] = apply(DM0,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
        tmp = ifelse(obsdata[,,idx1]>=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tmin32" | varout=="tmin28"){
        tmp = ifelse(obsdata[,,idx1]<=TH,1,0)
        yearly_obs[,,y] = apply(tmp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="rx1day"){
        yearly_obs[,,y] = apply(obsdata[,,idx1],c(1,2),max,na.rm=TRUE)
      }
    }
    
    climoobs = apply(yearly_obs,c(1,2),mean,na.rm=TRUE)
    rcp85GCMclimo[,,i] = ifelse(is.na(obsdata[,,1])==FALSE,climoobs,NA)
    
    ###
    # location grab
    
    #J17 = obsdata[inputlocations$R[1],inputlocations$C[1],]
    #J27 = obsdata[inputlocations$R[2],inputlocations$C[2],]
    #BCRAGD = obsdata[inputlocations$R[3],inputlocations$C[3],]
    #histdataframe = data.frame(as.character(filedates),J17,J27,BCRAGD)
    #names(histdataframe) = c("dates","J17","J27","BCRAGD")
    #rcp85GCMlocationlist[[i]] = histdataframe
  }
  
  message("Finished for ",GCMs[i]," / ",length(rcp85filesGCM))
}

#####
# save files

save(list=c("rcp85GCMclimo","rcp85climo","rcp85locationlist","rcp85GCMlocationlist"),file=paste("/home/woot0002/EAA/",varout,"_CMIP6_ssp585climos_",futureperiod[1],"-",futureperiod[2],"_EAAv2_v03062024.Rdata",sep=""))

#####
# Error Calcs

DSerror = histclimo
GCMerror = histGCMclimo

for(i in 1:length(GCMs)){
  DSerror[,,i] = histclimo[,,i]-obsclimo
  GCMerror[,,i] = histGCMclimo[,,i]-obsclimo
}

#####
# projected change RCP45

DSchange_rcp45 = rcp45climo-histclimo
GCMchange_rcp45 = rcp45GCMclimo-histGCMclimo

#####
# projected change RCP85

DSchange_rcp85 = rcp85climo-histclimo
GCMchange_rcp85 = rcp85GCMclimo-histGCMclimo

#####
# Plot climatology

range(c(histclimo,histGCMclimo,rcp45climo,rcp45GCMclimo,rcp85climo,rcp85GCMclimo,obsclimo),na.rm=TRUE)
diffcolorbar = colorramp(c(histclimo,histGCMclimo,rcp45climo,rcp45GCMclimo,rcp85climo,rcp85GCMclimo,obsclimo),colorchoice=rawclimoramp,type="raw",Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(0,3))

pdf(paste("/home/woot0002/EAA/",varout,"_climo_EAAv2_CMIP6.pdf",sep=""),width=7,height=5,onefile=TRUE)

testsfc=list(x=lon,y=lat,z=obsclimo)
surface(testsfc,type="I",main=paste("Daymetv4 ",varout," 1980-2014",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)

for(i in 1:length(GCMs)){
testsfc=list(x=lon,y=lat,z=histclimo[,,i])
surface(testsfc,type="I",main=paste(GCMs[i],"_EDQMv3_Daymetv4\n ",varout," 1980-2014",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)
}

for(i in 1:length(rcp45files)){
  testsfc=list(x=lon,y=lat,z=rcp45climo[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_EDQMv3_Daymetv4_SSP245\n ",varout," 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

for(i in 1:length(GCMs)){
  testsfc=list(x=lon,y=lat,z=rcp85climo[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_EDQMv3_Daymetv4_SSP585\n ",varout," 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

for(i in 1:length(GCMs)){
  testsfc=list(x=lon,y=lat,z=histGCMclimo[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"\n ",varout," 1980-2014",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

for(i in 1:length(rcp45files)){
  testsfc=list(x=lon,y=lat,z=rcp45GCMclimo[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_SSP245\n ",varout," 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

for(i in 1:length(GCMs)){
  testsfc=list(x=lon,y=lat,z=rcp85GCMclimo[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_SSP585\n ",varout," 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

dev.off()

#####
# Plot error

range(c(DSerror,GCMerror),na.rm=TRUE)
diffcolorbar = colorramp(c(DSerror,GCMerror),colorchoice=errorclimoramp,type="difference",Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(-4,4))

pdf(paste("/home/woot0002/EAA/",varout,"_climo_error_EAAv2_CMIP6_v103062024.pdf",sep=""),width=7,height=5,onefile=TRUE)

for(i in 1:length(GCMs)){
  testsfc=list(x=lon,y=lat,z=DSerror[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_EDQMv3_Daymetv4\n ",varout," error 1980-2014",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

testsfc=list(x=lon,y=lat,z=apply(DSerror,c(1,2),mean,na.rm=TRUE))
surface(testsfc,type="I",main=paste(varout," mean error 1980-2014",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)

for(i in 1:length(GCMs)){
  testsfc=list(x=lon,y=lat,z=GCMerror[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"\n ",varout," error 1980-2014",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

testsfc=list(x=lon,y=lat,z=apply(GCMerror,c(1,2),mean,na.rm=TRUE))
surface(testsfc,type="I",main=paste(varout," GCM mean error 1980-2014",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)

dev.off()

#####
# Plot Change RCP 4.5

range(c(DSchange_rcp45,GCMchange_rcp45,DSchange_rcp85,GCMchange_rcp85),na.rm=TRUE)
diffcolorbar = colorramp(c(DSchange_rcp45,GCMchange_rcp45,DSchange_rcp85,GCMchange_rcp85),colorchoice=changeclimoramp,type="difference",Blimit=20,use_fixed_scale = TRUE,fixed_scale = c(-80,80))

pdf(paste("/home/woot0002/EAA/",varout,"_climo_ssp245_",futureperiod[1],"-",futureperiod[2],"_EAAv2_CMIP6_v03062024_reduced.pdf",sep=""),width=7,height=5,onefile=TRUE)

for(i in 1:length(rcp45files)){
  testsfc=list(x=lon,y=lat,z=DSchange_rcp45[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_EDQMv3_Daymetv4_SSP245\n ",varout," projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

testsfc=list(x=lon,y=lat,z=apply(DSchange_rcp45,c(1,2),mean,na.rm=TRUE))
surface(testsfc,type="I",main=paste(varout," mean projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)

for(i in 1:length(rcp45files)){
  testsfc=list(x=lon,y=lat,z=GCMchange_rcp45[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_SSP245\n ",varout," projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
  message("Plotted ",i," / ",length(rcp45files))
}

testsfc=list(x=lon,y=lat,z=apply(GCMchange_rcp45,c(1,2),mean,na.rm=TRUE))
surface(testsfc,type="I",main=paste(varout," GCM mean projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)

dev.off()

#####
# Plot Change RCP 8.5

pdf(paste("/home/woot0002/EAA/",varout,"_climo_ssp585_",futureperiod[1],"-",futureperiod[2],"_EAAv2_CMIP6_v03062024_reduced.pdf",sep=""),width=7,height=5,onefile=TRUE)

for(i in 1:length(GCMs)){
  testsfc=list(x=lon,y=lat,z=DSchange_rcp85[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_EDQMv3_Daymetv4_SSP585\n ",varout," projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

testsfc=list(x=lon,y=lat,z=apply(DSchange_rcp85,c(1,2),mean,na.rm=TRUE))
surface(testsfc,type="I",main=paste(varout," mean projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)

for(i in 1:length(GCMs)){
  testsfc=list(x=lon,y=lat,z=GCMchange_rcp85[,,i])
  surface(testsfc,type="I",main=paste(GCMs[i],"_SSP585 \n ",varout," projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
}

testsfc=list(x=lon,y=lat,z=apply(GCMchange_rcp85,c(1,2),mean,na.rm=TRUE))
surface(testsfc,type="I",main=paste(varout," GCM mean projected change 2070-2099",sep=""),xlab="Longitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)


dev.off()

#########
# Error calcs, min, max, range, RMSE

# DS
apply(DSerror,3,range,na.rm=TRUE)
apply(DSerror,3,mean,na.rm=TRUE)
sqrt(apply(DSerror^2,3,mean,na.rm=TRUE)) #RMSE

# GCM
apply(GCMerror,3,range,na.rm=TRUE)
apply(GCMerror,3,mean,na.rm=TRUE)
sqrt(apply(GCMerror^2,3,mean,na.rm=TRUE)) #RMSE


#########
# Change Range RCP45 (OR SSP245)

# DS
apply(DSchange_rcp45,3,range,na.rm=TRUE)
apply(DSchange_rcp45,3,mean,na.rm=TRUE)

# GCM
apply(GCMchange_rcp45,3,range,na.rm=TRUE)
apply(GCMchange_rcp45,3,mean,na.rm=TRUE)

#########
# Change Range RCP85 (or SSP585)

# DS
apply(DSchange_rcp85,3,range,na.rm=TRUE)
apply(DSchange_rcp85,3,mean,na.rm=TRUE)

# GCM
apply(GCMchange_rcp85,3,range,na.rm=TRUE)
apply(GCMchange_rcp85,3,mean,na.rm=TRUE)


#########
# Boxplots - entire region

# adjust to tabulate values and projected changes

groupname = rep("daymetv4",length(obsclimo))
period = rep("obs",length(obsclimo))
type=rep("obs",length(obsclimo))
allclimodata = data.frame(groupname,period,type,as.vector(obsclimo))
names(allclimodata) = c("groupname","period","type","vardata")

for(s in 1:3){
  
  endpt = 5
  
  if(s==1){
    periodname="hist"
    dataused = histclimo
    GCMdataused = histGCMclimo
  } 
  if(s==2){
    periodname="ssp245"
    dataused = rcp45climo
    GCMdataused = rcp45GCMclimo
  } 
  if(s==3){
    periodname="ssp585"
    dataused = rcp85climo
    GCMdataused = rcp85GCMclimo
  } 
  
  
  for(i in 1:endpt){
    tmp = dataused[,,i]
    groupname = rep(paste(GCMs[i],"EDQM","daymetv4",periodname,sep="_"),length(tmp))
    period = rep(periodname,length(tmp))
    type=rep("DS",length(tmp))
    tmpframe = data.frame(groupname,period,type,as.vector(tmp))
    names(tmpframe) = c("groupname","period","type","vardata")
    allclimodata=rbind(allclimodata,tmpframe)
    
    tmp = GCMdataused[,,i]
    groupname = rep(paste(GCMs[i],periodname,sep="_"),length(tmp))
    period = rep(periodname,length(tmp))
    type=rep("RAW",length(tmp))
    tmpframe = data.frame(groupname,period,type,as.vector(tmp))
    names(tmpframe) = c("groupname","period","type","vardata")
    allclimodata=rbind(allclimodata,tmpframe)
  }
  message("Finished calcs for ",periodname)
}

### boxplots with ggplot2

allclimodata$groupname = as.character(allclimodata$groupname)

allclimodata$groupname <- factor(as.character(allclimodata$groupname),
                                 levels = c("daymetv4",                          
                                             "CanESM5_EDQM_daymetv4_hist","CanESM5_EDQM_daymetv4_ssp245","CanESM5_EDQM_daymetv4_ssp585",
                                             "CanESM5_hist", "CanESM5_ssp245","CanESM5_ssp585",
                                             "EC-Earth3_EDQM_daymetv4_hist", "EC-Earth3_EDQM_daymetv4_ssp245","EC-Earth3_EDQM_daymetv4_ssp585",
                                             "EC-Earth3_hist","EC-Earth3_ssp245","EC-Earth3_ssp585",
                                             "KACE-1-0-G_EDQM_daymetv4_hist","KACE-1-0-G_EDQM_daymetv4_ssp245","KACE-1-0-G_EDQM_daymetv4_ssp585",
                                             "KACE-1-0-G_hist","KACE-1-0-G_ssp245","KACE-1-0-G_ssp585",
                                             "KIOST-ESM_EDQM_daymetv4_hist","KIOST-ESM_EDQM_daymetv4_ssp245","KIOST-ESM_EDQM_daymetv4_ssp585",
                                             "KIOST-ESM_hist","KIOST-ESM_ssp245","KIOST-ESM_ssp585",
                                             "MPI-ESM1-2-HR_EDQM_daymetv4_hist","MPI-ESM1-2-HR_EDQM_daymetv4_ssp245","MPI-ESM1-2-HR_EDQM_daymetv4_ssp585",
                                             "MPI-ESM1-2-HR_hist","MPI-ESM1-2-HR_ssp245","MPI-ESM1-2-HR_ssp585"),ordered = TRUE)

ggplot(allclimodata, aes(x=groupname, y=vardata, fill=period)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("gray90","gray60","blue","red")) +
  ggtitle(paste(varout," observed and simulated values",sep=""))+xlab("")+ylab(unitsin)
ggsave(paste("/home/woot0002/EAA/climovalues_fulldomain_EAA_CMIP6_",varout,".pdf",sep=""),device = "pdf",width=9,height=7)

######
# similarly to above but for projected changes

for(s in 1:2){
  
 endpt = 5
  
  if(s==1){
    periodname="ssp245"
    dataused = DSchange_rcp45
    GCMdataused = GCMchange_rcp45
  } 
  if(s==2){
    periodname="ssp585"
    dataused = DSchange_rcp85
    GCMdataused = GCMchange_rcp85
  } 
  
  for(i in 1:endpt){
    tmp = dataused[,,i]
    groupname = rep(paste(GCMs[i],"EDQM","daymetv4",periodname,sep="_"),length(tmp))
    period = rep(periodname,length(tmp))
    type=rep("DS",length(tmp))
    tmpframe = data.frame(groupname,period,type,as.vector(tmp))
    names(tmpframe) = c("groupname","period","type","vardata")
    if(i==1 & s==1){
      allclimodata=tmpframe
    } else {
      allclimodata=rbind(allclimodata,tmpframe)
    }
    
    tmp = GCMdataused[,,i]
    groupname = rep(paste(GCMs[i],periodname,sep="_"),length(tmp))
    period = rep(periodname,length(tmp))
    type=rep("RAW",length(tmp))
    tmpframe = data.frame(groupname,period,type,as.vector(tmp))
    names(tmpframe) = c("groupname","period","type","vardata")
    allclimodata=rbind(allclimodata,tmpframe)
  }
  message("Finished calcs for ",periodname)
}

allclimodata$groupname<-as.character(allclimodata$groupname)

allclimodata$groupname <- factor(as.character(allclimodata$groupname),
                                 levels = c("CanESM5_EDQM_daymetv4_ssp245","CanESM5_EDQM_daymetv4_ssp585",
                                            "CanESM5_ssp245","CanESM5_ssp585",
                                            "EC-Earth3_EDQM_daymetv4_ssp245","EC-Earth3_EDQM_daymetv4_ssp585",
                                            "EC-Earth3_ssp245","EC-Earth3_ssp585",
                                            "KACE-1-0-G_EDQM_daymetv4_ssp245","KACE-1-0-G_EDQM_daymetv4_ssp585",
                                            "KACE-1-0-G_ssp245","KACE-1-0-G_ssp585",
                                            "KIOST-ESM_EDQM_daymetv4_ssp245","KIOST-ESM_EDQM_daymetv4_ssp585",
                                            "KIOST-ESM_ssp245","KIOST-ESM_ssp585",
                                            "MPI-ESM1-2-HR_EDQM_daymetv4_ssp245","MPI-ESM1-2-HR_EDQM_daymetv4_ssp585",
                                            "MPI-ESM1-2-HR_ssp245","MPI-ESM1-2-HR_ssp585"),ordered = TRUE)

ggplot(allclimodata, aes(x=groupname, y=vardata, fill=period)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("blue","red")) +
  ggtitle(paste(varout," projected change 2070-2099",sep=""))+xlab("")+ylab(unitsin)
ggsave(paste("/home/woot0002/EAA/projchange_fulldomain_EAA_CMIP6_",varout,".pdf",sep=""),device = "pdf",width=9,height=7)

#########
# Location Plots - Boxplots

###
# rearrange data files

alllocationdata = obsdataframe 
alllocationdata$group = "daymetv4"
alllocationdata$period = "obs"
alllocationdata$type = "obs"

for(i in 1:5){
  tmp = histlocationlist[[i]]
  tmp$group = paste(GCMs[i],"_EDQM_daymetv4_hist",sep="")
  tmp$period="hist"
  tmp$type="DS"
  alllocationdata = rbind(alllocationdata,tmp)
  
  tmp = histGCMlocationlist[[i]]
  tmp$group = paste(GCMs[i],"_hist",sep="")
  tmp$period="hist"
  tmp$type="RAW"
  alllocationdata = rbind(alllocationdata,tmp)
}

for(i in 1:5){
  tmp = rcp45locationlist[[i]]
  tmp$group = paste(GCMs[i],"_EDQM_daymetv4_ssp245",sep="")
  tmp$period="ssp245"
  tmp$type="DS"
  alllocationdata = rbind(alllocationdata,tmp)
  
  tmp = rcp45GCMlocationlist[[i]]
  tmp$group = paste(GCMs[i],"_ssp245",sep="")
  tmp$period="ssp245"
  tmp$type="RAW"
  alllocationdata = rbind(alllocationdata,tmp)
}

for(i in 1:5){
  tmp = rcp85locationlist[[i]]
  tmp$group = paste(GCMs[i],"_EDQM_daymetv4_ssp585",sep="")
  tmp$period="ssp585"
  tmp$type="DS"
  alllocationdata = rbind(alllocationdata,tmp)
  
  tmp = rcp85GCMlocationlist[[i]]
  tmp$group = paste(GCMs[i],"_ssp585",sep="")
  tmp$period="ssp585"
  tmp$type="RAW"
  alllocationdata = rbind(alllocationdata,tmp)
}

#########
#
alllocationdata$group<-as.character(alllocationdata$group)

alllocationdata$group<- factor(as.character(alllocationdata$group),
                                 levels = c(c("daymetv4",                          
                                              "CanESM5_EDQM_daymetv4_hist","CanESM5_EDQM_daymetv4_ssp245","CanESM5_EDQM_daymetv4_ssp585",
                                              "CanESM5_hist", "CanESM5_ssp245","CanESM5_ssp585",
                                              "EC-Earth3_EDQM_daymetv4_hist", "EC-Earth3_EDQM_daymetv4_ssp245","EC-Earth3_EDQM_daymetv4_ssp585",
                                              "EC-Earth3_hist","EC-Earth3_ssp245","EC-Earth3_ssp585",
                                              "KACE-1-0-G_EDQM_daymetv4_hist","KACE-1-0-G_EDQM_daymetv4_ssp245","KACE-1-0-G_EDQM_daymetv4_ssp585",
                                              "KACE-1-0-G_hist","KACE-1-0-G_ssp245","KACE-1-0-G_ssp585",
                                              "KIOST-ESM_EDQM_daymetv4_hist","KIOST-ESM_EDQM_daymetv4_ssp245","KIOST-ESM_EDQM_daymetv4_ssp585",
                                              "KIOST-ESM_hist","KIOST-ESM_ssp245","KIOST-ESM_ssp585",
                                              "MPI-ESM1-2-HR_EDQM_daymetv4_hist","MPI-ESM1-2-HR_EDQM_daymetv4_ssp245","MPI-ESM1-2-HR_EDQM_daymetv4_ssp585",
                                              "MPI-ESM1-2-HR_hist","MPI-ESM1-2-HR_ssp245","MPI-ESM1-2-HR_ssp585")),ordered = TRUE)

plotresults = TRUE
if(plotresults==TRUE){

ggplot(alllocationdata, aes(x=group, y=J17, fill=period)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("gray90","gray60","blue","red")) +
  ggtitle(paste(varname," observed and simulated values",sep=""))+xlab("")+ylab(unitsin)
ggsave(paste("/home/woot0002/EAA/dailyvalues_J17_EAA_CMIP6_",varout,".pdf",sep=""),device = "pdf",width=9,height=7)

ggplot(alllocationdata, aes(x=group, y=J27, fill=period)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("gray90","gray60","blue","red")) +
  ggtitle(paste(varname," observed and simulated values",sep=""))+xlab("")+ylab(unitsin)
ggsave(paste("/home/woot0002/EAA/dailyvalues_J27_EAA_CMIP6_",varout,".pdf",sep=""),device = "pdf",width=9,height=7)

ggplot(alllocationdata, aes(x=group, y=BCRAGD, fill=period)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("gray90","gray60","blue","red")) +
  ggtitle(paste(varname," observed and simulated values",sep=""))+xlab("")+ylab(unitsin)
ggsave(paste("/home/woot0002/EAA/dailyvalues_BCRAGD_EAA_CMIP6_",varout,".pdf",sep=""),device = "pdf",width=9,height=7)
}
#####
# calculation location climo

groups = unique(alllocationdata$group)

climoout = NULL

for(i in 1:length(groups)){
  
  tmp = subset(alllocationdata,group==groups[i])
  
  if(varout=="tasmax" | varout=="tasmin"){
    tmpyearly = aggregate(tmp[,2:4],by=list(year=as.numeric(substr(tmp[,1],1,4))),mean,na.rm=TRUE)
  }
  if(varout=="pr"){
    tmpyearly = aggregate(tmp[,2:4],by=list(year=as.numeric(substr(tmp[,1],1,4))),sum,na.rm=TRUE)
  }
  
  if(varout=="tmax85" | varout=="tmax90" | varout=="tmax95" | varout=="tmax100" | varout=="mdrn" | varout=="r1mm" | varout=="pr25" | varout=="pr50"){
    tmp[,2:4] = ifelse(tmp[,2:4]>=TH,1,0)
    tmpyearly = aggregate(tmp[,2:4],by=list(year=as.numeric(substr(tmp[,1],1,4))),sum,na.rm=TRUE)
  }
  if(varout=="tmin32" | varout=="tmin28"){
    tmp[,2:4] = ifelse(tmp[,2:4]<=TH,1,0)
    tmpyearly = aggregate(tmp[,2:4],by=list(year=as.numeric(substr(tmp[,1],1,4))),sum,na.rm=TRUE)
  }
  if(varout=="rx1day"){
    tmpyearly = aggregate(tmp[,2:4],by=list(year=as.numeric(substr(tmp[,1],1,4))),max,na.rm=TRUE)
  }
  
  tmpclimo = apply(tmpyearly[,2:4],2,mean,na.rm=TRUE)
  tmpclimo = t(data.frame(tmpclimo))
  climoout = rbind(climoout,tmpclimo)
  
  message("Finished ",groups[i])
}

climoout = data.frame(climoout)

climoout$groups = groups

write.table(climoout,paste("/home/woot0002/EAA/",varout,"_locationclimo_CMIP6.csv",sep=""),sep=",",row.names=FALSE)

