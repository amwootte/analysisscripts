
setwd("/home/woot0002/3to5/")
library(ncdf4)
library(maps)
library(maptools)

load("Domain3to5.Rdata")
londom = lon
latdom = lat

GCMs = c("CCSM4","MIROC5","MPI-ESM-LR")
exps = c("r6i1p1","r1i1p1","r1i1p1")

#####
# Gather historical data first

LOCAHprcptot = LOCAHrnnmm = LOCAHrx1day = list()

for(i in 1:length(GCMs)){
  
  MHfile = paste("/data4/data/DS_proj/LOCA/pr/historical/pr_",GCMs[i],"_",exps[i],"_historical.nc",sep="")
  dates = seq(as.Date("1950-01-01"),as.Date("2005-12-31"),by="day")
  test = nc_open(MHfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  
  if(length(times)<length(dates)){
    dates = dates[-which(substr(dates,6,10)=="02-29")]
  }
  years = 1981:2005
  
  lonidx = which(lon>=min(londom) & lon<=max(londom))
  latidx = which(lat>=min(latdom) & lat<=max(latdom))
  array1 = array2 = array3 = array(NA,dim=c(length(lonidx),length(latidx),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  
  varname = paste("pr_",GCMs[i],"_",exps[i],"_historical",sep="")
  
  for(y in 1:length(years)){
    datesidx = which(as.numeric(substr(dates,1,4))==years[y])
    if(i==1){
    temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))
    } else {
      temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))*86400
    }
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
  
  LOCAHprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  LOCAHrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  LOCAHrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}

######
# Grab future data

LOCAFprcptot = LOCAFrnnmm = LOCAFrx1day = list()

for(i in 1:length(GCMs)){
  
  MFfile = paste("/data4/data/DS_proj/LOCA/pr/future/pr_",GCMs[i],"_",exps[i],"_rcp85-2061-2070.nc",sep="")
  test = nc_open(MFfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  
  years = 2070:2099
  lonidx = which(lon>=min(londom) & lon<=max(londom))
  latidx = which(lat>=min(latdom) & lat<=max(latdom))
  array1 = array2 = array3 = array(NA,dim=c(length(lonidx),length(latidx),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  varname = paste("pr_",GCMs[i],"_",exps[i],"_rcp85",sep="")
  
  for(y in 1:length(years)){
    
    if(years[y]==2070){
      MFfile = paste("/data4/data/DS_proj/LOCA/pr/future/pr_",GCMs[i],"_",exps[i],"_rcp85-2061-2070.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("2061-01-01"),as.Date("2070-12-31"),by="day")
      if(length(times)<length(dates)){
        dates = dates[-which(substr(dates,6,10)=="02-29")]
      }
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
      temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))
      nc_close(test)
    }
    
    if(years[y]>=2071 & years[y]<=2080){
      MFfile = paste("/data4/data/DS_proj/LOCA/pr/future/pr_",GCMs[i],"_",exps[i],"_rcp85-2071-2080.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("2071-01-01"),as.Date("2080-12-31"),by="day")
      if(length(times)<length(dates)){
        dates = dates[-which(substr(dates,6,10)=="02-29")]
      }
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
        temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))
        nc_close(test)
    }
    
    if(years[y]>=2081 & years[y]<=2090){
      MFfile = paste("/data4/data/DS_proj/LOCA/pr/future/pr_",GCMs[i],"_",exps[i],"_rcp85-2081-2090.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("2081-01-01"),as.Date("2090-12-31"),by="day")
      if(length(times)<length(dates)){
        dates = dates[-which(substr(dates,6,10)=="02-29")]
      }
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
        temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))
        nc_close(test)
    }
    
    if(years[y]>=2091 & years[y]<=2099){
      MFfile = paste("/data4/data/DS_proj/LOCA/pr/future/pr_",GCMs[i],"_",exps[i],"_rcp85-2091-2100.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("2091-01-01"),as.Date("2100-12-31"),by="day")
      if(length(times)<length(dates)){
        dates = dates[-which(substr(dates,6,10)=="02-29")]
      }
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
        temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))
      
        nc_close(test)
    }
  
    array1[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
    array1[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array1[,,y])
    array2[,,y] = apply(temp,c(1,2),max,na.rm=TRUE)
    array2[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array2[,,y])
    temprain = ifelse(temp>1,1,0)
    array3[,,y] = apply(temprain,c(1,2),sum,na.rm=TRUE)
    array3[,,y] = ifelse(is.na(temp[,,1])==TRUE,NA,array3[,,y])
    message("Calcs for year ",years[y]," complete")
  }
  
  LOCAFprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  LOCAFrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  LOCAFrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}


lonLOCA = lon[lonidx]
latLOCA = lat[latidx]

save(list=c("GCMs","lonLOCA","latLOCA","LOCAHprcptot","LOCAHrnnmm","LOCAHrx1day"),file="LOCAH3to5climo.Rdata")
save(list=c("GCMs","lonLOCA","latLOCA","LOCAFprcptot","LOCAFrnnmm","LOCAFrx1day"),file="LOCAF3to5climo.Rdata")





