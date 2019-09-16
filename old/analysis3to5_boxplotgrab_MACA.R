
setwd("/home/woot0002/3to5/")
library(ncdf4)
library(maps)
library(maptools)

load("Domain3to5.Rdata")
londom = lon
latdom = lat

GCMs = c("CCSM4","MIROC5")
exps = c("r6i1p1","r1i1p1")

#####
# Gather historical data first

MACAHprcptot = MACAHrnnmm = MACAHrx1day = list()

for(i in 1:length(GCMs)){
  #macav2livneh_pr_",GCMs[i],"_",exps[i],"_historical_1970_1989_CONUS_daily.nc
  MFfile = paste("/data4/data/DS_proj/MACA_LIVNEH/pr/historical/macav2livneh_pr_",GCMs[i],"_",exps[i],"_historical_1970_1989_CONUS_daily.nc",sep="")
  test = nc_open(MFfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  
  years = 1981:2005
  lonidx = which(lon>=min(londom) & lon<=max(londom))
  latidx = which(lat>=min(latdom) & lat<=max(latdom))
  array1 = array2 = array3 = array(NA,dim=c(length(lonidx),length(latidx),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  varname = "precipitation"
  
  for(y in 1:length(years)){
    
    if(years[y]<=1989){
      MFfile = paste("/data4/data/DS_proj/MACA_LIVNEH/pr/historical/macav2livneh_pr_",GCMs[i],"_",exps[i],"_historical_1970_1989_CONUS_daily.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("1970-01-01"),as.Date("1989-12-31"),by="day")
      if(length(times)<length(dates)){
        dates = dates[-which(substr(dates,6,10)=="02-29")]
      }
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
      temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))
      nc_close(test)
    }
    
    if(years[y]>=1990 & years[y]<=2005){
      MFfile = paste("/data4/data/DS_proj/MACA_LIVNEH/pr/historical/macav2livneh_pr_",GCMs[i],"_",exps[i],"_historical_1990_2005_CONUS_daily.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("1990-01-01"),as.Date("2005-12-31"),by="day")
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
  
  MACAHprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  MACAHrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  MACAHrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}

######
# Grab future data

MACAFprcptot = MACAFrnnmm = MACAFrx1day = list()

for(i in 1:length(GCMs)){
  #macav2livneh_pr_",GCMs[i],"_",exps[i],"_historical_1970_1989_CONUS_daily.nc
  MFfile = paste("/data4/data/DS_proj/MACA_LIVNEH/pr/future/macav2livneh_pr_",GCMs[i],"_",exps[i],"_rcp85_2066_2085_CONUS_daily.nc",sep="")
  test = nc_open(MFfile)
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  
  years = 2070:2099
  lonidx = which(lon>=min(londom) & lon<=max(londom))
  latidx = which(lat>=min(latdom) & lat<=max(latdom))
  array1 = array2 = array3 = array(NA,dim=c(length(lonidx),length(latidx),length(years))) # array1 is prcptotal, array2 is rx1day, array3 is rnnmm
  varname = "precipitation"
  
  for(y in 1:length(years)){
    
    if(years[y]<=2085){
      MFfile = paste("/data4/data/DS_proj/MACA_LIVNEH/pr/future/macav2livneh_pr_",GCMs[i],"_",exps[i],"_rcp85_2066_2085_CONUS_daily.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("2066-01-01"),as.Date("2085-12-31"),by="day")
      if(length(times)<length(dates)){
        dates = dates[-which(substr(dates,6,10)=="02-29")]
      }
      datesidx = which(as.numeric(substr(dates,1,4))==years[y])
      temp = ncvar_get(test,varname,start=c(lonidx[1],latidx[1],datesidx[1]),count=c(length(lonidx),length(latidx),length(datesidx)))
      nc_close(test)
    }
    
    if(years[y]>=2086 & years[y]<=2099){
      MFfile = paste("/data4/data/DS_proj/MACA_LIVNEH/pr/future/macav2livneh_pr_",GCMs[i],"_",exps[i],"_rcp85_2086_2099_CONUS_daily.nc",sep="")
      test = nc_open(MFfile)
      times = ncvar_get(test,"time")
      dates = seq(as.Date("2086-01-01"),as.Date("2099-12-31"),by="day")
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
  
  MACAFprcptot[[i]] = apply(array1,c(1,2),mean,na.rm=TRUE)
  MACAFrnnmm[[i]] = apply(array3,c(1,2),mean,na.rm=TRUE)
  MACAFrx1day[[i]] = apply(array2,c(1,2),mean,na.rm=TRUE)
  message("Calcs for GCM ",GCMs[i]," complete")
}


lonMACA = lon[lonidx]
latMACA = lat[latidx]

save(list=c("GCMs","lonMACA","latMACA","MACAHprcptot","MACAHrnnmm","MACAHrx1day"),file="MACAH3to5climo.Rdata")
save(list=c("GCMs","lonMACA","latMACA","MACAFprcptot","MACAFrnnmm","MACAFrx1day"),file="MACAF3to5climo.Rdata")





