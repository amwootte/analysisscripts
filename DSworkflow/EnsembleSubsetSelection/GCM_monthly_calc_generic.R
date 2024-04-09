rm(list=ls(all=TRUE))
gc()

### monthly climatologies are for historical and projected for each GCM. 

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
library(PCICt)

GCMs = c("ACCESS-CM2","ACCESS-ESM1-5","CanESM5","CMCC-ESM2","EC-Earth3-CC","EC-Earth3",
         "EC-Earth3-Veg-LR","FGOALS-g3","GFDL-CM4","GFDL-ESM4","INM-CM4-8","INM-CM5-0",
         "IPSL-CM6A-LR","KACE-1-0-G","KIOST-ESM","MIROC6","MPI-ESM1-2-HR","MPI-ESM1-2-LR",
         "MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","TaiESM1")

histperiod = c(1980,2014)
projperiod = c(2070,2099)
varname="pr"
unitsout="mm"
rcpout = "ssp585"

runhist = TRUE
runproj = TRUE

filelisthist = filelistproj = c()

for(i in 1:length(GCMs)){
  
  histcommand = paste("ls /data4/data/CMIP6/",varname,"/historical/*",GCMs[i],"*.nc",sep="") # data files for GCMs are taken from ESGF repository and held in local repository for analysis
  projcommand = paste("ls /data4/data/CMIP6/",varname,"/",rcpout,"/*",GCMs[i],"*.nc",sep="")
  
  filelisthist = c(filelisthist,system(histcommand,intern=TRUE))
  filelistproj = c(filelistproj,system(projcommand,intern=TRUE))
  
}

tab1 = do.call("rbind",strsplit(filelisthist,split="_",fixed=TRUE))
tab2 = do.call("rbind",strsplit(tab1[,ncol(tab1)],split="-",fixed=TRUE))
tab2[,ncol(tab2)]=substr(tab2[,ncol(tab2)],1,8)
filelist1 = cbind(tab1[,c(3,4,5,6)],tab2[,c(1,2)])
filelist1 = as.data.frame(filelist1)
names(filelist1) = c("model","scen","experiment","grid","starttime","endtime")
filelist1$path = filelisthist
filelist1$startyear = as.numeric(substr(filelist1$starttime,1,4))
filelist1$endyear = as.numeric(substr(filelist1$endtime,1,4))
fileshist = filelist1


tab1 = do.call("rbind",strsplit(filelistproj,split="_",fixed=TRUE))
tab2 = do.call("rbind",strsplit(tab1[,ncol(tab1)],split="-",fixed=TRUE))
tab2[,ncol(tab2)]=substr(tab2[,ncol(tab2)],1,8)
filelist1 = cbind(tab1[,c(3,4,5,6)],tab2[,c(1,2)])
filelist1 = as.data.frame(filelist1)
names(filelist1) = c("model","scen","experiment","grid","starttime","endtime")
filelist1$path = filelistproj
filelist1$startyear = as.numeric(substr(filelist1$starttime,1,4))
filelist1$endyear = as.numeric(substr(filelist1$endtime,1,4))
filesproj = filelist1

#####
# historical monthly calculations

if(runhist==TRUE){

histyears = histperiod[1]:histperiod[2]

for(i in 1:length(GCMs)){
  
  tmplist = fileshist[which(fileshist$model==GCMs[i]),]
  filesub = which(tmplist$endyear>=histperiod[1])
  tmplist2 = tmplist[filesub,]
  
  fileorder = order(tmplist2$endyear)
  tmplist2=tmplist2[fileorder,]
  
  filesuse = unique(tmplist2$path)
  
  tmplist3 = NULL
  for(j in 1:length(filesuse)){
    tmpidx = which(tmplist2$path==filesuse[j])
    tplist = tmplist2[tmpidx,]
    tmplist3 = rbind(tmplist3,tplist[1,])
  }
  
  for(k in 1:length(filesuse)){
   
    nctest = nc_open(tmplist3$path[k]) 
    
    times = ncvar_get(nctest,"time")
    
    if(k==1){
      
      lon = ncvar_get(nctest,"lon")
      lat = ncvar_get(nctest,"lat")
      calendarinfo = nctest$var[[varname]]$dim[[3]]$calendar
      varnames = c()
      for(m in 1:length(nctest$var)){
        varnames[m] = nctest$var[[m]]$name
      }
      varidx = which(varnames==varname)
      timeunits = nctest$var[[varidx]]$dim[[3]]$units
      latunits = nctest$var[[varidx]]$dim[[2]]$units
      lonunits = nctest$var[[varidx]]$dim[[1]]$units
      
      indexdate = substr(timeunits,12,21)
      origindate = as.PCICt(indexdate,cal=calendarinfo) # PCICt package helps with weird calendars
      
      
      timeunits = paste(substr(timeunits,1,22),"12:00:00",sep=" ")
      
      
      #datesall = seq(as.Date(paste(tmplist3$startyear[which(tmplist3$startyear==min(tmplist3$startyear))],"-01-01",sep="")),as.Date(paste(tmplist3$endyear[which(tmplist3$endyear==max(tmplist3$endyear))],"-12-31",sep="")),by="day")
      yearmonall = seq(as.Date(paste(histyears[1],"-01-15",sep="")),as.Date(paste(histyears[length(histyears)],"-12-15",sep="")),by="month")
      timeoutmonall = c()
      
      outputarray1 = outputarray2 = array(NA,dim=c(length(lon),length(lat),length(yearmonall)))
      
    }
    
    dates = origindate+times*86400
    
    if(as.numeric(substr(dates[1],1,4))<as.numeric(tmplist3$startyear[k])){
      dates = seq(as.Date(paste(substr(tmplist3$starttime[k],1,4),substr(tmplist3$starttime[k],5,6),substr(tmplist3$starttime[k],7,8),sep="-")),
                  as.Date(paste(substr(tmplist3$endtime[k],1,4),substr(tmplist3$endtime[k],5,6),substr(tmplist3$endtime[k],7,8),sep="-")),
                  by="day")
      if(length(dates)>length(times)){
        
        if(length(which(substr(dates,6,10)=="02-29"))>0){
        dates = dates[-which(substr(dates,6,10)=="02-29")]
        }
        if(length(dates)>length(times)){
          startingdate = tmplist3$starttime[k]
          endingdate = tmplist3$endtime[k]
          seq5 <- seq(as.numeric(startingdate), as.numeric(endingdate))
          d5 <- seq5[ seq5 %/% 100 %% 100 %in% 1:12 & seq5 %% 100 %in% 1:30]
          dates = as.Date(as.character(d5),format="%Y%m%d")
        }
      }
    }
    
    if(k>1 & substr(dates[1],1,4)==substr(indexdate,1,4)){
      tunits = nctest$var[[varidx]]$dim[[3]]$units
      tunits = paste(substr(tunits,1,22),"12:00:00",sep=" ")
      idate = as.Date(substr(tunits,12,21))
      dates=idate+times
      times = dates-as.Date("1976-01-01")
    }
    
    if(any(is.na(dates)==TRUE)==TRUE){
      dates = dates[-which(is.na(dates)==TRUE)]
    }
    
    yearmon = unique(substr(dates,1,7))
    yearsinfile = as.numeric(substr(yearmon,1,4))
    if(any(yearsinfile<histperiod[1]) | any(yearsinfile>histperiod[length(histperiod)])){
      if(any(yearsinfile<histperiod[1])){
        yearmon = yearmon[-which(yearsinfile<histperiod[1])]
      }
    }
    
    
    timeoutmon= times[which(as.character(substr(dates,1,10)) %in% paste(yearmon,"01",sep="-"))]
    timeoutmonall = c(timeoutmonall,timeoutmon)
    
    vardata = ncvar_get(nctest,varname) #ncvar_get(nctest,nctest$var[[varidx]]$name) 
    
    if(mean(vardata,na.rm=TRUE)<1 & varname=="pr"){
      vardata = vardata*86400
    }
    nc_close(nctest)
    
    for(y in 1:length(yearmon)){
      
      outidx = which(substr(yearmonall,1,7)==yearmon[y])
      dateidx = which(substr(dates,1,7)==yearmon[y])
      
      if(varname=="pr"){
        outputarray1[,,outidx] = apply(vardata[,,dateidx],c(1,2),sum,na.rm=TRUE)
      }
      if(varname=="tasmax" | varname=="tasmin"){
        outputarray1[,,outidx] = apply(vardata[,,dateidx],c(1,2),mean,na.rm=TRUE)
      }
      #outputarray2[,,outidx] = apply(vardata[,,dateidx],c(1,2),calcrollsum,size=5)
      message("Year ",y," / ",length(yearmon))
    }
    message("Finished with file ",k," / ", length(filesuse))
  }
  
  #write historical files
  outyearidx = which(as.numeric(substr(yearmonall,1,4))>=histperiod[1] & as.numeric(substr(yearmonall,1,4))<=histperiod[2])
  
  outputarray1 = outputarray1[,,outyearidx]
  
  yearmonall = yearmonall[outyearidx]
  timeoutmonall = timeoutmonall[outyearidx]
  
  filenameout1 = paste("/home/woot0002/monthlyclimo/CMIP6/",varname,"_",tmplist3$model[1],"_",tmplist3$experiment[1],"_",histperiod[1],"-",histperiod[2],"_mon.nc",sep="")
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits, lat)
  dimT <- ncdim_def("time",timeunits,timeoutmonall,calendar=calendarinfo)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varname,unitsout,longname=varname, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout1 ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, outputarray1) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  
}
}
#####
# future calculations

if(runproj==TRUE){
projyears = projperiod[1]:projperiod[2]

for(i in 1:length(GCMs)){
  
  tmplist = filesproj[which(filesproj$model==GCMs[i]),]
  filesub = which(tmplist$endyear>=projperiod[1])
  tmplist2 = tmplist[filesub,]
  
  fileorder = order(tmplist2$endyear)
  tmplist2=tmplist2[fileorder,]
  
  filesuse = unique(tmplist2$path)
  
  tmplist3 = NULL
  for(j in 1:length(filesuse)){
    tmpidx = which(tmplist2$path==filesuse[j])
    tplist = tmplist2[tmpidx,]
    tmplist3 = rbind(tmplist3,tplist[1,])
  }
  
  if(any(tmplist3$startyear>projperiod[2])){
    filesout = which(tmplist3$startyear>projperiod[2])
    tmplist3 = tmplist3[-filesout,]
  }
  filesuse = unique(tmplist3$path)
  
  for(k in 1:length(filesuse)){
    
    nctest = nc_open(tmplist3$path[k]) 
    
    times = ncvar_get(nctest,"time")
    
    if(k==1){
      
      lon = ncvar_get(nctest,"lon")
      lat = ncvar_get(nctest,"lat")
      calendarinfo = nctest$var[[varname]]$dim[[3]]$calendar
      
      varnames = c()
      for(m in 1:length(nctest$var)){
        varnames[m] = nctest$var[[m]]$name
      }
      varidx = which(varnames==varname)
      timeunits = nctest$var[[varidx]]$dim[[3]]$units
      latunits = nctest$var[[varidx]]$dim[[2]]$units
      lonunits = nctest$var[[varidx]]$dim[[1]]$units
      
      timeunits = paste(substr(timeunits,1,22),"12:00:00",sep=" ")
      indexdate = substr(timeunits,12,21)
      origindate = as.PCICt(indexdate,cal=calendarinfo) # PCICt package helps with weird calendars
      
      timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
      
      #datesall = seq(as.Date(paste(tmplist3$startyear[which(tmplist3$startyear==min(tmplist3$startyear))],"-01-01",sep="")),as.Date(paste(tmplist3$endyear[which(tmplist3$endyear==max(tmplist3$endyear))],"-12-31",sep="")),by="day")
      yearmonall = seq(as.Date(paste(projyears[1],"-01-15",sep="")),as.Date(paste(projyears[length(projyears)],"-12-15",sep="")),by="month")
      timeoutmonall = c()
      
      outputarray1 = array(NA,dim=c(length(lon),length(lat),length(yearmonall)))
      
    }
    
    dates = origindate+times*86400
    
    if(as.numeric(substr(dates[1],1,4))<as.numeric(tmplist3$startyear[k])){
      message("Date correction applied")
      dates = seq(as.Date(paste(substr(tmplist3$starttime[k],1,4),substr(tmplist3$starttime[k],5,6),substr(tmplist3$starttime[k],7,8),sep="-")),
                  as.Date(paste(substr(tmplist3$endtime[k],1,4),substr(tmplist3$endtime[k],5,6),substr(tmplist3$endtime[k],7,8),sep="-")),
                  by="day")
      if(length(dates)>length(times)){
        if(length(which(substr(dates,6,10)=="02-29"))>0){
          dates = dates[-which(substr(dates,6,10)=="02-29")]
        }
        if(length(dates)>length(times)){
          startingdate = tmplist3$starttime[k]
          endingdate = tmplist3$endtime[k]
          seq5 <- seq(as.numeric(startingdate), as.numeric(endingdate))
          d5 <- seq5[ seq5 %/% 100 %% 100 %in% 1:12 & seq5 %% 100 %in% 1:30]
          dates = as.Date(as.character(d5),format="%Y%m%d")
        }
      }
    }
    
    if(k>1 & substr(dates[1],1,4)==substr(indexdate,1,4)){
      tunits = nctest$var[[varidx]]$dim[[3]]$units
      tunits = paste(substr(tunits,1,22),"12:00:00",sep=" ")
      idate = as.Date(substr(tunits,12,21))
      dates=idate+times
      times = dates-as.Date("1976-01-01")
    }
    
    if(any(is.na(dates)==TRUE)==TRUE){
      dates = dates[-which(is.na(dates)==TRUE)]
    }
    
    
    yearmon = unique(substr(dates,1,7))
    
    timeoutmon= times[which(as.character(substr(dates,1,10)) %in% paste(yearmon,"01",sep="-"))]
    timeoutmonall = c(timeoutmonall,timeoutmon)
    
    vardata = ncvar_get(nctest,varname) 
    
    if(mean(vardata,na.rm=TRUE)<1 & varname=="pr"){
      vardata = vardata*86400
    }
    nc_close(nctest)
    
    for(y in 1:length(yearmon)){
      
      outidx = which(substr(yearmonall,1,7)==yearmon[y])
      dateidx = which(substr(dates,1,7)==yearmon[y])
      if(varname=="pr"){
        outputarray1[,,outidx] = apply(vardata[,,dateidx],c(1,2),sum,na.rm=TRUE)
      }
      if(varname=="tasmax" | varname=="tasmin"){
        outputarray1[,,outidx] = apply(vardata[,,dateidx],c(1,2),mean,na.rm=TRUE)
      }
      message("Finished calcs for yearmon ",y," / ",length(yearmon))
    }
    message("Finished with file ",k," / ", length(filesuse))
  }
  
  #write historical files
  outyearidx = which(as.numeric(substr(yearmonall,1,4))>=projperiod[1] & as.numeric(substr(yearmonall,1,4))<=projperiod[2])
  
  outputarray1 = outputarray1[,,outyearidx]
  
  yearmonall = yearmonall[outyearidx]
  timeoutmonall = timeoutmonall[outyearidx]
  
  
  filenameout1 = paste("/home/woot0002/monthlyclimo/CMIP6/",varname,"_",tmplist3$model[1],"_",tmplist3$experiment[1],"_",rcpout,"_",projperiod[1],"-",projperiod[2],"_mon.nc",sep="")
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits, lat)
  dimT <- ncdim_def("time",timeunits,timeoutmonall,calendar=calendarinfo)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varname,unitsout,longname=varname, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout1 ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, outputarray1) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  
}
}
