#######################
#
# Monthly climo calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(PCICt)

var = "tasmin"

GCMtable = read.table("/home/woot0002/tasmin_GCMfilelist_LOCA.csv",sep=",",header=TRUE)
GCMtable$GCMs = paste(GCMtable$model,GCMtable$ensemble,sep="_")
GCMs = unique(paste(GCMtable$model,GCMtable$ensemble,sep="_"))

dates = seq(as.Date("2070-01-01"),as.Date("2099-12-31"),by="day")
monyear = dates[which(substr(dates,9,10)=="01")]
years = 2070:2099
months = 1:12

for(g in 1:length(GCMs)){
  
  tmptable = GCMtable[which(GCMtable$GCMs==GCMs[g] & GCMtable$variable==var),]
  
  for(y in 1:length(years)){
    tmptable$startyear = as.numeric(substr(tmptable$time,1,4))
    tmptable$endyear = as.numeric(substr(tmptable$time,10,13))
    
    if(nrow(tmptable)>1){
      periodcheck = c()
      for(i in 1:nrow(tmptable)){
        period = tmptable$startyear[i]:tmptable$endyear[i]
        periodcheck[i] = ifelse(years[y] %in% period,1,0)
      }
      tmptable2 = tmptable[which(periodcheck==1),]
    } else {
      tmptable2 = tmptable
    }
    
    if(nrow(tmptable2)<2){
    startdate = substr(tmptable2$time,1,8)
    enddate = substr(tmptable2$time,10,18)
    
    filedates = seq(as.Date(startdate,"%Y%m%d"),as.Date(enddate,"%Y%m%d"),by="day")
    
    test = nc_open(as.character(tmptable2$local_file[1]))
    time = ncvar_get(test,"time")
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    
    if(y==1){
      tmp = tmp2 = array(NA,dim=c(length(lon),length(lat),length(monyear)))
      if(var=="pr"){
        tmp3 = tmp4 = tmp5 = array(NA,dim=c(length(lon),length(lat),length(monyear)))
      }
    }
    
    if(length(time)<length(filedates)){
      filedates2 = filedates[-which(substr(filedates,6,10)=="02-29")]
      if(length(time)<length(filedates2)){
        filedates2=seq(as.PCICt(startdate,format="%Y%m%d",cal="360_day"),as.PCICt(enddate,format="%Y%m%d",cal="360_day"),by="day")
      }
    } else {
      filedates2 = filedates
    }
    
    yearidx = which(as.numeric(substr(filedates2,1,4))==years[y])
    
    vardata = ncvar_get(test,var,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    if(var=="pr"){
      vardata = vardata*86400
    }
    nc_close(test)
    tmpdates= filedates2[yearidx]
    } else { # if the dates are weird in some files.
      startdate = substr(tmptable2$time,1,8)
      enddate = substr(tmptable2$time,10,18)
      
      filedatesv1 = seq(as.Date(startdate[1],"%Y%m%d"),as.Date(enddate[1],"%Y%m%d"),by="day")
      filedatesv2 = seq(as.Date(startdate[2],"%Y%m%d"),as.Date(enddate[2],"%Y%m%d"),by="day")
      
      test = nc_open(as.character(tmptable2$local_file[1]))
      time = ncvar_get(test,"time")
      lon = ncvar_get(test,"lon")
      lat = ncvar_get(test,"lat")
      
      if(length(time)<length(filedatesv1)){
        filedates2v1 = filedatesv1[-which(substr(filedatesv1,6,10)=="02-29")]
        filedates2v2 = filedatesv2[-which(substr(filedatesv2,6,10)=="02-29")]
        if(length(time)<length(filedates2v1)){
          filedates2v1=seq(as.PCICt(startdate[1],format="%Y%m%d",cal="360_day"),as.PCICt(enddate[1],format="%Y%m%d",cal="360_day"),by="day")
          filedates2v2=seq(as.PCICt(startdate[2],format="%Y%m%d",cal="360_day"),as.PCICt(enddate[2],format="%Y%m%d",cal="360_day"),by="day")
        }
      } else {
        filedates2v1 = filedatesv1
        filedates2v2 = filedatesv2
      }
      nc_close(test)
      
      if(y==1){
        tmp = tmp2 = array(NA,dim=c(length(lon),length(lat),length(monyear)))
        if(var=="pr"){
          tmp3 = tmp4 = tmp5 = array(NA,dim=c(length(lon),length(lat),length(monyear)))
        }
      }
      
      if(filedates2v1[1]>filedates2v2[1]){
        filedates2 = c(filedates2v2,filedates2v1)
        v1 = as.character(tmptable2$local_file[2])
        f1 = filedates2v2
        yearidx1 = which(as.numeric(substr(f1,1,4))==years[y])
        
        v2 = as.character(tmptable2$local_file[1])
        f2 = filedates2v1
        yearidx2 = which(as.numeric(substr(f2,1,4))==years[y])
      }
      
      if(filedates2v1[1]<filedates2v2[1]){
        filedates2 = c(filedates2v1,filedates2v2)
        v1 = as.character(tmptable2$local_file[1])
        f1 = filedates2v1
        yearidx1 = which(as.numeric(substr(f1,1,4))==years[y])
        
        v2 = as.character(tmptable2$local_file[2])
        f2 = filedates2v2
        yearidx2 = which(as.numeric(substr(f2,1,4))==years[y])
      }
      
      yearidx = which(as.numeric(substr(filedates2,1,4))==years[y])
      tmpdates= filedates2[yearidx]
      
      vardata = array(NA,dim=c(length(lon),length(lat),length(tmpdates)))
      
      test=nc_open(v1)
      vardata1 = ncvar_get(test,var,start=c(1,1,yearidx1[1]),count=c(-1,-1,length(yearidx1)))
      nc_close(test)
      
      test=nc_open(v2)
      vardata2 = ncvar_get(test,var,start=c(1,1,yearidx2[1]),count=c(-1,-1,length(yearidx2)))
      nc_close(test)
      
      vardata[,,1:length(yearidx1)]=vardata1
      vardata[,,(length(yearidx1)+1):length(tmpdates)]=vardata2
      rm(vardata1)
      rm(vardata2)
      
      if(var=="pr"){
        vardata = vardata*86400
      }
    }
    
    for(m in 1:12){
      monidx = m +12*(y-1)
      mondatesidx = which(as.numeric(substr(tmpdates,6,7))==m)
      if(var=="pr"){
        tmp[,,monidx] = apply(vardata[,,mondatesidx],c(1,2),sum,na.rm=TRUE)
      pr50 = ifelse(vardata>=50.4,1,0)
      r1mm = ifelse(vardata>=1,1,0)
      tmp2[,,monidx] = apply(pr50[,,mondatesidx],c(1,2),sum,na.rm=TRUE)
      tmp3[,,monidx] = apply(r1mm[,,mondatesidx],c(1,2),sum,na.rm=TRUE)
      tmp4[,,monidx] = apply(vardata[,,mondatesidx],c(1,2),max,na.rm=TRUE)
      tmp5[,,monidx] = apply(vardata[,,mondatesidx],c(1,2),calcrollsum,size=5)
      }
      if(var=="tasmax"){
        tmp[,,monidx] = apply(vardata[,,mondatesidx]-273.15,c(1,2),mean,na.rm=TRUE)
        tmax95 = ifelse((vardata-273.15)>=35,1,0)
        tmp2[,,monidx] = apply(tmax95[,,mondatesidx],c(1,2),sum,na.rm=TRUE)
      }
      if(var=="tasmin"){
        tmp[,,monidx] = apply(vardata[,,mondatesidx]-273.15,c(1,2),mean,na.rm=TRUE)
        tmin32 = ifelse((vardata-273.15)<=0,1,0)
        tmp2[,,monidx] = apply(tmin32[,,mondatesidx],c(1,2),sum,na.rm=TRUE)
      }
    }
    message("Finished calcs for year ",y," / ",length(years))
  }

  var1mon = var2mon = array(NA,dim=c(length(lon),length(lat),12))
  if(var=="pr"){
    var3mon = var4mon = var5mon = array(NA,dim=c(length(lon),length(lat),12))
  } 
  
  for(m in 1:12){
    idx = which(as.numeric(substr(monyear,6,7))==m)
    var1mon[,,m] = apply(tmp[,,idx],c(1,2),mean,na.rm=TRUE)
    var2mon[,,m] = apply(tmp2[,,idx],c(1,2),mean,na.rm=TRUE)
    if(var=="pr"){
      var3mon[,,m] = apply(tmp3[,,idx],c(1,2),mean,na.rm=TRUE)
      var4mon[,,m] = apply(tmp4[,,idx],c(1,2),mean,na.rm=TRUE)
      var5mon[,,m] = apply(tmp5[,,idx],c(1,2),mean,na.rm=TRUE)
    }
  }
  
  newfilename = paste(var,GCMs[g],"CMIP5_projclimo_mon.nc",sep="_")
  
  dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  dimT <- ncdim_def("time","month",1:12)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  
  if(var=="pr"){
    units1="mm"
    lname1 = "Monthly Precipitation Climo"
    varname1 = "prclimo"
    units2="days"
    lname2 = "Monthly Heavy Rain Days (>=50.4mm) Climo"
    varname2 = "pr50climo"
    units3="days"
    lname3 = "Monthly Rain Days (>=1mm) Climo"
    varname3 = "r1mmclimo"
    units4="mm"
    lname4 = "Monthly Max 1-day Rainfall Climo"
    varname4 = "rx1dayclimo"
    units5="mm"
    lname5 = "Monthly Max 5-day Rainfall Climo"
    varname5 = "rx5dayclimo"
  } 
  
  if(var=="tasmax"){
    units1="C"
    lname1 = "Monthly High Temperature Climo"
    varname1 = "tmaxclimo"
    units2="days"
    lname2 = "Monthly Hot Days (tmax>=95F) Climo"
    varname2 = "tmax95climo"
  }
  
  if(var=="tasmin"){
    units1="C"
    lname1 = "Monthly Low Temperature Climo"
    varname1 = "tminclimo"
    units2="days"
    lname2 = "Monthly Cold Days (tmin<=32F) Climo"
    varname2 = "tmin32climo"
  }
  
  var1d <- ncvar_def(varname1,units=units1,longname=lname1, list(dimX,dimY,dimT), mv )
  var2d <- ncvar_def(varname2,units=units2,longname=lname2, list(dimX,dimY,dimT), mv )
  if(var=="pr"){
    var3d <- ncvar_def(varname3,units=units3,longname=lname3, list(dimX,dimY,dimT), mv )
    var4d <- ncvar_def(varname4,units=units4,longname=lname4, list(dimX,dimY,dimT), mv )
    var5d <- ncvar_def(varname5,units=units5,longname=lname5, list(dimX,dimY,dimT), mv )
  }
  
  #######
  # Create netcdf file
  if(var!="pr"){
    nc <- nc_create(paste("/home/woot0002/GCMs/",newfilename,sep="") ,  list(var1d,var2d) )
    # Write some data to the file
    ncvar_put(nc, var1d, var1mon) # no start or count: write all values\
    ncvar_put(nc, var2d, var2mon) # no start or count: write all values\
    # close ncdf
    nc_close(nc)
  } else {
    nc <- nc_create(paste("/home/woot0002/GCMs/",newfilename,sep="") ,  list(var1d,var2d,var3d,var4d,var5d) )
    # Write some data to the file
    ncvar_put(nc, var1d, var1mon) # no start or count: write all values\
    ncvar_put(nc, var2d, var2mon) # no start or count: write all values\
    ncvar_put(nc, var3d, var3mon) # no start or count: write all values\
    ncvar_put(nc, var4d, var4mon) # no start or count: write all values\
    ncvar_put(nc, var5d, var5mon) # no start or count: write all values\
    # close ncdf
    nc_close(nc)
  }
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(files))
  message("Process time: ",ptmend[3]," secs")
  
}

