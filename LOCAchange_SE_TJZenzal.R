source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(RPushbullet)

var = "FBI"
tempres = "ann"
period = c(2070,2099)

filesinhist = system(paste("ls /home/woot0002/TJZenzal/",var,"_*historical_",tempres,"_SE.nc",sep=""),intern=TRUE)
if(var!="frd" & var!="heatwaves" & var!="FLI" & var!="FBI"){
filesinrcp45 = system(paste("ls /home/woot0002/TJZenzal/",var,"_*rcp45_",tempres,"_",period[1],"-",period[2],"_SE.nc",sep=""),intern=TRUE)
filesinrcp85 = system(paste("ls /home/woot0002/TJZenzal/",var,"_*rcp85_",tempres,"_",period[1],"-",period[2],"_SE.nc",sep=""),intern=TRUE)
} else{
  filesinrcp45 = system(paste("ls /home/woot0002/TJZenzal/",var,"_*rcp45_",tempres,"_SE.nc",sep=""),intern=TRUE)
  filesinrcp85 = system(paste("ls /home/woot0002/TJZenzal/",var,"_*rcp85_",tempres,"_SE.nc",sep=""),intern=TRUE)
}
  
######
# climo calcs historical

GCMs_hist = c()

for(i in 1:length(filesinhist)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesinhist[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  GCMs_hist = c(GCMs_hist,filesplit2[2])
  
  nctest = nc_open(filesinhist[i])
  
  if(i==1){
    dataunits = nctest$var[[1]]$units
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    time=ncvar_get(nctest,"time")
    timeunits = nctest$var[[1]]$dim[[3]]$units
    latunits = nctest$var[[1]]$dim[[2]]$units
    lonunits = nctest$var[[1]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    
    if(tempres=="seas"){
      histdatout = array(NA,dim=c(length(lon),length(lat),4,length(filesinhist)))
    }
    if(tempres=="ann"){
      histdatout = array(NA,dim=c(length(lon),length(lat),length(filesinhist)))
    }
    
  }
  
  testdat = ncvar_get(nctest,var)
  if(var=="tasmax" | var=="tasmin"){
    if(testdat[140,140,1]>200){
      testdat = testdat-273.15
    }
  }
  
  if(mean(testdat,na.rm=TRUE)<1 & var=="pr"){
    testdat = testdat*86400
  }
  
  nc_close(nctest)
  
  if(tempres=="seas"){
    for(s in 1:4){
      if(s==1){
        idxout = which(as.numeric(substr(dates,6,7))==1)
      }
      if(s==2){
        idxout = which(as.numeric(substr(dates,6,7))==3)
      }
      if(s==3){
        idxout = which(as.numeric(substr(dates,6,7))==6)
      }
      if(s==4){
        idxout = which(as.numeric(substr(dates,6,7))==9)
      }
      tmp = apply(testdat[,,idxout],c(1,2),mean,na.rm=TRUE)
      histdatout[,,s,i]=ifelse(is.na(testdat[,,1])==FALSE,tmp,NA)
    }
  }
  if(tempres=="ann"){
    tmp = apply(testdat,c(1,2),mean,na.rm=TRUE)
    histdatout[,,i]=ifelse(is.na(testdat[,,1])==FALSE,tmp,NA)
  }
  
  message("finished with file ",i," / ",length(filesinhist))
  
}

######
# climo calcs rcp45

GCMs_rcp45 = c()

for(i in 1:length(filesinrcp45)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesinrcp45[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  GCMs_rcp45 = c(GCMs_rcp45,filesplit2[2])
  
  nctest = nc_open(filesinrcp45[i])
  
  if(i==1){
    dataunits = nctest$var[[1]]$units
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    time=ncvar_get(nctest,"time")
    timeunits = nctest$var[[1]]$dim[[3]]$units
    latunits = nctest$var[[1]]$dim[[2]]$units
    lonunits = nctest$var[[1]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    
    if(tempres=="seas"){
      rcp45datout = array(NA,dim=c(length(lon),length(lat),4,length(filesinrcp45)))
    }
    if(tempres=="ann"){
      rcp45datout = array(NA,dim=c(length(lon),length(lat),length(filesinrcp45)))
    }
    
  }
  
  testdat = ncvar_get(nctest,var)
  if(var=="tasmax" | var=="tasmin"){
    if(testdat[140,140,1]>200){
      testdat = testdat-273.15
    }
  }
  
  if(mean(testdat,na.rm=TRUE)<1 & var=="pr"){
    testdat = testdat*86400
  }
  
  nc_close(nctest)
  
  if(tempres=="seas"){
    for(s in 1:4){
      if(s==1){
        idxout = which(as.numeric(substr(dates,6,7))==1)
      }
      if(s==2){
        idxout = which(as.numeric(substr(dates,6,7))==3)
      }
      if(s==3){
        idxout = which(as.numeric(substr(dates,6,7))==6)
      }
      if(s==4){
        idxout = which(as.numeric(substr(dates,6,7))==9)
      }
      tmp = apply(testdat[,,idxout],c(1,2),mean,na.rm=TRUE)
      rcp45datout[,,s,i]=ifelse(is.na(testdat[,,1])==FALSE,tmp,NA)
    }
  }
  if(tempres=="ann"){
    tmp = apply(testdat,c(1,2),mean,na.rm=TRUE)
    rcp45datout[,,i]=ifelse(is.na(testdat[,,1])==FALSE,tmp,NA)
  }
  message("finished with file ",i," / ",length(filesinrcp45))
}


######
# climo calcs rcp85

GCMs_rcp85 = c()

for(i in 1:length(filesinrcp85)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesinrcp85[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  GCMs_rcp85 = c(GCMs_rcp85,filesplit2[2])
  
  nctest = nc_open(filesinrcp85[i])
  
  if(i==1){
    dataunits = nctest$var[[1]]$units
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    time=ncvar_get(nctest,"time")
    timeunits = nctest$var[[1]]$dim[[3]]$units
    latunits = nctest$var[[1]]$dim[[2]]$units
    lonunits = nctest$var[[1]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    
    if(tempres=="seas"){
      rcp85datout = array(NA,dim=c(length(lon),length(lat),4,length(filesinrcp85)))
    }
    if(tempres=="ann"){
      rcp85datout = array(NA,dim=c(length(lon),length(lat),length(filesinrcp85)))
    }
    
  }
  
  testdat = ncvar_get(nctest,var)
  if(var=="tasmax" | var=="tasmin"){
    if(testdat[140,140,1]>200){
      testdat = testdat-273.15
    }
  }
  
  if(mean(testdat,na.rm=TRUE)<1 & var=="pr"){
    testdat = testdat*86400
  }
  
  nc_close(nctest)
  
  if(tempres=="seas"){
    for(s in 1:4){
      if(s==1){
        idxout = which(as.numeric(substr(dates,6,7))==1)
      }
      if(s==2){
        idxout = which(as.numeric(substr(dates,6,7))==3)
      }
      if(s==3){
        idxout = which(as.numeric(substr(dates,6,7))==6)
      }
      if(s==4){
        idxout = which(as.numeric(substr(dates,6,7))==9)
      }
      tmp = apply(testdat[,,idxout],c(1,2),mean,na.rm=TRUE)
      rcp85datout[,,s,i]=ifelse(is.na(testdat[,,1])==FALSE,tmp,NA)
    }
  }
  if(tempres=="ann"){
      tmp = apply(testdat,c(1,2),mean,na.rm=TRUE)
      rcp85datout[,,i]=ifelse(is.na(testdat[,,1])==FALSE,tmp,NA)
    }
  
  message("finished with file ",i," / ",length(filesinrcp85))
}

######
# subset down if needed

if(length(GCMs_hist)!=length(GCMs_rcp45) | length(GCMs_hist)!=length(GCMs_rcp85) | length(GCMs_rcp45)!=length(GCMs_rcp85)){
  if(length(GCMs_hist)<length(GCMs_rcp45) & length(GCMs_hist)<length(GCMs_rcp85)){
    rcp45idx = which(GCMs_rcp45 %in% GCMs_hist)
    rcp85idx = which(GCMs_rcp85 %in% GCMs_hist)
    if(tempres=="seas"){
      rcp45datout = rcp45datout[,,,rcp45idx]
      rcp85datout = rcp85datout[,,,rcp85idx]
    } 
    if(tempres=="ann"){
      rcp45datout = rcp45datout[,,rcp45idx]
      rcp85datout = rcp85datout[,,rcp85idx]
    } 
    GCMs_rcp45 = GCMs_rcp45[rcp45idx]
    GCMs_rcp85 = GCMs_rcp45[rcp85idx]
    GCMnums = rcp45idx
  }
} else {
  GCMnums = 1:length(GCMs_hist)
}

#####
# projectedchange

rcp45change = rcp45datout-histdatout
rcp85change = rcp85datout-histdatout

#####
# ensemble means
if(tempres=="seas"){
  rcp45change_ensmean = apply(rcp45change,c(1,2,3),mean,na.rm=TRUE)  
  rcp85change_ensmean = apply(rcp85change,c(1,2,3),mean,na.rm=TRUE)
  rcp45datout_ensmean = apply(rcp45datout,c(1,2,3),mean,na.rm=TRUE)
  rcp85datout_ensmean = apply(rcp85datout,c(1,2,3),mean,na.rm=TRUE)
  histdatout_ensmean = apply(histdatout,c(1,2,3),mean,na.rm=TRUE)
}

if(tempres=="ann"){
  rcp45change_ensmean = apply(rcp45change,c(1,2),mean,na.rm=TRUE)  
  rcp85change_ensmean = apply(rcp85change,c(1,2),mean,na.rm=TRUE)
  rcp45datout_ensmean = apply(rcp45datout,c(1,2),mean,na.rm=TRUE)
  rcp85datout_ensmean = apply(rcp85datout,c(1,2),mean,na.rm=TRUE)
  histdatout_ensmean = apply(histdatout,c(1,2),mean,na.rm=TRUE)
}

#####
# netcdf writing

#####

if(tempres=="ann"){
  filenameout1 = paste("/home/woot0002/TJZenzal/",var,"_",period[1],"-",period[2],"_ann_climo_SE.nc",sep="")
  filenameout2 = paste("/home/woot0002/TJZenzal/",var,"_",period[1],"-",period[2],"_ann_climo_SE_ensmean.nc",sep="")
  filenameout3 = paste("/home/woot0002/TJZenzal/",var,"_",period[1],"-",period[2],"_ann_climo_SE_GCMs.csv",sep="")
  
  GCMtab= data.frame(GCMnums,GCMs_hist)
  names(GCMtab) = c("number","model")
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits, lat)
  dimG <- ncdim_def("model","GCM number",GCMnums)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1da <- ncvar_def(paste(var,"_ensmean_rcp45",sep=""),dataunits,longname=paste(var,"_ensmean_rcp45",sep=""), list(dimX,dimY), mv ,compression=9)
  var1db <- ncvar_def(paste(var,"_ensmean_rcp85",sep=""),dataunits,longname=paste(var,"_ensmean_rcp85",sep=""), list(dimX,dimY), mv ,compression=9)
  var2da <- ncvar_def(paste(var,"_ensmean_rcp45change",sep=""),dataunits,longname=paste(var,"_ensmean_rcp45change",sep=""), list(dimX,dimY), mv ,compression=9)
  var2db <- ncvar_def(paste(var,"_ensmean_rcp85change",sep=""),dataunits,longname=paste(var,"_ensmean_rcp85change",sep=""), list(dimX,dimY), mv ,compression=9)
  var3d <- ncvar_def(paste(var,"_ensmean_hist",sep=""),dataunits,longname=paste(var,"_ensmean_hist",sep=""), list(dimX,dimY), mv ,compression=9)
  
  var4da <- ncvar_def(paste(var,"_rcp45",sep=""),dataunits,longname=paste(var,"_allmodels_rcp45",sep=""), list(dimX,dimY,dimG), mv ,compression=9)
  var4db <- ncvar_def(paste(var,"_rcp85",sep=""),dataunits,longname=paste(var,"_allmodels_rcp85",sep=""), list(dimX,dimY,dimG), mv ,compression=9)
  var5da <- ncvar_def(paste(var,"_rcp45change",sep=""),dataunits,longname=paste(var,"_allmodels_rcp45change",sep=""), list(dimX,dimY,dimG), mv ,compression=9)
  var5db <- ncvar_def(paste(var,"_rcp85change",sep=""),dataunits,longname=paste(var,"_allmodels_rcp85change",sep=""), list(dimX,dimY,dimG), mv ,compression=9)
  var6d <- ncvar_def(paste(var,"_hist",sep=""),dataunits,longname=paste(var,"_allmodels_hist",sep=""), list(dimX,dimY,dimG), mv ,compression=9)
  
  
  #######
  # Create netcdf file var4da,var4db,var5da,var5db,var6d,
  
  nc <- nc_create(filenameout1 ,  list(var4da,var4db,var5da,var5db,var6d) )
  ncvar_put(nc, var4da, rcp45datout) # no start or count: write all values\
  ncvar_put(nc, var4db, rcp85datout) # no start or count: write all values\
  ncvar_put(nc, var5da, rcp45change) # no start or count: write all values\
  ncvar_put(nc, var5db, rcp85change) # no start or count: write all values\
  ncvar_put(nc, var6d, histdatout) # no start or count: write all values\
  nc_close(nc)
  
  nc2 <- nc_create(filenameout2 ,  list(var1da,var1db,var2da,var2db,var3d) )
  ncvar_put(nc2, var1da, rcp45datout_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var1db, rcp85datout_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var2da, rcp45change_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var2db, rcp85change_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var3d, histdatout_ensmean) # no start or count: write all values\
  nc_close(nc)
  
  write.table(GCMtab,file=filenameout3,row.names=FALSE,sep=",")
  
  gc()
  ptmend = proc.time()-ptm
  
}

if(tempres=="seas"){
  filenameout = paste("/home/woot0002/TJZenzal/",var,"_",period[1],"-",period[2],"_seas_climo_SE.nc",sep="")
  filenameout2 = paste("/home/woot0002/TJZenzal/",var,"_",period[1],"-",period[2],"_seas_climo_SE_ensmean.nc",sep="")
  filenameout3 = paste("/home/woot0002/TJZenzal/",var,"_",period[1],"-",period[2],"_seas_climo_SE_GCMs.csv",sep="")

  GCMtab= data.frame(GCMnums,GCMs_hist)
  names(GCMtab) = c("number","model")
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits, lat)
  dimT <- ncdim_def("seas","",1:4)
  dimG <- ncdim_def("model","GCM number",GCMnums)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1da <- ncvar_def(paste(var,"_ensmean_rcp45",sep=""),dataunits,longname=paste(var,"_ensmean_rcp45",sep=""), list(dimX,dimY,dimT), mv ,compression=9)
  var1db <- ncvar_def(paste(var,"_ensmean_rcp85",sep=""),dataunits,longname=paste(var,"_ensmean_rcp85",sep=""), list(dimX,dimY,dimT), mv ,compression=9)
  var2da <- ncvar_def(paste(var,"_ensmean_rcp45change",sep=""),dataunits,longname=paste(var,"_ensmean_rcp45change",sep=""), list(dimX,dimY,dimT), mv ,compression=9)
  var2db <- ncvar_def(paste(var,"_ensmean_rcp85change",sep=""),dataunits,longname=paste(var,"_ensmean_rcp85change",sep=""), list(dimX,dimY,dimT), mv ,compression=9)
  var3d <- ncvar_def(paste(var,"_ensmean_hist",sep=""),dataunits,longname=paste(var,"_ensmean_hist",sep=""), list(dimX,dimY,dimT), mv ,compression=9)
  
  var4da <- ncvar_def(paste(var,"_rcp45",sep=""),dataunits,longname=paste(var,"_allmodels_rcp45",sep=""), list(dimX,dimY,dimT,dimG), mv ,compression=9)
  var4db <- ncvar_def(paste(var,"_rcp85",sep=""),dataunits,longname=paste(var,"_allmodels_rcp85",sep=""), list(dimX,dimY,dimT,dimG), mv ,compression=9)
  var5da <- ncvar_def(paste(var,"_rcp45change",sep=""),dataunits,longname=paste(var,"_allmodels_rcp45change",sep=""), list(dimX,dimY,dimT,dimG), mv ,compression=9)
  var5db <- ncvar_def(paste(var,"_rcp85change",sep=""),dataunits,longname=paste(var,"_allmodels_rcp85change",sep=""), list(dimX,dimY,dimT,dimG), mv ,compression=9)
  var6d <- ncvar_def(paste(var,"_hist",sep=""),dataunits,longname=paste(var,"_allmodels_hist",sep=""), list(dimX,dimY,dimT,dimG), mv ,compression=9)
  
  
  #######
  # Create netcdf file var4da,var4db,var5da,var5db,var6d,
  
  nc <- nc_create(filenameout ,  list(var4da,var4db,var5da,var5db,var6d) )
  ncvar_put(nc, var4da, rcp45datout) # no start or count: write all values\
  ncvar_put(nc, var4db, rcp85datout) # no start or count: write all values\
  ncvar_put(nc, var5da, rcp45change) # no start or count: write all values\
  ncvar_put(nc, var5db, rcp85change) # no start or count: write all values\
  ncvar_put(nc, var6d, histdatout) # no start or count: write all values\
  nc_close(nc)
  
  nc2 <- nc_create(filenameout2 ,  list(var1da,var1db,var2da,var2db,var3d) )
  ncvar_put(nc2, var1da, rcp45datout_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var1db, rcp85datout_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var2da, rcp45change_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var2db, rcp85change_ensmean) # no start or count: write all values\
  ncvar_put(nc2, var3d, histdatout_ensmean) # no start or count: write all values\
  nc_close(nc)
  
  write.table(GCMtab,file=filenameout3,row.names=FALSE,sep=",")
  
  gc()
  ptmend = proc.time()-ptm
  
}








