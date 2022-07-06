source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(RPushbullet)

var = "pr"
tempres = "monthly"
varout = "pr"

futureperiod = c(2070,2099)
scen = "rcp85"

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",var,"/future/*",scen,"*.nc",sep=""),intern=TRUE)
idxout=grep("2100",filesin)
if(length(idxout)>0){
  filesin = filesin[-idxout]
}

if(var=="tasmax" | var=="tasmin"){
filesplit = do.call(rbind,strsplit(filesin,"_",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,ncol(filesplit)],".",fixed=TRUE))
filesplit3 = cbind(filesplit[,c(3,4,5)])
filesplit3 = data.frame(filesplit3)
names(filesplit3)=c("GCM","exp","scen")
year = as.numeric(filesplit2[,1])
filesplit3 = cbind(filesplit3,year)

idxin = which(filesplit3$year>=futureperiod[1] & filesplit3$year<=futureperiod[2])
filesin = filesin[idxin]

}

if(var=="pr"){
  filesplit = do.call(rbind,strsplit(filesin,"_",fixed=TRUE))
  filesplit2 = do.call(rbind,strsplit(filesplit[,ncol(filesplit)],".",fixed=TRUE))
  filesplit2a = do.call(rbind,strsplit(filesplit2[,1],"-",fixed=TRUE))
  
  filesplit3 = cbind(filesplit[,c(3,4)])
  filesplit3 = data.frame(filesplit3)
  names(filesplit3)=c("GCM","exp")
  
  filesplit2a = data.frame(filesplit2a)
  filesplit2a[,2] = as.numeric(as.character(filesplit2a[,2]))
  filesplit2a[,3] = as.numeric(as.character(filesplit2a[,3]))
  names(filesplit2a) = c("scen","startyear","endyear")
  
  filesplit3 = cbind(filesplit3,filesplit2a)
  
  filesplit3$idx = NA
  for(j in 1:nrow(filesplit3)){
    filesplit3$idx[j] = ifelse(any(filesplit3$startyear[j]:filesplit3$endyear[j] %in% futureperiod[1]:futureperiod[2])==TRUE,1,0)
  }
  
  idxin = which(filesplit3$idx==1)
  filesin = filesin[idxin]
  filesplit3 = filesplit3[idxin,]
}

GCMs = as.character(unique(filesplit3$GCM))

latbnds = c(24,42)
lonbnds = c(-112,-72)

for(i in 1:length(GCMs)){
  
  ptm = proc.time()
  
  filesused = filesin[grep(paste(GCMs[i],"_",sep=""),filesin)]
  timeout=c()
  datesout=c()
  for(f in 1:length(filesused)){
    
    filesplit = strsplit(filesused[f],split="/",fixed=TRUE)
    filesplit2 = strsplit(filesplit[[1]][[length(filesplit[[1]])]],"_",fixed=TRUE)
    filesplit2 = filesplit2[[length(filesplit2)]]
    
    varin = paste(filesplit2[1],filesplit2[2],filesplit2[3],scen,sep="_")
    varout2 = paste(varout,filesplit2[2],filesplit2[3],scen,sep="_")
    
    nctest = nc_open(filesused[f])
    
    if(i==1){
      dataunits = nctest$var[[1]]$units
      lon = ncvar_get(nctest,"lon")
      lat = ncvar_get(nctest,"lat")
      timeunits = nctest$var[[1]]$dim[[3]]$units
      latunits = nctest$var[[1]]$dim[[2]]$units
      lonunits = nctest$var[[1]]$dim[[1]]$units
      
      timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
      
      LON = rep(lon,each=length(lat))
      LAT = rep(lat,length(lon))
      R = rep(1:length(lon),each=length(lat))
      C = rep(1:length(lat),length(lon))
      modelgrid = data.frame(R,C,LON,LAT)
      names(modelgrid) = c("R","C","lon","lat")
      if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
      
      idxuse = which(modelgrid$lat>=latbnds[1] & modelgrid$lat<=latbnds[2] & modelgrid$lon>=lonbnds[1] & modelgrid$lon<=lonbnds[2])
      modelgriduse = modelgrid[idxuse,]
      
      lonstart = min(modelgriduse$R)
      latstart = min(modelgriduse$C)
      loncount = length(unique(modelgriduse$R))
      latcount = length(unique(modelgriduse$C))
      
      lonout = unique(modelgriduse$lon)
      latout = unique(modelgriduse$lat)
      
    }
    time=ncvar_get(nctest,"time")
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    
    datesout = c(datesout,as.character(dates))
    timeout= c(timeout,time[which(substr(dates,9,10)=="01")])
    
    if(f==1){
      yearmon = seq(as.Date(paste(futureperiod[1],"-01-01",sep="")),as.Date(paste(futureperiod[2],"-12-01",sep="")),by="month")
      yearmon = substr(yearmon,1,7)
      testdatout = array(NA,dim=c(length(lonout),length(latout),length(yearmon)))
    }
    testdat = ncvar_get(nctest,varin,start=c(lonstart,latstart,1),count=c(loncount,latcount,-1))
    if(mean(testdat,na.rm=TRUE)<1 & var=="pr"){
      testdat=testdat*86400
    }
    nc_close(nctest)
    
    if(var=="tasmin" | var=="tasmax"){
      if(f==1){
        I = 1
      }
      ymon = unique(substr(dates,1,7))
      for(y in 1:length(ymon)){
        idxuse = which(substr(dates,1,7)==ymon[y])
        if(varout=="tasmax" | varout=="tasmin"){
          tmp1 = apply(testdat[,,idxuse],c(1,2),mean,na.rm=TRUE)
        }
        if(varout=="tmax100"){
          td = ifelse(testdat[,,idxuse]>=37.7778,1,0)
          tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
        }
        
        testdatout[,,I]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
        I=I+1
      }
    }
    
    if(var=="pr"){
      if(f==1){
        I = 1
      }
      dates = dates[-length(dates)]
      ymon = unique(substr(dates,1,7))
      ymon = ymon[which(as.numeric(substr(ymon,1,4))>=futureperiod[1] & as.numeric(substr(ymon,1,4))<=futureperiod[2])]
      
      for(y in 1:length(ymon)){
        idxuse = which(substr(dates,1,7)==ymon[y])
        if(varout=="pr"){
          tmp1 = apply(testdat[,,idxuse],c(1,2),sum,na.rm=TRUE)
        }
        testdatout[,,I]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
        I=I+1
      }
    }
    message("finished with file ",f," / ",length(filesused))
  }
  
  if(tempres=="daily"){
    filenameout = paste("/home/woot0002/RCMES/data/LOCA/",varout2,"_SE.nc",sep="")
  }
  if(tempres=="monthly"){
    filenameout = paste("/home/woot0002/RCMES/data/LOCA/",varout2,"_mon_",futureperiod[1],"-",futureperiod[2],"_SE.nc",sep="")
  }
  
  if(var=="pr"){
    datesoutused = datesout[which(substr(datesout,9,10)=="01")]
    timeout = timeout[which(as.numeric(substr(datesoutused,1,4))>=futureperiod[1] & as.numeric(substr(datesoutused,1,4))<=futureperiod[2])]
    
  }
  
  dimX <- ncdim_def( "lon", lonunits, lonout)
  dimY <- ncdim_def( "lat", latunits, latout)
  dimT <- ncdim_def("time",timeunits,timeout)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varout,dataunits,longname=varout, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, testdatout) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(testdat)
  gc()
  ptmend = proc.time()-ptm
  
  pbPost(type="note",title="Future LOCA trimmed file written out",body=paste("Trimmed file written for ",tempres," ",varin," it took ",ptmend[3]," secs.",sep=""),email="amwootte@ou.edu")
  
}

  
 

