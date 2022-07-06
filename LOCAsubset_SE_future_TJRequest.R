source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(multiApply)
#library(RPushbullet)

var = "tasmax"
tempres = "annual"
varout = "FLI"

futureperiod = c(2070,2099)
scen = "rcp85"

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",var,"/future/*",scen,"*.nc",sep=""),intern=TRUE)
idxout=grep("2006-2100",filesin)
if(length(idxout)>0){
  filesin = filesin[-idxout]
}

if(varout=="heatwaves"){
  filesin = filesin[-grep("HadGEM2-ES",filesin)]
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

latbnds = c(25,37)
lonbnds = c(-107.90,-88.770)

for(i in 1:length(GCMs)){
  
  ptm = proc.time()
  
  filesused = filesin[grep(paste(GCMs[i],"_",sep=""),filesin)]
  timeout=c()
  datesout=c()
  for(f in 1:length(filesused)){ #
    
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
    
    timeout= c(timeout,time)
    
    if(f==1){
      
      if(var=="pr"){
        filesplittab = do.call("rbind",strsplit(filesused,split="-",fixed=TRUE))
        datesall = seq(as.Date(paste(filesplittab[1,(ncol(filesplittab)-1)],"-01-01",sep="")),as.Date(paste(substr(filesplittab[nrow(filesplittab),ncol(filesplittab)],1,4),"-12-31",sep="")),by="day")
      }
      if(var=="tasmax"){
        filesplittab = do.call("rbind",strsplit(filesused,split="_",fixed=TRUE))
        datesall = seq(as.Date(paste(substr(filesplittab[1,ncol(filesplittab)],1,4),"-01-01",sep="")),as.Date(paste(substr(filesplittab[nrow(filesplittab),ncol(filesplittab)],1,4),"-12-31",sep="")),by="day")
      }
      if(var=="tasmin"){
        filesplittab = do.call("rbind",strsplit(filesused,split="_",fixed=TRUE))
        datesall = seq(as.Date(paste(substr(filesplittab[1,ncol(filesplittab)],1,4),"-01-01",sep="")),as.Date(paste(substr(filesplittab[nrow(filesplittab),ncol(filesplittab)],1,4),"-12-31",sep="")),by="day")
      }
      testdatout = array(NA,dim=c(length(lonout),length(latout),length(datesall)))
      if(varout=="heatwaves"){
        testdatout2 = array(NA,dim=c(length(lonout),length(latout),length(datesall)))
      }
    }
    testdat = ncvar_get(nctest,varin,start=c(lonstart,latstart,1),count=c(loncount,latcount,-1))
    if(mean(testdat,na.rm=TRUE)<1 & var=="pr"){
      testdat=testdat*86400
    }
    if(var=="tasmax" | var=="tasmin"){
      if(testdat[140,140,1]>200){
        testdat = testdat-273.15
      }
    }
    nc_close(nctest)
    
    if(varout=="heatwaves"){
    file2 = paste("/data4/data/DS_proj/LOCA/tasmin/future/tasmin",filesplittab[f,3],filesplittab[f,4],filesplittab[f,5],filesplittab[f,6],sep="_")
    varin2 = paste("tasmin",filesplit2[2],filesplit2[3],filesplit2[4],sep="_")
    
    nctest = nc_open(file2)
    testdat2 = ncvar_get(nctest,varin2,start=c(lonstart,latstart,1),count=c(loncount,latcount,-1))
    nc_close(nctest)
    if(testdat2[140,140,1]>200){
      testdat2 = testdat2-273.15
    }
    
    if(f==1){
      
      tmaxq95file = system(paste("ls /home/woot0002/TJZenzal/tmaxq95",filesplittab[f,3],"*_historical_climo_SE.nc",sep="_"),intern=TRUE)
      tminq95file = system(paste("ls /home/woot0002/TJZenzal/tminq95",filesplittab[f,3],"*_historical_climo_SE.nc",sep="_"),intern=TRUE)
    
      nctest = nc_open(tmaxq95file)
      tmaxq95 = ncvar_get(nctest,"tmaxq95")
      nc_close(nctest)
    
      nctest = nc_open(tminq95file)
      tminq95 = ncvar_get(nctest,"tminq95")
      nc_close(nctest)
    }
    }
    
    if(var=="pr"){
      dateidx = which(as.numeric(substr(datesall,1,4))>=as.numeric(filesplittab[f,(ncol(filesplittab)-1)]) & as.numeric(substr(datesall,1,4))<=as.numeric(substr(filesplittab[f,ncol(filesplittab)],1,4)))
    }
    if(var=="tasmax" | var=="tasmin"){
      dateidx = which(as.numeric(substr(datesall,1,4))==as.numeric(substr(filesplittab[f,ncol(filesplittab)],1,4)))
    }
    
    if((length(dateidx)-dim(testdat)[3])==1){
      dateidx = dateidx[-length(dateidx)]
    }
    if((length(dateidx)-dim(testdat)[3])== -1){
      dateidx = c(dateidx[1]-1,dateidx)
    }
    
    testdatout[,,dateidx]= testdat
    if(varout=="heatwaves"){
      testdatout2[,,dateidx]=testdat2
    }
    rm(testdat)
    message("finished with file ",f," / ",length(filesused))
    
}
  
  datesidx2 = which(as.numeric(substr(datesall,1,4))>=futureperiod[1] & as.numeric(substr(datesall,1,4))<=futureperiod[2])
  testdatout = testdatout[,,datesidx2]
  datesout = datesall[datesidx2]
  
  if(varout=="heatwaves"){
    testdatout2 = testdatout2[,,datesidx2]
  }
  
  if(tempres=="annual"){
    years = unique(substr(datesout,1,4))
    tmp = array(NA,dim=c(loncount,latcount,length(years)))
    for(y in 1:length(years)){
      idxmon = which(substr(datesout,1,4)==years[y])
      
      if(varout=="pr"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),sum,na.rm=TRUE)  
      } 
      if(varout=="rx1day"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),max,na.rm=TRUE)
      }
      if(varout=="rx5day"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),calcrollsum,size=5)
      }
      if(varout=="tasmax" | varout=="tasmin"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="frd"){
        tmp1 = array(NA,dim=c(loncount,latcount))
        for(r in 1:lonout){
          for(c in 1:latout){
            message("Working on r ",r," and c ",c)
            tmp1[r,c] = lastfreeze(testdatout[r,c,idxmon],startdate=startdate,inputtimes=timeout[idxmon],threshold=0)
          }
        }
      }
      if(varout=="FLI" | varout=="FBI"){
        latmat = matrix(NA,nrow=length(lonout),ncol=length(latout))
        for(r in 1:nrow(latmat)){
          latmat[r,] = latout
        }
        data=list(TMAX = testdatout[,,idxmon],TMIN=testdatout2[,,idxmon],lat=latmat)
        if(varout=="FLI"){
          tmp1 <- Apply(data, target = list(3, 3, NULL), findfli)$output1  
        }
        if(varout=="FBI"){
          tmp1 <- Apply(data, target = list(3, 3, NULL), findfbi)$output1  
        }
      }
      
      if(varout=="heatwaves"){
        tmp1 = array(NA,dim=c(loncount,latcount))
        for(r in 1:length(lonout)){
          for(c in 1:length(latout)){
            tmp1[r,c] = heatwave.calc(testdatout[r,c,idxmon],testdatout2[r,c,idxmon],tmaxq95[r,c],tminq95[r,c])
            #message("Calcs complete for r ",r," and c ",c)
          }
        }
      }
      if(varout=="tmax100"){
        td = ifelse(testdatout[,,idxmon]>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      
      tmp[,,y]=ifelse(is.na(testdatout[,,1])==TRUE,NA,tmp1)
      message("Finished with yearmon ",y," / ",length(years))
    }
    testdat=tmp
  }
  
  if(tempres=="monthly"){
    yearmon = unique(substr(datesout,1,7))
    tmp = array(NA,dim=c(loncount,latcount,length(yearmon)))
    for(y in 1:length(yearmon)){
      idxmon = which(substr(datesout,1,7)==yearmon[y])
      
      if(varout=="pr"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),sum,na.rm=TRUE)  
      } 
      if(varout=="rx1day"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),max,na.rm=TRUE)
      }
      if(varout=="rx5day"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),calcrollsum,size=5)
      }
      if(varout=="tasmax" | varout=="tasmin"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="tmax100"){
        td = ifelse(testdatout[,,idxmon]>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      
      tmp[,,y]=ifelse(is.na(testdatout[,,1])==TRUE,NA,tmp1)
      message("Finished with yearmon ",y," / ",length(yearmon))
    }
    testdat=tmp
  }
  
  if(tempres=="seasonal"){
    years = as.numeric(substr(datesout,1,4))
    seas = rep(NA,length(datesout))
    seas[which(as.numeric(substr(datesout,6,7))>=12 | as.numeric(substr(datesout,6,7))<=2)] = 1
    seas[which(as.numeric(substr(datesout,6,7))>=3 & as.numeric(substr(datesout,6,7))<=5)] = 2
    seas[which(as.numeric(substr(datesout,6,7))>=6 & as.numeric(substr(datesout,6,7))<=8)] = 3
    seas[which(as.numeric(substr(datesout,6,7))>=9 & as.numeric(substr(datesout,6,7))<=11)] = 4
    
    yearseasall = paste(years,seas,sep="-")
    yearseas = unique(yearseasall)
    
    tmp = array(NA,dim=c(loncount,latcount,length(yearseas)))
    for(y in 1:length(yearseas)){
      idxmon = which(yearseasall==yearseas[y])
      
      if(varout=="pr"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),sum,na.rm=TRUE)  
      } 
      if(varout=="rx1day"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),max,na.rm=TRUE)
      }
      if(varout=="rx5day"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),calcrollsum,size=5)
      }
      if(varout=="mdrn"){
        tp = ifelse(testdatout[,,idxmon]>=0.254,1,0)
        tmp1 = apply(tp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="tasmax" | varout=="tasmin"){
        tmp1 = apply(testdatout[,,idxmon],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="tmax100"){
        td = ifelse(testdatout[,,idxmon]>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      
      tmp[,,y]=ifelse(is.na(testdatout[,,1])==TRUE,NA,tmp1)
      message("Finished with yearmon ",y," / ",length(yearseas))
    }
    testdat=tmp
  }
  
  if(tempres=="daily"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_day_",futureperiod[1],"-",futureperiod[2],"_SE.nc",sep="")
  }
  if(tempres=="monthly"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_mon_",futureperiod[1],"-",futureperiod[2],"_SE.nc",sep="")
  }
  if(tempres=="annual"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_ann_SE.nc",sep="")
    timeout = timeout[which(substr(datesout,6,10)=="01-01")]
  }
  if(tempres=="seasonal"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_seas_",futureperiod[1],"-",futureperiod[2],"_SE.nc",sep="")
    
    if(i==1){
      whichlist = c()
      yearsu = unique(years)
      monthsin = c(1,3,6,9)
      I=1
      for(j in 1:length(yearsu)){
        for(m in 1:length(monthsin)){
          whichlist[I] = which(as.numeric(substr(datesout,1,4))==yearsu[j] & as.numeric(substr(datesout,6,7))==monthsin[m] & as.numeric(substr(datesout,9,10))==1)
          I=I+1
        }
      }
    }
    timeout = timeout[datesidx2]
    timeout = timeout[whichlist]
  }
  
  #if(var=="pr"){
  #  datesoutused = datesout[which(substr(datesout,9,10)=="01")]
  #  timeout = timeout[which(as.numeric(substr(datesoutused,1,4))>=futureperiod[1] & as.numeric(substr(datesoutused,1,4))<=futureperiod[2])]
  #}
  
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
  ncvar_put(nc, var1d, testdat) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(testdat)
  gc()
  ptmend = proc.time()-ptm
  
  #pbPost(type="note",title="Future LOCA trimmed file written out",body=paste("Trimmed file written for ",tempres," ",varin," it took ",ptmend[3]," secs.",sep=""),email="amwootte@ou.edu")
}

  

 

