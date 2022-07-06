source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(ggplot2)

var = "tasmax"

BINLIMIT=20
colorchoicediff = "bluetored"
diffbartype = "difference"
period = "sea"

histfiles = system(paste("ls /home/woot0002/RCMES/data/LOCA/",var,"*historical_mon_SE.nc",sep=""),intern=TRUE)
rcp45files = system(paste("ls /home/woot0002/RCMES/data/LOCA/",var,"*rcp45_mon_*_SE.nc",sep=""),intern=TRUE)
rcp85files = system(paste("ls /home/woot0002/RCMES/data/LOCA/",var,"*rcp85_mon_*_SE.nc",sep=""),intern=TRUE)

####
# make sure equal number of GCMs historical and future.

filesplit = do.call(rbind,strsplit(histfiles,split="/",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,ncol(filesplit)],split="_",fixed=TRUE))
filesplit2 = data.frame(filesplit2[,1:4])
names(filesplit2)=c("var","GCM","exp","scen")
histfilelist = filesplit2

filesplit = do.call(rbind,strsplit(rcp45files,split="/",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,ncol(filesplit)],split="_",fixed=TRUE))
filesplit2a = do.call(rbind,strsplit(filesplit2[,6],split="-",fixed=TRUE))
startyear = as.numeric(filesplit2a[,1])
endyear = as.numeric(filesplit2a[,2])
filesplit2 = data.frame(filesplit2[,1:4],startyear,endyear)
names(filesplit2)=c("var","GCM","exp","scen","startyear","endyear")
rcp45filelist = filesplit2

filesplit = do.call(rbind,strsplit(rcp85files,split="/",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,ncol(filesplit)],split="_",fixed=TRUE))
filesplit2a = do.call(rbind,strsplit(filesplit2[,6],split="-",fixed=TRUE))
startyear = as.numeric(filesplit2a[,1])
endyear = as.numeric(filesplit2a[,2])
filesplit2 = data.frame(filesplit2[,1:4],startyear,endyear)
names(filesplit2)=c("var","GCM","exp","scen","startyear","endyear")
rcp85filelist = filesplit2

if(nrow(histfilelist)<nrow(rcp85filelist) | nrow(histfilelist)<nrow(rcp85filelist)){
  rcp85idx = which(rcp85filelist$GCM %in% histfilelist$GCM)
  rcp85files = rcp85files[rcp85idx]
  rcp85filelist=rcp85filelist[rcp85idx,]
  rcp45idx = which(rcp45filelist$GCM %in% histfilelist$GCM)
  rcp45files = rcp45files[rcp45idx]
  rcp45filelist=rcp45filelist[rcp45idx,]
}

if(nrow(histfilelist)>nrow(rcp85filelist) | nrow(histfilelist)>nrow(rcp85filelist)){
  if(nrow(rcp85filelist)<nrow(rcp45filelist)){
    rcp45idx = which(rcp45filelist$GCM %in% rcp85filelist$GCM)
    rcp45files = rcp45files[rcp45idx]
    rcp45filelist=rcp45filelist[rcp45idx,]
    histidx = which(histfilelist$GCM %in% rcp85filelist$GCM)
    histfiles = histfiles[histidx]
    histfilelist=histfilelist[histidx,]
  }
  if(nrow(rcp85filelist)>nrow(rcp45filelist)){
    rcp85idx = which(rcp85filelist$GCM %in% rcp45filelist$GCM)
    rcp85files = rcp85files[rcp85idx]
    rcp85filelist=rcp85filelist[rcp85idx,]
    histidx = which(histfilelist$GCM %in% rcp45filelist$GCM)
    histfiles = histfiles[histidx]
    histfilelist=histfilelist[histidx,]
  }
}

################
# Determine historical climatology

histdates = seq(as.Date("1976-01-15"),as.Date("2005-12-15"),by="month")
years = unique(as.numeric(substr(histdates,1,4)))

for(i in 1:length(histfiles)){
  
  nctest = nc_open(histfiles[i])
  
  if(i==1){
    lon=ncvar_get(nctest,"lon")
    lat=ncvar_get(nctest,"lat")
    if(period=="ann"){
      histclimo = array(NA,dim=c(length(lon),length(lat),length(histfiles)))
    }
    if(period=="sea"){
      histclimo = array(NA,dim=c(length(lon),length(lat),length(histfiles),4))
    }
  }
  
  vardat = ncvar_get(nctest,var)
  nc_close(nctest)
  
  if(period=="ann"){
    tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      yearidx = which(as.numeric(substr(histdates,1,4))==years[y])
      if(var=="tasmax" | var=="tasmin"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),mean,na.rm=TRUE)  
      }
      if(var=="pr" | var=="tmax100"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),sum,na.rm=TRUE)  
      }
      if(var=="rx1day"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),max,na.rm=TRUE)  
      }
    }
    tmp2 = apply(tmp,c(1,2),mean,na.rm=TRUE)
    tmp2 = ifelse(is.na(vardat[,,1])==FALSE,tmp2,NA)
    histclimo[,,i] = tmp2
  }
  
  if(period=="sea"){
    for(s in 1:4){
    tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      if(s==1){
        yearidx = which(as.numeric(substr(histdates,1,4))==years[y] & (as.numeric(substr(histdates,6,7))<=2 | as.numeric(substr(histdates,6,7))==12))
      }
      if(s==2){
        yearidx = which(as.numeric(substr(histdates,1,4))==years[y] & as.numeric(substr(histdates,6,7))>=3 & as.numeric(substr(histdates,6,7))<=5)
      }
      if(s==3){
        yearidx = which(as.numeric(substr(histdates,1,4))==years[y] & as.numeric(substr(histdates,6,7))>=6 & as.numeric(substr(histdates,6,7))<=8)
      }
      if(s==4){
        yearidx = which(as.numeric(substr(histdates,1,4))==years[y] & as.numeric(substr(histdates,6,7))>-9 & as.numeric(substr(histdates,6,7))<=11)
      }
      if(var=="tasmax" | var=="tasmin"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),mean,na.rm=TRUE)  
      }
      if(var=="pr" | var=="tmax100"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),sum,na.rm=TRUE)  
      }
      if(var=="rx1day"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),max,na.rm=TRUE)  
      }
    }
    tmp2 = apply(tmp,c(1,2),mean,na.rm=TRUE)
    tmp2 = ifelse(is.na(vardat[,,1])==FALSE,tmp2,NA)
    histclimo[,,i,s] = tmp2
    }
  }
  
  message("Finished Calcs for ",i," / ",length(histfiles))
}

if(var=="tasmax"){
  meanvals = apply(histclimo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>200)
  if(length(fixidx)>0){
    histclimo[,,fixidx,] = histclimo[,,fixidx,]-273.15  
  }
}

if(var=="tmax100"){
  meanvals = apply(histclimo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>300)
  if(length(fixidx)>0){
    histclimo[,,fixidx] = NA  
  }
}

if(var=="pr" | var=="rx1day"){
  meanvals = apply(histclimo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals<1)
  if(length(fixidx)>0){
    histclimo[,,fixidx,] = histclimo[,,fixidx,]*86400  
  }
  #histclimo = histclimo*10
}


################
# Determine rcp45 climatology

projdates = seq(as.Date("2036-01-15"),as.Date("2065-12-15"),by="month")
years = unique(as.numeric(substr(projdates,1,4)))

for(i in 1:length(rcp45files)){
  ptm = proc.time()
  nctest = nc_open(rcp45files[i])
  
  if(i==1){
    lon=ncvar_get(nctest,"lon")
    lat=ncvar_get(nctest,"lat")
    if(period=="ann"){
      rcp45climo = array(NA,dim=c(length(lon),length(lat),length(rcp45files)))
    }
    if(period=="sea"){
      rcp45climo = array(NA,dim=c(length(lon),length(lat),length(rcp45files),4))
    }
  }
  
  vardat = ncvar_get(nctest,var)
  nc_close(nctest)
  
  if(period=="ann"){
    tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      yearidx = which(as.numeric(substr(projdates,1,4))==years[y])
      if(var=="tasmax" | var=="tasmin"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),mean,na.rm=TRUE)  
      }
      if(var=="pr" | var=="tmax100"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),sum,na.rm=TRUE)  
      }
      if(var=="rx1day"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),max,na.rm=TRUE)  
      }
    }
    tmp2 = apply(tmp,c(1,2),mean,na.rm=TRUE)
    tmp2 = ifelse(is.na(vardat[,,1])==FALSE,tmp2,NA)
    rcp45climo[,,i] = tmp2
  }
  
  if(period=="sea"){
    for(s in 1:4){
      tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
      for(y in 1:length(years)){
        if(s==1){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & (as.numeric(substr(projdates,6,7))<=2 | as.numeric(substr(projdates,6,7))==12))
        }
        if(s==2){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & as.numeric(substr(projdates,6,7))>=3 & as.numeric(substr(projdates,6,7))<=5)
        }
        if(s==3){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & as.numeric(substr(projdates,6,7))>=6 & as.numeric(substr(projdates,6,7))<=8)
        }
        if(s==4){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & as.numeric(substr(projdates,6,7))>-9 & as.numeric(substr(projdates,6,7))<=11)
        }
        if(var=="tasmax" | var=="tasmin"){
          tmp[,,y] = apply(vardat[,,yearidx],c(1,2),mean,na.rm=TRUE)  
        }
        if(var=="pr" | var=="tmax100"){
          tmp[,,y] = apply(vardat[,,yearidx],c(1,2),sum,na.rm=TRUE)  
        }
        if(var=="rx1day"){
          tmp[,,y] = apply(vardat[,,yearidx],c(1,2),max,na.rm=TRUE)  
        }
      }
      tmp2 = apply(tmp,c(1,2),mean,na.rm=TRUE)
      tmp2 = ifelse(is.na(vardat[,,1])==FALSE,tmp2,NA)
      rcp45climo[,,i,s] = tmp2
    }
  }
  ptmend = proc.time()-ptm
  message("Finished Calcs for ",i," / ",length(rcp45files)," - done in ",ptmend[3]," secs")
}


if(var=="tasmax"){
  meanvals = apply(rcp45climo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>200)
  if(length(fixidx)>0){
  rcp45climo[,,fixidx,] = rcp45climo[,,fixidx,]-273.15
  }
}

if(var=="tmax100"){
  meanvals = apply(rcp45climo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>300)
  if(length(fixidx)>0){
    rcp45climo[,,fixidx] = NA  
  }
}

if(var=="pr" | var=="rx1day"){
  meanvals = apply(rcp45climo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>6000)
  if(length(fixidx)>0){
    rcp45climo[,,fixidx,] = rcp45climo[,,fixidx,]/86400
  }
  #rcp45climo = rcp45climo*10
}


################
# Determine rcp85 climatology

for(i in 1:length(rcp85files)){
  
  nctest = nc_open(rcp85files[i])
  
  if(i==1){
    lon=ncvar_get(nctest,"lon")
    lat=ncvar_get(nctest,"lat")
    if(period=="ann"){
      rcp85climo = array(NA,dim=c(length(lon),length(lat),length(rcp85files)))
    }
    if(period=="sea"){
      rcp85climo = array(NA,dim=c(length(lon),length(lat),length(rcp85files),4))
    }
  }
  
  vardat = ncvar_get(nctest,var)
  nc_close(nctest)
  
  if(period=="ann"){
    tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)){
      yearidx = which(as.numeric(substr(projdates,1,4))==years[y])
      if(var=="tasmax" | var=="tasmin"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),mean,na.rm=TRUE)  
      }
      if(var=="pr" | var=="tmax100"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),sum,na.rm=TRUE)  
      }
      if(var=="rx1day"){
        tmp[,,y] = apply(vardat[,,yearidx],c(1,2),max,na.rm=TRUE)  
      }
    }
    tmp2 = apply(tmp,c(1,2),mean,na.rm=TRUE)
    tmp2 = ifelse(is.na(vardat[,,1])==FALSE,tmp2,NA)
    rcp85climo[,,i] = tmp2
  }
  
  if(period=="sea"){
    for(s in 1:4){
      tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
      for(y in 1:length(years)){
        if(s==1){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & (as.numeric(substr(projdates,6,7))<=2 | as.numeric(substr(projdates,6,7))==12))
        }
        if(s==2){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & as.numeric(substr(projdates,6,7))>=3 & as.numeric(substr(projdates,6,7))<=5)
        }
        if(s==3){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & as.numeric(substr(projdates,6,7))>=6 & as.numeric(substr(projdates,6,7))<=8)
        }
        if(s==4){
          yearidx = which(as.numeric(substr(projdates,1,4))==years[y] & as.numeric(substr(projdates,6,7))>-9 & as.numeric(substr(projdates,6,7))<=11)
        }
        if(var=="tasmax" | var=="tasmin"){
          tmp[,,y] = apply(vardat[,,yearidx],c(1,2),mean,na.rm=TRUE)  
        }
        if(var=="pr" | var=="tmax100"){
          tmp[,,y] = apply(vardat[,,yearidx],c(1,2),sum,na.rm=TRUE)  
        }
        if(var=="rx1day"){
          tmp[,,y] = apply(vardat[,,yearidx],c(1,2),max,na.rm=TRUE)  
        }
      }
      tmp2 = apply(tmp,c(1,2),mean,na.rm=TRUE)
      tmp2 = ifelse(is.na(vardat[,,1])==FALSE,tmp2,NA)
      rcp85climo[,,i,s] = tmp2
    }
  }
  
  message("Finished Calcs for ",i," / ",length(rcp85files))
}


if(var=="tasmax"){
  meanvals = apply(rcp85climo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>200)
  if(length(fixidx)>0){
    rcp85climo[,,fixidx,] = rcp85climo[,,fixidx,]-273.15
  }
}


if(var=="tmax100"){
  meanvals = apply(rcp85climo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>300)
  if(length(fixidx)>0){
    rcp85climo[,,fixidx] = NA  
  }
}


if(var=="pr" | var=="rx1day"){
  meanvals = apply(rcp85climo,3,mean,na.rm=TRUE)
  fixidx = which(meanvals>6000)
  if(length(fixidx)>0){
    rcp85climo[,,fixidx,] = rcp85climo[,,fixidx,]/86400
  }
  #rcp85climo = rcp85climo*10
}


#######
# Projected Changes

rcp45change = rcp45climo-histclimo
rcp85change = rcp85climo-histclimo

if(period=="ann"){
  rcp45ensmean = apply(rcp45change,c(1,2),mean,na.rm=TRUE)
  rcp85ensmean = apply(rcp85change,c(1,2),mean,na.rm=TRUE)
} 

if(period=="sea"){
  rcp45ensmean = apply(rcp45change,c(1,2,4),mean,na.rm=TRUE)
  rcp85ensmean = apply(rcp85change,c(1,2,4),mean,na.rm=TRUE)
}
#########
# plot ensemble mean maps - annual

if(period=="ann"){

plotfilename = paste("/home/woot0002/EnsMean_",var,"_LOCA_absolute_",period,"_change.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=20,height=7)

if(var=="tmax100"){
  rcp45idxlow = which(rcp45ensmean<0)
  rcp85idxlow = which(rcp85ensmean<0)
  if(length(rcp45idxlow)>0){
    rcp45ensmean = ifelse(rcp45ensmean<0,0.01,rcp45ensmean)
  }
  if(length(rcp85idxlow)>0){
    rcp85ensmean = ifelse(rcp85ensmean<0,0.01,rcp85ensmean)
  }
}

ensmeandat = c(rcp45ensmean,rcp85ensmean)
ensmeandat = ensmeandat[-which(is.na(ensmeandat)==TRUE)]

diffcolorbar = colorramp(ensmeandat,colorchoice=colorchoicediff,Blimit=BINLIMIT)

par(mfrow=c(1,2))
testsfc1 = list(x=lon,y=lat,z=rcp45ensmean)
surface(testsfc1,type="I",main="RCP 4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
map("state",add=TRUE)
map("world",add=TRUE)
text(-93,28.45,labels=paste("MAX = ",round(max(rcp45ensmean,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp45ensmean,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-93,26.45,labels=paste("MIN = ",round(min(rcp45ensmean,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

testsfc1 = list(x=lon,y=lat,z=rcp85ensmean)
surface(testsfc1,type="I",main="RCP 8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
map("state",add=TRUE)
map("world",add=TRUE)
text(-93,28.45,labels=paste("MAX = ",round(max(rcp85ensmean,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp85ensmean,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
text(-93,26.45,labels=paste("MIN = ",round(min(rcp85ensmean,na.rm=TRUE),1),sep=""),cex=1,pos = 4)

if(var=="tasmax"){
  
  diffcolorbar = colorramp(c(rcp45ensmean,rcp85ensmean)*9/5,colorchoice=colorchoicediff,Blimit=BINLIMIT)
  
  par(mfrow=c(1,2))
  testsfc1 = list(x=lon,y=lat,z=rcp45ensmean*9/5)
  surface(testsfc1,type="I",main="RCP 4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
  map("state",add=TRUE)
  map("world",add=TRUE)
  text(-93,28.45,labels=paste("MAX = ",round(max(rcp45ensmean*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp45ensmean*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,26.45,labels=paste("MIN = ",round(min(rcp45ensmean*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
  testsfc1 = list(x=lon,y=lat,z=rcp85ensmean*9/5)
  surface(testsfc1,type="I",main="RCP 8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
  map("state",add=TRUE)
  map("world",add=TRUE)
  text(-93,28.45,labels=paste("MAX = ",round(max(rcp85ensmean*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp85ensmean*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,26.45,labels=paste("MIN = ",round(min(rcp85ensmean*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
}

if(var=="pr" | var=="rx1day"){
  
  diffcolorbar = colorramp(c(rcp45ensmean,rcp85ensmean)/25.4,colorchoice=colorchoicediff,Blimit=BINLIMIT)
  
  par(mfrow=c(1,2))
  testsfc1 = list(x=lon,y=lat,z=rcp45ensmean/25.4)
  surface(testsfc1,type="I",main="RCP 4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
  map("state",add=TRUE)
  map("world",add=TRUE)
  text(-93,28.45,labels=paste("MAX = ",round(max(rcp45ensmean/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp45ensmean/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,26.45,labels=paste("MIN = ",round(min(rcp45ensmean/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
  testsfc1 = list(x=lon,y=lat,z=rcp85ensmean/25.4)
  surface(testsfc1,type="I",main="RCP 8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
  map("state",add=TRUE)
  map("world",add=TRUE)
  text(-93,28.45,labels=paste("MAX = ",round(max(rcp85ensmean/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp85ensmean/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,26.45,labels=paste("MIN = ",round(min(rcp85ensmean/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
}


dev.off()

#######
# Gather data for locations

locdat = read.table("/home/woot0002/LOCAlocations.csv",sep=",",header=TRUE)

###
# create model grid

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

###
# get cells to use
locdat$modlon = NA
locdat$modlat = NA
locdat$R = NA
locdat$C = NA

for(l in 1:nrow(locdat)){
  pointarea = distfunc(locdat$lon[l],locdat$lat[l],modelgrid)
  locdat$R[l] = pointarea$R[1]
  locdat$C[l] = pointarea$C[1]
  locdat$modlon[l] = pointarea$lon[1]
  locdat$modlat[l] = pointarea$lat[1]
}

###
# Create table and boxplots

  vardatatable = NULL
  GCM = rep(rcp45filelist$GCM,2)
  for(l in 1:nrow(locdat)){
    loc = locdat$city[l]
    scen = rep(c("rcp45","rcp85"),each=length(unique(GCM)))
    vardat = c(rcp45change[locdat$R[l],locdat$C[l],],rcp85change[locdat$R[l],locdat$C[l],])
    tmp = data.frame(loc,scen,GCM,vardat)
    vardatatable = rbind(vardatatable,tmp)
  }
  
  ensdatatable_rcp45 = aggregate(vardatatable$vardat[which(vardatatable$scen=="rcp45")],by=list(loc=vardatatable$loc[which(vardatatable$scen=="rcp45")]),mean,na.rm=TRUE)
  names(ensdatatable_rcp45) = c("loc","vardat")
  ensdatatable_rcp45$scen = "rcp45"
  ensdatatable_rcp85 = aggregate(vardatatable$vardat[which(vardatatable$scen=="rcp85")],by=list(loc=vardatatable$loc[which(vardatatable$scen=="rcp85")]),mean,na.rm=TRUE)
  names(ensdatatable_rcp85) = c("loc","vardat")
  ensdatatable_rcp85$scen = "rcp85"
  ensdatatable = rbind(ensdatatable_rcp45,ensdatatable_rcp85)
  
  ggplot(vardatatable, aes(x=scen, y=vardat,fill=scen)) + 
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + 
    scale_fill_manual(values=c("blue","red")) + 
    geom_point(data=ensdatatable, mapping=aes(x=scen,y=vardat),size=2.5,shape=23,fill="black") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed") +
    ggtitle("2036-2065 Projected Change")+xlab("Scenario")+ylab("deg C") +facet_wrap(facets=vars(loc),nrow=3,ncol=3)
  outfile = paste("/home/woot0002/Boxplots",var,"_LOCA_absolute_",period,"_change.pdf",sep="")
  ggsave(outfile,width=10,height=10)
  
  if(var=="tasmax"){
    vardatatable$vardat = vardatatable$vardat*9/5
    ensdatatable$vardat = ensdatatable$vardat*9/5
    ggplot(vardatatable, aes(x=scen, y=vardat,fill=scen)) + 
      stat_boxplot(geom ='errorbar', width = 0.5) +
      geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + 
      scale_fill_manual(values=c("blue","red")) + 
      geom_point(data=ensdatatable, mapping=aes(x=scen,y=vardat),size=2.5,shape=23,fill="black") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed") +
      ggtitle("2036-2065 Projected Change")+xlab("Scenario")+ylab("deg F") +facet_wrap(facets=vars(loc),nrow=3,ncol=3)
    outfile = paste("/home/woot0002/Boxplots",var,"_LOCA_absolute_",period,"_change_F.pdf",sep="")
    ggsave(outfile,width=10,height=10)
  }
  
  if(var=="pr" | var=="rx1day"){
    vardatatable$vardat = vardatatable$vardat/25.4
    ensdatatable$vardat = ensdatatable$vardat/25.4
    ggplot(vardatatable, aes(x=scen, y=vardat,fill=scen)) + 
      stat_boxplot(geom ='errorbar', width = 0.5) +
      geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + 
      scale_fill_manual(values=c("blue","red")) + 
      geom_point(data=ensdatatable, mapping=aes(x=scen,y=vardat),size=2.5,shape=23,fill="black") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed") +
      ggtitle("2036-2065 Projected Change")+xlab("Scenario")+ylab("deg C") +facet_wrap(facets=vars(loc),nrow=3,ncol=3)
    outfile = paste("/home/woot0002/Boxplots",var,"_LOCA_absolute_",period,"_change_in.pdf",sep="")
    ggsave(outfile,width=10,height=10)
  }

}

#########
# plot ensemble mean maps - seasonal

if(period=="sea"){
  
  plotfilename = paste("/home/woot0002/EnsMean_",var,"_LOCA_absolute_",period,"_change.pdf",sep="") 
  
  pdf(plotfilename,onefile=TRUE,width=20,height=14)
  par(mfrow=c(2,2))
  diffcolorbar = colorramp(c(rcp45ensmean,rcp85ensmean),colorchoice=colorchoicediff,Blimit=BINLIMIT)
  for(s in 1:4){
  testsfc1 = list(x=lon,y=lat,z=rcp45ensmean[,,s])
  surface(testsfc1,type="I",main=paste("RCP 4.5 season=",s,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
  map("state",add=TRUE)
  text(-93,28.45,labels=paste("MAX = ",round(max(rcp45ensmean[,,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp45ensmean[,,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,26.45,labels=paste("MIN = ",round(min(rcp45ensmean[,,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  for(s in 1:4){
  testsfc1 = list(x=lon,y=lat,z=rcp85ensmean[,,s])
  surface(testsfc1,type="I",main=paste("RCP 8.5 season=",s,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
  map("state",add=TRUE)
  text(-93,28.45,labels=paste("MAX = ",round(max(rcp85ensmean[,,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp85ensmean[,,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-93,26.45,labels=paste("MIN = ",round(min(rcp85ensmean[,,s],na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  
  if(var=="tasmax"){
    
    diffcolorbar = colorramp(c(rcp45ensmean,rcp85ensmean)*9/5,colorchoice=colorchoicediff,Blimit=BINLIMIT)
    
    for(s in 1:4){
    testsfc1 = list(x=lon,y=lat,z=rcp45ensmean[,,s]*9/5)
    surface(testsfc1,type="I",main=paste("RCP 4.5 season=",s,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
    map("state",add=TRUE)
    text(-93,28.45,labels=paste("MAX = ",round(max(rcp45ensmean[,,s]*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp45ensmean[,,s]*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,26.45,labels=paste("MIN = ",round(min(rcp45ensmean[,,s]*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    }
    
    for(s in 1:4){
    testsfc1 = list(x=lon,y=lat,z=rcp85ensmean[,,s]*9/5)
    surface(testsfc1,type="I",main=paste("RCP 8.5 season=",s,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
    map("state",add=TRUE)
    text(-93,28.45,labels=paste("MAX = ",round(max(rcp85ensmean[,,s]*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp85ensmean[,,s]*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,26.45,labels=paste("MIN = ",round(min(rcp85ensmean[,,s]*9/5,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    }
  }
  
  if(var=="pr" | var=="rx1day"){
    
    diffcolorbar = colorramp(c(rcp45ensmean,rcp85ensmean)/25.4,colorchoice=colorchoicediff,Blimit=BINLIMIT)
    
    for(s in 1:4){
    testsfc1 = list(x=lon,y=lat,z=rcp45ensmean[,,s]/25.4)
    surface(testsfc1,type="I",main=paste("RCP 4.5 season=",s,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
    map("state",add=TRUE)
    text(-93,28.45,labels=paste("MAX = ",round(max(rcp45ensmean[,,s]/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp45ensmean[,,s]/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,26.45,labels=paste("MIN = ",round(min(rcp45ensmean[,,s]/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    } 
    
    for(s in 1:4){
    testsfc1 = list(x=lon,y=lat,z=rcp85ensmean[,,s]/25.4)
    surface(testsfc1,type="I",main=paste("RCP 8.5 season=",s,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude",ylim=c(25,40))
    map("state",add=TRUE)
    text(-93,28.45,labels=paste("MAX = ",round(max(rcp85ensmean[,,s]/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,27.45,labels=paste("MEAN = ",round(mean(rcp85ensmean[,,s]/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-93,26.45,labels=paste("MIN = ",round(min(rcp85ensmean[,,s]/25.4,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }
  }
  
  dev.off()
}





