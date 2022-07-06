
library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("/data2/3to5/I35/scripts/colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maptools)
library(ggplot2)
library(zoo)
library(lars)
library(RPushbullet)

stationfile = read.table("/home/woot0002/stationsNMWSC.csv",sep=",",header=TRUE)

# load DS files

prfile_hist = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_*historical*.nc",intern=TRUE))

tasmaxfile_hist = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_*historical*.nc",intern=TRUE))

tasminfile_hist = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_*historical*.nc",intern=TRUE))

prfile_proj = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_day*rcp*.nc",intern=TRUE))

tasmaxfile_proj = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_day*rcp*.nc",intern=TRUE))

tasminfile_proj = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_day*rcp*.nc",intern=TRUE))

#####
# Create file breakdown tables

filebreakdown = do.call(rbind,strsplit(prfile_proj,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(prfile_hist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

####
# find data point to use in all.

nctest = nc_open(prfile_hist[1])
lon = ncvar_get(nctest,"lon")-360
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

###
# create model grid

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360


####
# start station loop

for(f in 1:nrow(stationfile)){

  ptm = proc.time()
SITEID = as.character(stationfile$Site_ID[f])
loclon = stationfile$Longitude[f] # gauging station location 
loclat = stationfile$Latitude[f]

  ###
  # get cells to use
  
  pointarea = distfunc(loclon,loclat,modelgrid)
  
  ######
  # Start a loop through historical files.
  histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  
  for(i in 1:length(prfile_hist)){
    
    ###
    # Get PRCP, TMAX, and TMIN data
    
    nctest = nc_open(prfile_hist[i])
    PRCP = ncvar_get(nctest,"pr",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))*86400
    nc_close(nctest)
    
    nctest = nc_open(tasmaxfile_hist[i])
    TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
    nc_close(nctest)
    
    nctest = nc_open(tasminfile_hist[i])
    TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
    nc_close(nctest)
    
    print(length(PRCP))
    print(length(TMAX))
    print(length(TMIN))
    
    if(length(PRCP)<length(histdates)){
      idxin = which(substr(histdates,6,10)!="02-29")
        tmp = rep(NA,length(histdates))
        tmp[idxin] = PRCP
        PRCP=tmp
    }
    
    if(length(TMAX)<length(histdates)){
      idxin = which(substr(histdates,6,10)!="02-29")
      tmp = rep(NA,length(histdates))
      tmp[idxin] = TMAX
      TMAX=tmp
    }
    
    if(length(TMIN)<length(histdates)){
      idxin = which(substr(histdates,6,10)!="02-29")
      tmp = rep(NA,length(histdates))
      tmp[idxin] = TMIN
      TMIN=tmp
    }
    
    
    #####
    # unit conversion
    
    PRCP = PRCP/25.4
    TMAX = TMAX-273.15
    TMAX = (TMAX*9/5)+32
    TMIN = TMIN-273.15
    TMIN = (TMIN*9/5)+32
    DATE = histdates
    
    if(i==1){
      prcpdat = data.frame(DATE,PRCP)
      tmaxdat = data.frame(DATE,TMAX)
      tmindat = data.frame(DATE,TMIN)
    } else {
      prcpdat = cbind(prcpdat,PRCP)
      tmaxdat = cbind(prcpdat,TMAX)
      tmindat = cbind(prcpdat,TMIN)
    }
    
    names(prcpdat)[i+1] = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")
    names(tmaxdat)[i+1] = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")
    names(tmindat)[i+1] = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")
    
    message("Finished gathering results for file ",i," / ",length(prfile_hist))
  }
  
  write.table(prcpdat,file=paste("/home/woot0002/NMWSC/",SITEID,"_PRCP_1981-2005.csv",sep=""),sep=",",row.names=FALSE)
  write.table(tmaxdat,file=paste("/home/woot0002/NMWSC/",SITEID,"_TMAX_1981-2005.csv",sep=""),sep=",",row.names=FALSE)
  write.table(tmindat,file=paste("/home/woot0002/NMWSC/",SITEID,"_TMIN_1981-2005.csv",sep=""),sep=",",row.names=FALSE)
  
  ######
  # Start a loop through historical files.
  projdates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
  
  for(i in 1:length(prfile_proj)){
    
    ###
    # Get PRCP, TMAX, and TMIN data
    
    nctest = nc_open(prfile_proj[i])
    PRCP = ncvar_get(nctest,"pr",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))*86400
    nc_close(nctest)
    
    nctest = nc_open(tasmaxfile_proj[i])
    TMAX = ncvar_get(nctest,"tasmax",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
    nc_close(nctest)
    
    nctest = nc_open(tasminfile_proj[i])
    TMIN = ncvar_get(nctest,"tasmin",start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
    nc_close(nctest)
    
    print(length(PRCP))
    print(length(TMAX))
    print(length(TMIN))
    
    if(length(PRCP)<length(projdates)){
      idxin = which(substr(projdates,6,10)!="02-29")
      tmp = rep(NA,length(projdates))
      tmp[idxin] = PRCP
      PRCP=tmp
    }
    
    if(length(TMAX)<length(projdates)){
      idxin = which(substr(projdates,6,10)!="02-29")
      tmp = rep(NA,length(projdates))
      tmp[idxin] = TMAX
      TMAX=tmp
    }
    
    if(length(TMIN)<length(projdates)){
      idxin = which(substr(projdates,6,10)!="02-29")
      tmp = rep(NA,length(projdates))
      tmp[idxin] = TMIN
      TMIN=tmp
    }
    
    #####
    # unit conversion
    
    PRCP = PRCP/25.4
    TMAX = TMAX-273.15
    TMAX = (TMAX*9/5)+32
    TMIN = TMIN-273.15
    TMIN = (TMIN*9/5)+32
    DATE = projdates
    
    if(i==1){
      prcpdat = data.frame(DATE,PRCP)
      tmaxdat = data.frame(DATE,TMAX)
      tmindat = data.frame(DATE,TMIN)
    } else {
      prcpdat = cbind(prcpdat,PRCP)
      tmaxdat = cbind(prcpdat,TMAX)
      tmindat = cbind(prcpdat,TMIN)
    }
    
    names(prcpdat)[i+1] = paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_")
    names(tmaxdat)[i+1] = paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_")
    names(tmindat)[i+1] = paste(projfilebreakdown$scen[i],projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_")
    
    message("Finished gathering results for file ",i," / ",length(prfile_proj))
  }
  
  write.table(prcpdat,file=paste("/home/woot0002/NMWSC/",SITEID,"_PRCP_2006-2099.csv",sep=""),sep=",",row.names=FALSE)
  write.table(tmaxdat,file=paste("/home/woot0002/NMWSC/",SITEID,"_TMAX_2006-2099.csv",sep=""),sep=",",row.names=FALSE)
  write.table(tmindat,file=paste("/home/woot0002/NMWSC/",SITEID,"_TMIN_2006-2099.csv",sep=""),sep=",",row.names=FALSE)
  
  ptmend = proc.time()-ptm
  
  pbPost(type="note",title="NMWSC files written out",body=paste("Data files written for ",SITEID," and it took ",ptmend[3]," secs.",sep=""),email="amwootte@ou.edu")
  
}