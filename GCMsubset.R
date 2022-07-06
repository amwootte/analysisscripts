source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(RPushbullet)

var = "pr"

LOCAfiles = system(paste("ls /home/woot0002/trimmedfiles/LOCA/",var,"_*.nc",sep=""),intern=TRUE)

domain = "day"
experiment = "historical"
latbnds = c(24,42)
lonbnds = c(-112,-88)

infosplit = do.call(rbind,strsplit(LOCAfiles,"_",fixed=TRUE))

#https://climatedata.oscer.ou.edu/climatedata_updated-oscer_path.csv

GCMtable = read.table("/home/woot0002/climatedata_full_CMIP5_external.csv",sep=",",header=TRUE,colClasses = "character")
GCMtable = subset(GCMtable,variable==var)

idx1 = which(GCMtable$domain == domain)

GCMtable = GCMtable[idx1,]
idx2= which(GCMtable$experiment == experiment)

GCMtable = GCMtable[idx2,]

idxs = c()

for(r in 1:nrow(infosplit)){
  
  idxout = which(GCMtable$model == infosplit[r,2] & GCMtable$ensemble == infosplit[r,3])
  idxs = c(idxs,idxout)
  
}

GCMtable2 = GCMtable[idxs,]

startyear = as.numeric(substr(GCMtable2$time,1,4))
endyear = as.numeric(substr(GCMtable2$time,10,13))

years = 1976:2005
checkidx = c()

for(i in 1:length(startyear)){
  fileperiod = startyear[i]:endyear[i]
  check = any(years %in% fileperiod)
  if(check==TRUE) checkidx = c(checkidx,i)
}

GCMtable3=GCMtable2[checkidx,]

#GCMtable3 = subset(GCMtable3,latest=="1.0")

combos = unique(paste(GCMtable3$variable,GCMtable3$model,GCMtable3$ensemble,GCMtable3$time,sep="_"))
GCMtable3$combos = paste(GCMtable3$variable,GCMtable3$model,GCMtable3$ensemble,GCMtable3$time,sep="_")

idxsout = c()
for(c in 1:length(combos)){
  idx = which(GCMtable3$combos==combos[c])
  if(length(idx)>1){
    idxsout = c(idxsout,idx[1])
  } else {
    idxsout = c(idxsout,idx)
  }
}

GCMtable4 = GCMtable3[idxsout,]
GCMexp = unique(paste(GCMtable4$model,GCMtable4$ensemble,sep="_"))

for(i in 1:length(GCMexp)){
  
  idxuse = which(paste(GCMtable4$model,GCMtable4$ensemble,sep="_")==GCMexp[i])
  filestmp = GCMtable4[idxuse,]
  
  ptm = proc.time()
  
  filesplit = strsplit(filestmp$SCCASC_climatedata_path,split="/",fixed=TRUE)
  filestartyears = fileendyears = c()
  for(f in 1:nrow(filestmp)){
    filesplittmp = strsplit(filesplit[[f]][length(filesplit[[f]])],"_",fixed=TRUE)[[1]]
    fileyears = substr(filesplittmp[length(filesplittmp)],1,nchar(filesplittmp[length(filesplittmp)])-3)
    filestartyears[f] = as.numeric(substr(fileyears,1,4))
    fileendyears[f] = as.numeric(substr(fileyears,10,13))
  }
  
  varin = paste(var,GCMexp,"historical",sep="_")
  
  for(y in 1:length(years)){
    
    fileidx = which(filestartyears<=years[y] & fileendyears>=years[y])
    nctest = nc_open(filestmp$SCCASC_climatedata_path[fileidx])
    
    
  }
  
  
  
  if(i==1){
    dataunits = nctest$var[[1]]$units
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    time=ncvar_get(nctest,"time")
    timeunits = nctest$var[[1]]$dim[[3]]$units
    latunits = nctest$var[[1]]$dim[[2]]$units
    lonunits = nctest$var[[1]]$dim[[1]]$units
    
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    timesin = which(as.numeric(substr(dates,1,4))>=1976)
    timestart = timesin[1]
    timecount = length(timesin)
    
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
  
  testdat = ncvar_get(nctest,varin,start=c(lonstart,latstart,timestart),count=c(loncount,latcount,timecount))
  nc_close(nctest)
  
  #####
  filenameout = paste("/home/woot0002/trimmedfiles/LOCA/",varin,".nc",sep="")
  
  dimX <- ncdim_def( "lon", lonunits, lonout)
  dimY <- ncdim_def( "lat", latunits, latout)
  dimT <- ncdim_def("time",timeunits,time[timesin])
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varin,dataunits,longname=varin, list(dimX,dimY,dimT), mv ,compression=9)
  
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
  
  pbPost(type="note",title="LOCA trimmed file written out",body=paste("Trimmed file written for ",varin," it took ",ptmend[3]," secs.",sep=""),email="amwootte@ou.edu")
  
}
