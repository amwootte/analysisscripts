source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
#library(RPushbullet)

var = "pr"
varfile = "prcp"

years = 1980:2014
difftype = "absolute"
varunits = "mm"
seasonin = "ann"

filesin = system(paste("ls /home/woot0002/RCMES/data/OBS/DAYMETv4/*",varfile,"*.nc",sep=""),intern=TRUE) # local repository for publicly available DayMETv4. In local repository, this is reprojected from the native coordinate system to geographic coordinates.

filesplit = do.call(rbind,strsplit(filesin,split="/",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,8],split="_",fixed=TRUE))
filesplit2 = data.frame(filesplit2)
filesplit2$year = as.numeric(substr(filesplit2[,ncol(filesplit2)],1,4))

# set map plot arguments
tempres="monthly"

fileidx = which(filesplit2$year>=years[1] & filesplit2$year<=years[length(years)])

filesin = filesin[fileidx]
yearmon = as.character(unique(filesplit2[fileidx,ncol(filesplit2)]))

timeout = c()
yearmonlist = c()
files.out = c()

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  nctest = nc_open(filesin[i])
  
  if(i==1){
    varnames = names(nctest$var)
    dataunits = nctest$var[[which(varnames==varfile)]]$units
    
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    timeunits = nctest$var[[which(varnames==varfile)]]$dim[[3]]$units
    latunits = nctest$var[[which(varnames==varfile)]]$dim[[2]]$units
    lonunits = nctest$var[[which(varnames==varfile)]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    indexdate = as.Date(substr(timeunits,12,21))
    
    loncount = length(lon)
    latcount = length(lat)
    
    #dataout = array(NA,dim=c(loncount,latcount,length(years)*12))
    
  }
  
  times=ncvar_get(nctest,"time")
  dates=indexdate+times
  
  if(any(substr(dates,6,10)=="02-29")==TRUE){
    tmp = c(as.character(dates[1:59]),as.character(dates[61:365]),paste(years[i],"-12-31",sep=""))
    dates2 = tmp
  } else {
    dates2 = dates
  }
  
  if(i==1){
    
    datesall = seq(as.Date(paste(years[1],"-01-15",sep="")),as.Date(paste(years[length(years)],"-12-15",sep="")),by="month")
    timeout = as.vector(datesall-indexdate)
    yearmonlist = substr(datesall,1,7)
    
  }
  
  nc_close(nctest)
  
  if(tempres=="monthly"){
    
    for(j in 1:12){
    
      nctest = nc_open(filesin[i])
      datesidx = which(as.numeric(substr(dates2,6,7))==j)
    
    testdat = ncvar_get(nctest,varfile,start=c(1,1,datesidx[1]),count=c(loncount,latcount,length(datesidx)))
    
      if(var=="pr"){
        tmp1 = apply(testdat,c(1,2),sum,na.rm=TRUE)
        tmp1 = ifelse(tmp1>=0,tmp1,NA)
      } 
    if(var=="r1mm"){
      tmpdat = ifelse(testdat>=1,1,0)
      tmp1 = apply(tmpdat,c(1,2),sum,na.rm=TRUE)
      tmp1 = ifelse(tmp1>=0,testdat[,,1],NA)
    } 
     if(var=="rx1day"){
       tmp1 = apply(testdat,c(1,2),max,na.rm=TRUE)
      }
      if(var=="rx5day"){
        tmp1 = apply(testdat,c(1,2),calcrollsum,size=5)
      }
      if(var=="tasmax" | var=="tasmin"){
        tmp1 = apply(testdat,c(1,2),mean,na.rm=TRUE)
      }
      if(var=="tmax100"){
        td = ifelse(testdat>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      idxout = j+12*(i-1)
      tmp1=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
      message("Finished calcs for month ",j)
      nc_close(nctest)
      rm(testdat)
      gc()
      
      
      yearinfile = years[i]
      file.out <- paste0("/home/woot0002/RCMES/scratch/Daymetslice.",yearinfile,".",j,".",idxout,".Rdata")
      files.out = c(files.out,file.out)
      print(file.out)
      save(tmp1, file=file.out)
      rm(tmp1)
      gc()
      message("Finished with month ",j," / 12")
    }
  }
  message("Finished with year ",i," / ",length(filesin))
}

####
# read in Rdata slices and write to netcdf.

for(k in 1:length(files.out)){
  load(files.out[k])
  idxout = k
  
if(k==1){
  
    dimX <- ncdim_def( "lon", "degrees W", lon)
    dimY <- ncdim_def( "lat", "degrees N", lat)
  
  dimT <- ncdim_def("time",timeunits,timeout)
  mv <- 1E20 # missing value to use
  
  fileout = paste("/home/woot0002/RCMES/data/OBS/",var,"_DAYMETv4_mon_CONUS_",years[1],"-",years[length(years)],".nc",sep="")
  
  var1d <- ncvar_def(var,varunits, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  nc <- nc_create(fileout,  var1d,force_v4 = TRUE )
} else {
  nc = nc_open(fileout,write=TRUE)
}
# Write some data to the file
ncvar_put(nc, var1d, tmp1,start=c(1,1,idxout),count=c(-1,-1,1)) # no start or count: write all values\
# close ncdf
nc_close(nc)
