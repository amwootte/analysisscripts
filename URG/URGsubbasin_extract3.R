#############
#
# NM-WSC project
#
# HRU test daily precipitation data pull

ptm= Sys.time()
.libPaths('/home/amwootte/R/x86_64-redhat-linux-gnu-library/3.3/')
library(iterators)
library(parallel)
library(foreach)
library(doParallel)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)

cl <- makePSOCKcluster(20,outfile="")#makePSOCKcluster(20)
registerDoParallel(cl)

##########
# Set input information

message("Setting inputs")

subbasin = "Blanco"
shapefile = "/home/amwootte/shapefiles/Blanco_trans/Blanco_trans"
shapedimension = "nhm_id"
load("/home/amwootte/weightlist.Rdata")
load("/home/amwootte/gridregion.Rdata")

DS = "PARM"
varname = "tasmin"
period = "future"
timeperiod = c(2006,2099)
scen = "rcp"

message("Inputs set!")

###########
# Get netcdf file list

message("Getting list of applicable netcdfs")
ncfiles = system(paste("ls /condo/climatedata2/3to5/I35/",varname,"/",DS,"/",varname,"_day*",scen,"*.nc",sep=""),intern=TRUE)
GCMs = rep(NA,length(ncfiles))
obs = rep(NA,length(ncfiles))
GCMs[grep("Bd1",ncfiles)] = "CCSM4"
GCMs[grep("Bd2",ncfiles)] = "MIROC5"
GCMs[grep("Bd3",ncfiles)] = "MPI-ESM-LR"
obs[grep("D0",ncfiles)] = "Daymet"
obs[grep("L0",ncfiles)] = "Livneh"
obs[grep("P0",ncfiles)] = "PRISM"

filesplit = do.call("rbind",strsplit(ncfiles,"_",fixed=TRUE))
scens = filesplit[,4]
rm(filesplit)

message("Gathered file list")

#nctest = nc_open(ncfiles[1]) # lat and lon are the same for all so get that here too
#lon = ncvar_get(nctest,"lon")
#lat = ncvar_get(nctest,"lat")
#nc_close(nctest)
#rm(nctest)
#message("Gathered common lat/lon grid")
###########
# open and reproject relevant shapefile

message("opening URG shapefile")
test = readShapePoly(shapefile)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

###########
# run calcs for all netcdfs

message("starting parallel calculations for all HUC regions")
#
foreach(n=1:length(ncfiles)) %dopar%  { #n in 1:length(ncfiles) - .packages=c('raster','ncdf4')
  .libPaths('/home/amwootte/R/x86_64-redhat-linux-gnu-library/3.3/')
  library(ncdf4)
  library(raster)
  library(sp)
  library(rgdal)
  
  if(period=="future"){
    dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
  } 
  
  if(period=="historical"){
    dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  }
  
  nctest = nc_open(ncfiles[n])
  times = ncvar_get(nctest,"time")
  
  timelength = length(times)
  #rm(times)
  
  if(length(dates)>timelength){
    datesin = dates[-which(substr(dates,6,10)=="02-29")]
    #rm(dates)
  } else {
    datesin = dates
    #rm(dates)
  }
  
  #timeidx = which(as.numeric(substr(datesin,1,4))>=timeperiod[1] & as.numeric(substr(datesin,1,4))<=timeperiod[2])
  timeidx = 1:length(datesin)
  #rm(datesin)
  
  ##########
  # Calculate area weighted mean for each HRU.
  
  prdaily = matrix(9999,nrow=length(timeidx),ncol=length(eval(parse(text=paste("test@data$",shapedimension,sep=""))))) # initialize matrix with output data.
  #
  #prdaily = matrix(9999,nrow=10,ncol=1021) # initialize matrix with output data.
  
  for(i in 1:length(timeidx)) { # loop over dates to create rasters from netcdf daily data. length(timeidx)
    
    if(varname=="pr"){
      tmpV = ncvar_get(nctest,varname,start=c(1,1,timeidx[i]),count=c(-1,-1,1))*86400
    } else {
      tmpV = ncvar_get(nctest,varname,start=c(1,1,timeidx[i]),count=c(-1,-1,1))
    }
    
    modrasV = raster(t(tmpV)[length(lat):1,])
    #rm(tmpV)
    if(all(lon>0)){
      extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
    } else {
      extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
    }
    areas = eval(parse(text=paste("test@data$",shapedimension,sep="")))
    
    tmp=c()
    
    for(a in 1:length(eval(parse(text=paste("test@data$",shapedimension,sep=""))))){ # total is 1021
      test.sub <- test[as.character(eval(parse(text=paste("test@data$",shapedimension,sep="")))) %in% areas[a], ] # in here are two important details. 1) column name in the data array, 2) item in the data array to subset by
      valueout=extract(modrasV,test.sub,weights=TRUE,fun=mean,na.rm=TRUE)  
      tmp[a] = valueout
    }
    
    prdaily[i,] = tmp
    #rm(tmp)
    #message("Finished Calcs for ",i," / ",length(timeidx))
    
  }
  
  nc_close(nctest)
  outfilename = paste(varname,"_daily_",GCMs[n],"_",DS,"_",obs[n],"_",scens[n],"_",timeperiod[1],"-",timeperiod[2],"_",subbasin,".txt",sep="")
  #outfilename = paste(varname,"_test_",n,".txt",sep="")
  message("Finished Calcs - writing out file")
  prdatfr = data.frame(prdaily)
  names(prdatfr) = paste(shapedimension,areas,sep="_")
  write.table(prdatfr,file=outfilename,sep=" ",row.names=FALSE)
  #rm(prdaily)
  message("Text file written out - ",outfilename)
}

message("Calculations complete!")
ptmend = Sys.time()-ptm
message("Calculation time: ",ptmend," seconds")

stopCluster(cl)
rm(cl)
