####
# Following the creation of individual monthly files, this script and another (GCM_projclimov2.R) combine GCM and obs files together and calculate a full ensemble mean from the GCM files that is also included in the output
# This facilitates the calculations for the ensemble subset selection. 

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

#####
# RCMES combine files - doing a regrid and combine

var = "pr"
OBSfile = paste("/home/woot0002/RCMES/data/OBS/",var,"_DAYMETv4_mon_CONUS_1980-2014.nc",sep="") # file created by Daymetclimov2.R

filelist = system(paste("ls /data4/data/CMIP6/monthlyclimo/",var,"*.nc",sep=""),intern=TRUE) # files created by GCM_monthly_calc_generic.R
GCMfiles = filelist[grep("1980-2014",filelist)]

climotype = "ann"

mask = TRUE
maskregion = "SGP-NCA"
maskfile = "/home/woot0002/RCMES/masks/SGP-NCA_DAYMETv4-1km_mask.nc" # This script can use a mask file provided to it
# mask file is based on DayMETv4 files and is the Daymetv4 grid in geographic coordinates. 
# Locations in the SGP-NCA region (Oklahoma, Kansas, and Texas) have values of 1, the other locations have values of 0.
# can also use the information below to help with combination.

if(mask==TRUE){
  if(maskregion=="NE-NCA"){
    xlimin = c(-86,-64)
    ylimin = c(36,48)
 
  }
  if(maskregion=="SW-NCA"){
    xlimin = c(-125,-102)
    ylimin = c(31,42)

  }
  if(maskregion=="SGP-NCA"){
    xlimin = c(-107,-93)
    ylimin = c(25,40)
  }
  if(maskregion=="C-PrEP"){
    xlimin = c(-110,-88)
    ylimin = c(25,43)
  }
  if(maskregion=="SE-NCA"){
    xlimin = c(-95,-75)
    ylimin = c(24,40)

  }
  if(maskregion=="MW-NCA"){
    xlimin = c(-99,-80)
    ylimin = c(36,50)

  }
  if(maskregion=="NGP-NCA"){
    xlimin = c(-118,-95)
    ylimin = c(40,49)

  }
  if(maskregion=="NW-NCA"){
    xlimin = c(-125,-111)
    ylimin = c(42,49)

  }
} 

####
# read OBS data

nctest = nc_open(OBSfile)
lat = ncvar_get(nctest,"lat")
lon = ncvar_get(nctest,"lon")
#times = ncvar_get(nctest,"time")

varnames = c()
#outputdata = list()
for(i in 1:length(nctest$var)){
  varnames[i]=nctest$var[[i]]$name
#  outputdata[[i]]=ncvar_get(nctest,varnames[i])
}

varunits = nctest$var[[1]]$units

dimunits=c()
for(i in 1:length(nctest$dim)){
  dimunits[i]=nctest$dim[[i]]$units
}

nc_close(nctest)

####
# Calculate climatology - monthly climatology calculation

histdates = seq(as.Date("1980-01-15"),as.Date("2014-12-15"),by="month")
years = 1980:2014
if(climotype=="monthly"){
  
for(m in 1:12){
for(y in 1:length(years)){
  nctest = nc_open(OBSfile)
  idx = which(as.numeric(substr(histdates,1,4))==years[y] & as.numeric(substr(histdates,6,7))==m)
  if(y==1){
    lonidx = which(lon>= min(lon) & lon<= max(lon))
    latidx = which(lat>= min(lat) & lat<= max(lat))
  }
  
  outputdata=ncvar_get(nctest,varnames[1],start=c(lonidx[1],latidx[1],idx[1]),count=c(length(lonidx),length(latidx),length(idx)))
  if(y==1 & m==1){
    obsclimo = array(NA,dim=c(length(lonidx),length(latidx),12))
  }
  
  if(y==1){
    obsyearly = array(NA,dim=c(length(lonidx),length(latidx),length(years)))
    obsyearly[,,y] = outputdata
  } else {
    obsyearly[,,y] = outputdata
  }
  
  nc_close(nctest)
  message("Done calcs for year ",years[y])
}
  obsclimo[,,m] = apply(obsyearly,c(1,2),mean,na.rm=TRUE)
  message("Done monthly climo calcs for month ",m)
}
}


####
# Calculate climatology - annual climatology calculation

histdates = seq(as.Date("1980-01-15"),as.Date("2014-12-15"),by="month")
years = 1980:2014
if(climotype=="annual"){
  
    for(y in 1:length(years)){
      nctest = nc_open(OBSfile)
      idx = which(as.numeric(substr(histdates,1,4))==years[y])
     
      if(y==1){
        lonidx = which(lon>= min(lon) & lon<= max(lon))
        latidx = which(lat>= min(lat) & lat<= max(lat))
      }
      
      outputdata=ncvar_get(nctest,varnames[1],start=c(lonidx[1],latidx[1],idx[1]),count=c(length(lonidx),length(latidx),length(idx)))
      
      if(y==1){
        obsyearly = array(NA,dim=c(length(lonidx),length(latidx),length(years)))
      }
      if(var=="pr" | var=="r1mm"){
        obsyearly[,,y] = apply(outputdata,c(1,2),sum,na.rm=TRUE)
      }
      if(var=="rx1day" | var=="rx5day"){
        obsyearly[,,y] = apply(outputdata,c(1,2),max,na.rm=TRUE)
      }
      if(var=="tasmax" | var=="tasmin"){
        obsyearly[,,y] = apply(outputdata,c(1,2),mean,na.rm=TRUE)
      }
      
      nc_close(nctest)
      message("Done calcs for year ",years[y])
    }
    obsclimo = apply(obsyearly,c(1,2),mean,na.rm=TRUE)
}


####
# Calculate climatology - single month climatology calculation

histdates = seq(as.Date("1980-01-15"),as.Date("2014-12-15"),by="month")
years = 1980:2014
if(climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"| climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP" | climotype=="OCT" | climotype=="NOV" | climotype=="DEC"){
  
  if(climotype=="JAN") m=1
  if(climotype=="FEB") m=2
  if(climotype=="MAR") m=3
  if(climotype=="APR") m=4
  if(climotype=="MAY") m=5
  if(climotype=="JUN") m=6
  if(climotype=="JUL") m=7
  if(climotype=="AUG") m=8
  if(climotype=="SEP") m=9
  if(climotype=="OCT") m=10
  if(climotype=="NOV") m=11
  if(climotype=="DEC") m=12
  
  for(y in 1:length(years)){
    nctest = nc_open(OBSfile)
    idx = which(as.numeric(substr(histdates,1,4))==years[y] & as.numeric(substr(histdates,6,7))==m)
    
    if(y==1){
      lonidx = which(lon>= min(lon) & lon<= max(lon))
      latidx = which(lat>= min(lat) & lat<= max(lat))
      obsyearly = array(NA,dim=c(length(lonidx),length(latidx),length(years)))
    }
    
    obsyearly[,,y]=ncvar_get(nctest,varnames[1],start=c(lonidx[1],latidx[1],idx[1]),count=c(length(lonidx),length(latidx),length(idx)))
    
    nc_close(nctest)
    message("Done calcs for year ",years[y])
  }
  obsclimo = apply(obsyearly,c(1,2),mean,na.rm=TRUE)
}

lonobs = lon[lonidx]
latobs = lat[latidx]

#####
# mask read in
if(mask==TRUE){
nctest = nc_open(maskfile)
maskinfo =ncvar_get(nctest,"mask",start=c(lonidx[1],latidx[1]),count=c(length(lonidx),length(latidx)))
nc_close(nctest)

if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"| climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP" | climotype=="OCT" | climotype=="NOV" | climotype=="DEC"){
  obsclimo = ifelse(maskinfo==1,obsclimo,NA)  
}
if(climotype=="monthly"){
  for(m in 1:12) obsclimo[,,m] = ifelse(maskinfo==1,obsclimo[,,m],NA)  
}
}

#####
# grab GCM files, calculate climo, and regrid

climodat = list()

for(f in 1:length(GCMfiles)){
  #tmp = ncvar_get(nctest,"ENS",start=c(1,1,1),count=c(-1,-1,1))
  
  nctest = nc_open(GCMfiles[f])
  lat = ncvar_get(nctest,"lat")
  lon=ncvar_get(nctest,"lon")-360
  times = ncvar_get(nctest,"time")
  
  varnames = nctest$var[[length(nctest$var)]]$name
  outputdata=ncvar_get(nctest,varnames)
  
  #varnames = c()
  #outputdata = list()
  #for(i in 1:length(nctest$var)){
  #  varnames[i]=nctest$var[[i]]$name
  #  outputdata[[i]]=ncvar_get(nctest,varnames[i])
  #}
  
  #checkdat = apply(outputdata[[1]],3,mean,na.rm=TRUE)
  #idxfix = which(checkdat>10000)
  #if(length(idxfix)>0){
  #  outputdata[[1]][,,idxfix] = outputdata[[1]][,,idxfix]/86400
  #}
  
  varunits = nctest$var[[1]]$units
  
  dimunits=c()
  for(i in 1:length(nctest$dim)){
    dimunits[i]=nctest$dim[[i]]$units
  }
  
  nc_close(nctest)
  
  if(dim(outputdata)[3]>300){
    Y= length(years)
  } else {
    Y=length(years)-1
  }
  
  ####
  # single month climotype
  
  if(climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"| climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP" | climotype=="OCT" | climotype=="NOV" | climotype=="DEC"){
    outtmp = array(NA,dim=c(length(lon),length(lat),Y))
    if(climotype=="JAN") m=1
    if(climotype=="FEB") m=2
    if(climotype=="MAR") m=3
    if(climotype=="APR") m=4
    if(climotype=="MAY") m=5
    if(climotype=="JUN") m=6
    if(climotype=="JUL") m=7
    if(climotype=="AUG") m=8
    if(climotype=="SEP") m=9
    if(climotype=="OCT") m=10
    if(climotype=="NOV") m=11
    if(climotype=="DEC") m=12
    
    for(y in 1:Y){
      yearidx = which(as.numeric(substr(histdates,1,4))==years[y] & as.numeric(substr(histdates,6,7))==m)
      outtmp[,,y] = ifelse(is.na(outputdata[,,yearidx])==FALSE,outputdata[,,yearidx],NA)
      outtmp[,,y] = ifelse(is.na(outputdata[,,1])==FALSE,outtmp[,,y],NA)
      message("Finished with year ",y," / ",length(years))
    }
    
    climotmp = apply(outtmp,c(1,2),mean,na.rm=TRUE)
    
    ### regrid
    
    int1 = interp.surface.grid(list(x=lon,y=lat,z=climotmp),grid.list=list(x=lonobs,y=latobs))
    newmatrix = matrix(int1[[3]],nrow=length(lonobs),ncol=length(latobs))
    if(var=="pr"){
      newmatrix = ifelse(newmatrix<0,0,newmatrix)
    }
    
    climodat[[f]] = ifelse(is.na(obsclimo)==FALSE,newmatrix,NA)
  }
  
  ####
  # annual climotype
  
  if(climotype=="annual"){
  outtmp = array(NA,dim=c(length(lon),length(lat),Y))
  for(y in 1:Y){
    yearidx = which(as.numeric(substr(histdates,1,4))==years[y])
    
    tmpin = ifelse(is.na(outputdata[,,yearidx])==FALSE,outputdata[,,yearidx],NA)
    
    if(var=="rx5day"){
      outtmp[,,y] = apply(tmpin,c(1,2),max,na.rm=FALSE)
    }
    if(var=="rx1day"){
      outtmp[,,y] = apply(tmpin,c(1,2),max,na.rm=FALSE)
    }
    if(var=="pr" | var=="r1mm"){
      outtmp[,,y] = apply(tmpin,c(1,2),sum,na.rm=FALSE)
    }
    if(var=="tasmax" | var=="tasmin"){
      outtmp[,,y] = apply(tmpin,c(1,2),mean,na.rm=FALSE)
    }
    outtmp[,,y] = ifelse(is.na(outputdata[,,1])==FALSE,outtmp[,,y],NA)
    message("Finished with year ",y," / ",length(years))
  }
  
  climotmp = apply(outtmp,c(1,2),mean,na.rm=TRUE)
  
  ### regrid
  
  int1 = interp.surface.grid(list(x=lon,y=lat,z=climotmp),grid.list=list(x=lonobs,y=latobs))
  newmatrix = matrix(int1[[3]],nrow=length(lonobs),ncol=length(latobs))
  if(var=="pr"){
    newmatrix = ifelse(newmatrix<0,0,newmatrix)
  }
  
  climodat[[f]] = ifelse(is.na(obsclimo)==FALSE,newmatrix,NA)
  }
  
  #####
  # monthly climatology type
  
  if(climotype=="monthly"){
    climotmp = array(NA,dim=c(length(lonobs),length(latobs),12))
    for(m in 1:12){
      monidx = which(as.numeric(substr(histdates,6,7))==m)
      
      tmpin = outputdata[,,monidx]
      tmpin2 = apply(tmpin,c(1,2),mean,na.rm=FALSE)
    
      tmpin2 = ifelse(is.na(outputdata[,,1])==FALSE,tmpin2,NA)
      message("tmpin2 value = ",mean(tmpin2,na.rm=TRUE))
      
      ### regrid
      int1 = interp.surface.grid(list(x=lon,y=lat,z=tmpin2),grid.list=list(x=lonobs,y=latobs))
      newmatrix = matrix(int1[[3]],nrow=length(lonobs),ncol=length(latobs))
      if(var=="pr"){
        newmatrix = ifelse(newmatrix<0,0,newmatrix)
      }
      message("climotmp mean value = ",mean(newmatrix,na.rm=TRUE))
      climotmp[,,m] = ifelse(is.na(obsclimo[,,1])==FALSE,newmatrix,NA)
      message("Finished with m ",m," / 12")
    }
  
    climodat[[f]] = climotmp
  }
  
  message("Finished with GCM file ",f," / ",length(GCMfiles))
}

####
# Ensemble mean calculation
if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"| climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP" | climotype=="OCT" | climotype=="NOV" | climotype=="DEC"){
  climodat2 = array(NA,dim=c(length(lonobs),length(latobs),length(GCMfiles)))
  for(f in 1:length(GCMfiles)){
    climodat2[,,f] = climodat[[f]]
  }
  ENSdat = apply(climodat2,c(1,2),mean,na.rm=TRUE)

} else {
  climodat2 = array(NA,dim=c(length(lonobs),length(latobs),12,length(GCMfiles)))
  for(f in 1:length(GCMfiles)){
    climodat2[,,,f] = climodat[[f]]
  }
  ENSdat = apply(climodat2,c(1,2,3),mean,na.rm=TRUE)
}

if(var=="tasmax" | var=="tasmin"){
  obsclimo = obsclimo+273.15
}


#####
# write netcdf

if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"| climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP" | climotype=="OCT" | climotype=="NOV" | climotype=="DEC"){

filesplit = do.call("rbind",strsplit(GCMfiles,"_",fixed=TRUE))
GCMs = paste(filesplit[,2],filesplit[,3],sep="_")

dimX <- ncdim_def( "lon", dimunits[1],lonobs)
dimY <- ncdim_def( "lat", dimunits[2],latobs)

mv <- 1E20 # missing value to use

var1d = ncvar_def("OBS",varunits, list(dimX,dimY), mv ,compression=9)
for(i in 1:length(GCMs)){
  assign(paste("var",(i+1),"d",sep=""),ncvar_def(GCMs[i],varunits, list(dimX,dimY), mv ,compression=9)) 
}
assign(paste("var",(i+2),"d",sep=""),ncvar_def("ENS",varunits, list(dimX,dimY), mv ,compression=9)) 

fileout = paste("/home/woot0002/RCMES/",var,"_CPREP2_climo_",maskregion,"_",climotype,".nc",sep="")

nc <- nc_create(fileout , eval(parse(text=paste("list(",paste("var",1:(i+2),"d",sep="",collapse=","),")",sep=""))),force_v4 = TRUE)


#nc <- nc_create(fileout , list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d,var11d,
#                               var12d,var13d,var14d,var15d,var16d,var17d,var18d,var19d,var20d,var21d,
#                               var22d,var23d,var24d,var25d,var26d,var27d,var28d,var29d,var30d,var31d),force_v4 = TRUE)

eval(parse(text=paste("ncvar_put(nc, var",1,"d, obsclimo)",sep="")))

for(j in 2:(length(GCMs)+1)){
  eval(parse(text=paste("ncvar_put(nc, var",j,"d, climodat[[",j-1,"]])",sep="")))
}

eval(parse(text=paste("ncvar_put(nc, var",j+1,"d, ENSdat)",sep="")))
nc_close(nc)

system(paste("chmod 777 ",fileout,sep=""))

} 

######
# monthly climo type write netcdf


if(climotype=="monthly"){
  
  filesplit = do.call("rbind",strsplit(GCMfiles,"_",fixed=TRUE))
  GCMs = paste(filesplit[,2],filesplit[,3],sep="_")
  
  dimX <- ncdim_def( "lon", dimunits[1],lonobs)
  dimY <- ncdim_def( "lat", dimunits[2],latobs)
  dimT <- ncdim_def("time","month",1:12)
  
  mv <- 1E20 # missing value to use
  
  var1d = ncvar_def("OBS",varunits, list(dimX,dimY,dimT), mv ,compression=9)
  for(i in 1:length(GCMs)){
    assign(paste("var",(i+1),"d",sep=""),ncvar_def(GCMs[i],varunits, list(dimX,dimY,dimT), mv ,compression=9)) 
  }
  assign(paste("var",(i+2),"d",sep=""),ncvar_def("ENS",varunits, list(dimX,dimY,dimT), mv ,compression=9)) 
  
  fileout = paste("/home/woot0002/RCMES/",var,"_CPREP2_climo_",maskregion,"_monthly.nc",sep="")
  
  nc <- nc_create(fileout , eval(parse(text=paste("list(",paste("var",1:(i+2),"d",sep="",collapse=","),")",sep=""))),force_v4 = TRUE)
  
  
  #nc <- nc_create(fileout , list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d,var11d,
  #                               var12d,var13d,var14d,var15d,var16d,var17d,var18d,var19d,var20d,var21d,
  #                               var22d,var23d,var24d,var25d,var26d,var27d,var28d,var29d,var30d,var31d),force_v4 = TRUE)
  
  eval(parse(text=paste("ncvar_put(nc, var",1,"d, obsclimo)",sep="")))
  
  for(j in 2:(length(GCMs)+1)){
    eval(parse(text=paste("ncvar_put(nc, var",j,"d, climodat[[",j-1,"]])",sep="")))
  }
  
  eval(parse(text=paste("ncvar_put(nc, var",j+1,"d, ENSdat)",sep="")))
  nc_close(nc)
  
  system(paste("chmod 777 ",fileout,sep=""))
  
} 

