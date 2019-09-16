####################
#
# CPO Analysis
# Monthly Observation Data Climatology

library(ncdf4)
library(maps)
library(fields)
library(sp)

histperiod = c(1980,2005)
histdates = seq(as.Date("1980-01-01"),as.Date("2005-12-31"),by="day")
histmonths = seq(as.Date("1980-01-15"),as.Date("2005-12-15"),by="month")
latrange = c(26,40)
lonrange = c(251,270)

filepath = "/data4/data/OBS/"
obsdat = "METDATA"

if(obsdat=="LIVNEH"){
  filelist = system(paste("ls ",filepath,obsdat,"/*.nc",sep=""),intern=TRUE)
  
split1 = do.call("rbind",strsplit(filelist,".",fixed=TRUE))
years = as.numeric(substr(split1[,2],1,4))
yearidx = which(years>=histperiod[1] & years<= histperiod[2])

for(y in 1:length(yearidx)){

nctest = nc_open(filelist[yearidx[1]])

if(y==1){
  lon = ncvar_get(nctest,"lon")
  lat = ncvar_get(nctest,"lat")

  if(lon[1]<0) lonrange=lonrange-360
  lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
  latidx = which(lat>= latrange[1] & lat<=latrange[2])
  
  TMAX = TMIN = PR = TMAX95 = TMIN32 = PR50 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
  
}

tmax = ncvar_get(nctest,"Tmax",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))
tmin = ncvar_get(nctest,"Tmin",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))
pr = ncvar_get(nctest,"Prec",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))

nc_close(nctest)

tmax95 = ifelse(tmax>=35,1,0)
tmin32 = ifelse(tmin<=0,1,0)
pr50 = ifelse(pr>=50.4,1,0)

TMAX[,,y] = apply(tmax,c(1,2),mean,na.rm=TRUE)
TMIN[,,y] = apply(tmin,c(1,2),mean,na.rm=TRUE)

PR[,,y] = apply(pr,c(1,2),sum,na.rm=TRUE)
TMAX95[,,y] = apply(tmax95,c(1,2),sum,na.rm=TRUE)
TMIN32[,,y] = apply(tmin32,c(1,2),sum,na.rm=TRUE)
PR50[,,y] = apply(pr50,c(1,2),sum,na.rm=TRUE)

message("Finished Calcs for ",histmonths[y])

}

}

if(obsdat=="METDATA"){
  filelisttmax = system(paste("ls ",filepath,obsdat,"/tasmax/*.nc",sep=""),intern=TRUE)
  filelisttmin = system(paste("ls ",filepath,obsdat,"/tasmin/*.nc",sep=""),intern=TRUE)
  filelistpr = system(paste("ls ",filepath,obsdat,"/pr/*.nc",sep=""),intern=TRUE)
  
  split1 = do.call("rbind",strsplit(filelisttmax,"_",fixed=TRUE))
  years = as.numeric(substr(split1[,2],1,4))
  yearidx = which(years>=histperiod[1] & years<= histperiod[2])
  
  for(y in 1:length(yearidx)){
    
    nctest = nc_open(filelisttmax[yearidx[1]])
    
    if(y==1){
      lon = ncvar_get(nctest,"lon")
      lat = ncvar_get(nctest,"lat")
      
      if(lon[1]<0) lonrange=lonrange-360
      lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
      latidx = which(lat>= latrange[1] & lat<=latrange[2])
      
      TMAX = TMIN = PR = TMAX95 = TMIN32 = PR50 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
      
    }
    
    tmax = ncvar_get(nctest,"air_temperature",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))-273.15
    nc_close(nctest)
    
    nctest = nc_open(filelisttmin[yearidx[1]])
    tmin = ncvar_get(nctest,"air_temperature",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))-273.15
    nc_close(nctest)
    
    nctest = nc_open(filelistpr[yearidx[1]])
    pr = ncvar_get(nctest,"precipitation_amount",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))
    nc_close(nctest)
    
    tmax95 = ifelse(tmax>=35,1,0)
    tmin32 = ifelse(tmin<=0,1,0)
    pr50 = ifelse(pr>=50.4,1,0)
    
    dates = histdates[which(as.numeric(substr(histdates,1,4))==years[yearidx[y]])]
    if(length(dates)>dim(tmin)[3]) dates = dates[-which(substr(dates,6,10)=="02-29")]
    
    tmaxmon = tminmon = prmon = tmax95mon = tmin32mon = pr50mon = array(NA,dim=c(length(lonidx),length(latidx),12))
    for(month in 1:12){
      monidx = which(as.numeric(substr(dates,6,7))==month)
      outidx = month+(12*(y-1))
      
      TMAX[,,month] = apply(tmax[,,monidx],c(1,2),mean,na.rm=TRUE)
      TMIN[,,month] = apply(tmin[,,monidx],c(1,2),mean,na.rm=TRUE)
      
      PR[,,month] = apply(pr[,,monidx],c(1,2),sum,na.rm=TRUE)
      TMAX95[,,month] = apply(tmax95[,,monidx],c(1,2),sum,na.rm=TRUE)
      TMIN32[,,month] = apply(tmin32[,,monidx],c(1,2),sum,na.rm=TRUE)
      PR50[,,month] = apply(pr50[,,monidx],c(1,2),sum,na.rm=TRUE)
    }

    message("Finished Calcs for ",years[yearidx[y]])
    
  }
  
}

if(obsdat=="DAYMET"){
  filelisttmax = system(paste("ls ",filepath,obsdat,"/na/tasmax/*.nc4",sep=""),intern=TRUE)
  filelisttmin = system(paste("ls ",filepath,obsdat,"/na/tasmin/*.nc4",sep=""),intern=TRUE)
  filelistpr = system(paste("ls ",filepath,obsdat,"/na/pr/*.nc4",sep=""),intern=TRUE)
  
  split1 = do.call("rbind",strsplit(filelisttmax,"_",fixed=TRUE))
  years = as.numeric(substr(split1[,4],1,4))
  #
  
  for(h in 1:length(histmonths)){
    
    yearidx = which(years== as.numeric(substr(histmonths[h],1,4)))
    dates = histdates[which(as.numeric(substr(histdates,1,4))==years[yearidx])]
    
    nctest = nc_open(filelisttmax[yearidx[1]])
    
    if(h==1){
      lon = ncvar_get(nctest,"lon")
      lat = ncvar_get(nctest,"lat")
      times = ncvar_get(nctest,"time")
      if(lon[1]<0) lonrange=lonrange-360
      lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
      latidx = which(lat>= latrange[1] & lat<=latrange[2])
      TMAX = TMIN = PR = TMAX95 = TMIN32 = PR50 = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
    }
    
    if(length(dates)>length(times)) dates = dates[-which(substr(dates,6,10)=="02-29")]
    dateidx = which(as.numeric(substr(dates,6,7))==as.numeric(substr(histmonths[h],6,7)))
    
    tmax = ncvar_get(nctest,"tmax",start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    nc_close(nctest)
    
    nctest = nc_open(filelisttmin[yearidx[1]])
    tmin = ncvar_get(nctest,"tmin",start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    nc_close(nctest)
    
    nctest = nc_open(filelistpr[yearidx[1]])
    pr = ncvar_get(nctest,"prcp",start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    nc_close(nctest)
    
    tmax95 = ifelse(tmax>=35,1,0)
    tmin32 = ifelse(tmin<=0,1,0)
    pr50 = ifelse(pr>=50.4,1,0)
    
      TMAX[,,h] = apply(tmax,c(1,2),mean,na.rm=TRUE)
      TMIN[,,h] = apply(tmin,c(1,2),mean,na.rm=TRUE)
      
      PR[,,h] = apply(pr,c(1,2),sum,na.rm=TRUE)
      TMAX95[,,h] = apply(tmax95,c(1,2),sum,na.rm=TRUE)
      TMIN32[,,h] = apply(tmin32,c(1,2),sum,na.rm=TRUE)
      PR50[,,h] = apply(pr50,c(1,2),sum,na.rm=TRUE)
    
    message("Finished Calcs for ",histmonths[h])
    
  }
  
}

##########
# monthly climo calcs
tmaxclimo = tminclimo = prclimo = tmax95climo = tmin32climo = pr50climo = array(NA,dim=c(length(lonidx),length(latidx),12))

for(month in 1:12){
  
  monthidx = which(as.numeric(substr(histmonths,6,7))==month)
  tmaxclimo[,,month] = apply(TMAX[,,monthidx],c(1,2),mean,na.rm=TRUE)
  tminclimo[,,month] = apply(TMIN[,,monthidx],c(1,2),mean,na.rm=TRUE)
  prclimo[,,month] = apply(PR[,,monthidx],c(1,2),mean,na.rm=TRUE)
  
  tmax95climo[,,month] = apply(TMAX95[,,monthidx],c(1,2),mean,na.rm=TRUE)
  tmin32climo[,,month] = apply(TMIN32[,,monthidx],c(1,2),mean,na.rm=TRUE)
  pr50climo[,,month] = apply(PR50[,,monthidx],c(1,2),mean,na.rm=TRUE)
  
  prclimo[,,month] = ifelse(is.na(tmaxclimo[,,month])==FALSE,prclimo[,,month],NA)
  pr50climo[,,month] = ifelse(is.na(tmaxclimo[,,month])==FALSE,pr50climo[,,month],NA)
  tmax95climo[,,month] = ifelse(is.na(tmaxclimo[,,month])==FALSE,tmax95climo[,,month],NA)
  tmin32climo[,,month] = ifelse(is.na(tmaxclimo[,,month])==FALSE,tmin32climo[,,month],NA)
  
  message("Finished Calcs for month ", month)
}

tmaxvar = apply(TMAX,c(1,2),var,na.rm=TRUE)
tminvar = apply(TMIN,c(1,2),var,na.rm=TRUE)
prvar = apply(PR,c(1,2),var,na.rm=TRUE)

tmax95var = apply(TMAX95,c(1,2),var,na.rm=TRUE)
tmin32var = apply(TMIN32,c(1,2),var,na.rm=TRUE)
pr50var = apply(PR50,c(1,2),var,na.rm=TRUE)

tmaxvar=ifelse(is.na(tmaxclimo[,,1])==FALSE,tmaxvar,NA)
tminvar=ifelse(is.na(tmaxclimo[,,1])==FALSE,tminvar,NA)
prvar=ifelse(is.na(tmaxclimo[,,1])==FALSE,prvar,NA)
tmax95var=ifelse(is.na(tmaxclimo[,,1])==FALSE,tmax95var,NA)
tmin32var=ifelse(is.na(tmaxclimo[,,1])==FALSE,tmin32var,NA)
pr50var=ifelse(is.na(tmaxclimo[,,1])==FALSE,pr50var,NA)

###########
# variance averages

tmaxvaravg=mean(tmaxvar,na.rm=TRUE)
tminvaravg=mean(tminvar,na.rm=TRUE)
prvaravg=mean(prvar,na.rm=TRUE)
tmax95varavg=mean(tmax95var,na.rm=TRUE)
tmin32varavg=mean(tmin32var,na.rm=TRUE)
pr50varavg=mean(pr50var,na.rm=TRUE)

save(list=c("tmaxvaravg","tminvaravg","prvaravg","tmax95varavg","tmin32varavg","pr50varavg"),file=paste(obsdat,"_variances.Rdata",sep=""))


###########
# write netcdfs out

dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
dimY <- ncdim_def( "lat", "degrees_north", lat[latidx])
dimT <- ncdim_def("month","month",1:12)

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

var1d <- ncvar_def("tmaxclimo","degrees_C", list(dimX,dimY,dimT), mv )
var2d <- ncvar_def("tminclimo","degrees_C", list(dimX,dimY,dimT), mv )
var3d <- ncvar_def("prclimo","mm", list(dimX,dimY,dimT), mv )
var4d <- ncvar_def("tmax95climo","degrees_C", list(dimX,dimY,dimT), mv )
var5d <- ncvar_def("tmin32climo","degrees_C", list(dimX,dimY,dimT), mv )
var6d <- ncvar_def("pr50climo","mm", list(dimX,dimY,dimT), mv )

#######
# Create netcdf file

nc <- nc_create(paste("/home/woot0002/",obsdat,"_histclimo.nc",sep="") ,  list(var1d,var2d,var3d,var4d,var5d,var6d) )

# Write some data to the file
ncvar_put(nc, var1d, tmaxclimo) # no start or count: write all values\
ncvar_put(nc, var2d, tminclimo)
ncvar_put(nc, var3d, prclimo)
ncvar_put(nc, var4d, tmax95climo) # no start or count: write all values\
ncvar_put(nc, var5d, tmin32climo)
ncvar_put(nc, var6d, pr50climo)

# close ncdf
nc_close(nc)


