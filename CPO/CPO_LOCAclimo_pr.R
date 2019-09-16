####################
#
# CPO Analysis
# Monthly Historical LOCA Climatology

library(ncdf4)
library(maps)
library(fields)
library(sp)

histperiod = c(1981,2005)
histdates = seq(as.Date("1950-01-01"),as.Date("2005-12-31"),by="day")
histmonths = seq(as.Date("1981-01-15"),as.Date("2005-12-15"),by="month")
latrange = c(26,40)
lonrange = c(251,270)

filepath = "/data4/data/DS_proj/LOCA/"

###########
#######
#  pr all GCMs

filelistpr1 = system(paste("ls ",filepath,"pr/historical/*historical.nc",sep=""),intern=TRUE)
split3 = do.call("rbind",strsplit(filelistpr1,"_",fixed=TRUE))

PRlist = list()
PR50list = list()
R1MMlist = list()
RX1DAYlist = list()
RX5DAYlist = list()

for(f in 1:length(filelistpr1)){
  nctest = nc_open(filelistpr1[f])
  times = ncvar_get(nctest,"time")
  
  if(length(times)==length(histdates)) dates=histdates
  if(length(times)<length(histdates)) dates = histdates[-which(substr(histdates,6,10)=="02-29")]
  
  if(f==1){
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    if(lon[1]<0 & lonrange[1]>0) lonrange=lonrange-360
    if(lon[1]>0 & lonrange[1]>0) lonrange=lonrange-360;lon=lon-360;
    lonidx = which(lon>= lonrange[1] & lon<=lonrange[2])
    latidx = which(lat>= latrange[1] & lat<=latrange[2])
  }
  nc_close(nctest)
  PR = PR50 = R1MM = RX1DAY = RX5DAY = array(NA,dim=c(length(lonidx),length(latidx),length(histmonths)))
  
  for(h in 1:length(histmonths)){
    dateidx = which(substr(dates,1,7)==substr(histmonths[h],1,7))
    nctest = nc_open(filelistpr1[f])
    pr = ncvar_get(nctest,paste("pr_",split3[f,3],"_",split3[f,4],"_historical",sep=""),start=c(lonidx[1],latidx[1],dateidx[1]),count=c(length(lonidx),length(latidx),length(dateidx)))
    if(max(pr,na.rm=TRUE)<1){
      pr = pr*86400
    }
    nc_close(nctest)
    
    pr50 = ifelse(pr>=50.4,1,0)
    r1mm = ifelse(pr>=1,1,0)
    
    PR[,,h] = apply(pr,c(1,2),sum,na.rm=TRUE)
    PR[,,h] = ifelse(is.na(pr[,,1])==FALSE,PR[,,h],NA)
    PR50[,,h] = apply(pr50,c(1,2),sum,na.rm=TRUE)
    PR50[,,h] = ifelse(is.na(pr[,,1])==FALSE,PR50[,,h],NA)
    R1MM[,,h] = apply(r1mm,c(1,2),sum,na.rm=TRUE)
    R1MM[,,h] = ifelse(is.na(pr[,,1])==FALSE,R1MM[,,h],NA)
    RX1DAY[,,h] = apply(pr,c(1,2),max,na.rm=TRUE)
    RX1DAY[,,h] = ifelse(is.na(pr[,,1])==FALSE,RX1DAY[,,h],NA)
    RX5DAY[,,h] = apply(pr,c(1,2),calcrollsum,size=5)
    RX5DAY[,,h] = ifelse(is.na(pr[,,1])==FALSE,RX5DAY[,,h],NA)
    message("Finished calcs for histmonth ",histmonths[h]," and file ",f)
  }
  
  prclimo = pr50climo = r1mmclimo = rx1dayclimo = rx5dayclimo = array(NA,dim=c(length(lonidx),length(latidx),12))
  for(month in 1:12){
    monidx = which(as.numeric(substr(histmonths,6,7))==month)
    prclimo[,,month] = apply(PR[,,monidx],c(1,2),mean,na.rm=TRUE)
    pr50climo[,,month] = apply(PR50[,,monidx],c(1,2),mean,na.rm=TRUE)
    r1mmclimo[,,month] = apply(R1MM[,,monidx],c(1,2),mean,na.rm=TRUE)
    rx1dayclimo[,,month] = apply(RX1DAY[,,monidx],c(1,2),mean,na.rm=TRUE)
    rx5dayclimo[,,month] = apply(RX5DAY[,,monidx],c(1,2),mean,na.rm=TRUE)
    message("Finished climo calcs for month ",month," and file ",f)
  }
  
  PRlist[[f]] = prclimo
  PR50list[[f]] = pr50climo
  R1MMlist[[f]] = r1mmclimo
  RX1DAYlist[[f]] = rx1dayclimo
  RX5DAYlist[[f]] = rx5dayclimo
  
  rm(PR)
  rm(PR50)
  rm(R1MM)
  rm(RX1DAY)
  rm(RX5DAY)
  rm(prclimo)
  rm(pr50climo)
  rm(r1mmclimo)
  rm(rx1dayclimo)
  rm(rx5dayclimo)
  
  message("Finished calcs for file ",filelistpr1[f])
  
}

###############
# netcdf creation

dimX <- ncdim_def( "lon", "degrees_east", lon[lonidx])
dimY <- ncdim_def( "lat", "degrees_north", lat[latidx])
dimT <- ncdim_def("month","month",1:12)
# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

for(i in 1:nrow(split3)){
  
  var1d <- ncvar_def("prclimo","mm", list(dimX,dimY,dimT), mv )
  var2d <- ncvar_def("r1mmclimo","days", list(dimX,dimY,dimT), mv )
  var3d <- ncvar_def("pr50climo","days", list(dimX,dimY,dimT), mv )
  var4d <- ncvar_def("rx1dayclimo","mm", list(dimX,dimY,dimT), mv )
  var5d <- ncvar_def("rx5dayclimo","mm", list(dimX,dimY,dimT), mv )
  
  #######
  # Create netcdf file
  
  nc <- nc_create(paste("/home/woot0002/pr_",split3[i,3],"_",split3[i,4],"_LOCA_histclimo.nc",sep="") ,  list(var1d,var2d,var3d,var4d,var5d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, PRlist[[i]]) # no start or count: write all values\
  ncvar_put(nc, var2d, R1MMlist[[i]]) 
  ncvar_put(nc, var3d, PR50list[[i]]) # no start or count: write all values\
  ncvar_put(nc, var4d, RX1DAYlist[[i]]) 
  ncvar_put(nc, var5d, RX5DAYlist[[i]]) 
  # close ncdf
  nc_close(nc)
  
}
