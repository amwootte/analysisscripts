#######################
#
# Monthly climo calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = "tasmin"

files = system(paste("ls /data2/3to5/I35/",var,"/*/",var,"_day*historical*.nc",sep=""),intern=TRUE)

filelist = do.call("rbind",strsplit(files,"/",fixed=TRUE))

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
monyear = dates[which(substr(dates,9,10)=="01")]
years = 1981:2005
months = 1:12

if(var=="pr"){ 
  endpt = length(files)-1
} else {
  endpt = length(files)
}
  
for(i in 1:endpt){
  
  ptm = proc.time()
  
  test = nc_open(files[i])
  times = ncvar_get(test,"time")
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  
  if(length(dates)>length(times)){
    histdates = dates[-which(substr(dates,6,10)=="02-29")]
  } else {
    histdates = dates
  }
  
  tmp = tmp2 = array(NA,dim=c(length(lon),length(lat),length(monyear)))
  if(var=="pr"){
    tmp3 = tmp4 = tmp5 = array(NA,dim=c(length(lon),length(lat),length(monyear)))
  }
  
  for(j in 1:length(monyear)){
    idx = which(as.numeric(substr(histdates,1,4))==as.numeric(substr(monyear[j],1,4)) & as.numeric(substr(histdates,6,7))==as.numeric(substr(monyear[j],6,7)))
    test = nc_open(files[i])
    if(var=="pr"){
    prtmp = ncvar_get(test,var,start=c(1,1,idx[1]),count=c(-1,-1,length(idx)))*86400
    nc_close(test)
    tmp[,,j] = apply(prtmp,c(1,2),sum,na.rm=TRUE)
    tmp[,,j] = ifelse(is.na(prtmp[,,1])==FALSE,tmp[,,j],NA)
    pr50 = ifelse(prtmp>=50.4,1,0)
    r1mm = ifelse(prtmp>=1,1,0)
    
    tmp2[,,j] = apply(pr50,c(1,2),sum,na.rm=TRUE)
    tmp2[,,j] = ifelse(is.na(prtmp[,,1])==FALSE,tmp2[,,j],NA)
    
    tmp3[,,j] = apply(r1mm,c(1,2),sum,na.rm=TRUE)
    tmp3[,,j] = ifelse(is.na(prtmp[,,1])==FALSE,tmp3[,,j],NA)
    
    tmp4[,,j] = apply(prtmp,c(1,2),max,na.rm=TRUE)
    tmp4[,,j] = ifelse(is.na(prtmp[,,1])==FALSE,tmp4[,,j],NA)
    
    tmp5[,,j] = apply(prtmp,c(1,2),calcrollsum,size=5)
    tmp5[,,j] = ifelse(is.na(prtmp[,,1])==FALSE,tmp5[,,j],NA)
    
    } else {
      prtmp = ncvar_get(test,var,start=c(1,1,idx[1]),count=c(-1,-1,length(idx)))-273.15
      nc_close(test)
      tmp[,,j] = apply(prtmp,c(1,2),mean,na.rm=TRUE)
      if(var=="tasmax"){
        prtmp2 = ifelse(prtmp>=35,1,0)
      }
      if(var=="tasmin"){
        prtmp2 = ifelse(prtmp<=0,1,0)
      }
      tmp2[,,j] = apply(prtmp2,c(1,2),sum,na.rm=TRUE)
      tmp2[,,j] = ifelse(is.na(prtmp[,,1])==FALSE,tmp2[,,j],NA)
    }
    message("Calcs complete for month-year ",monyear[j])
  }
  
  var1mon = var2mon = array(NA,dim=c(length(lon),length(lat),12))
  if(var=="pr"){
    var3mon = var4mon = var5mon = array(NA,dim=c(length(lon),length(lat),12))
  } 
  for(m in 1:12){
    idx = which(as.numeric(substr(monyear,6,7))==m)
    var1mon[,,m] = apply(tmp[,,idx],c(1,2),mean,na.rm=TRUE)
    var2mon[,,m] = apply(tmp2[,,idx],c(1,2),mean,na.rm=TRUE)
    if(var=="pr"){
      var3mon[,,m] = apply(tmp3[,,idx],c(1,2),mean,na.rm=TRUE)
      var4mon[,,m] = apply(tmp4[,,idx],c(1,2),mean,na.rm=TRUE)
      var5mon[,,m] = apply(tmp5[,,idx],c(1,2),mean,na.rm=TRUE)
    }
  }

  newfilename = paste(substr(filelist[i,7],1,nchar(filelist[i,7])-3),"_mon.nc",sep="")
  
  dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  dimT <- ncdim_def("time","month",1:12)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  
  if(var=="pr"){
    units1="mm"
    lname1 = "Monthly Precipitation Climo"
    varname1 = "prclimo"
    units2="days"
    lname2 = "Monthly Heavy Rain Days (>=50.4mm) Climo"
    varname2 = "pr50climo"
    units3="days"
    lname3 = "Monthly Rain Days (>=1mm) Climo"
    varname3 = "r1mmclimo"
    units4="mm"
    lname4 = "Monthly Max 1-day Rainfall Climo"
    varname4 = "rx1dayclimo"
    units5="mm"
    lname5 = "Monthly Max 5-day Rainfall Climo"
    varname5 = "rx5dayclimo"
  } 
  
  if(var=="tasmax"){
    units1="C"
    lname1 = "Monthly High Temperature Climo"
    varname1 = "tmaxclimo"
    units2="days"
    lname2 = "Monthly Hot Days (tmax>=95F) Climo"
    varname2 = "tmax95climo"
  }
  
  if(var=="tasmin"){
    units1="C"
    lname1 = "Monthly Low Temperature Climo"
    varname1 = "tminclimo"
    units2="days"
    lname2 = "Monthly Cold Days (tmin<=32F) Climo"
    varname2 = "tmin32climo"
  }
  
  var1d <- ncvar_def(varname1,units=units1,longname=lname1, list(dimX,dimY,dimT), mv )
  var2d <- ncvar_def(varname2,units=units2,longname=lname2, list(dimX,dimY,dimT), mv )
  if(var=="pr"){
    var3d <- ncvar_def(varname3,units=units3,longname=lname3, list(dimX,dimY,dimT), mv )
    var4d <- ncvar_def(varname4,units=units4,longname=lname4, list(dimX,dimY,dimT), mv )
    var5d <- ncvar_def(varname5,units=units5,longname=lname5, list(dimX,dimY,dimT), mv )
  }
  
  #######
  # Create netcdf file
  if(var!="pr"){
    nc <- nc_create(paste("/home/woot0002/monthlyclimo/",newfilename,sep="") ,  list(var1d,var2d) )
    # Write some data to the file
    ncvar_put(nc, var1d, var1mon) # no start or count: write all values\
    ncvar_put(nc, var2d, var2mon) # no start or count: write all values\
    # close ncdf
    nc_close(nc)
  } else {
    nc <- nc_create(paste("/home/woot0002/monthlyclimo/",newfilename,sep="") ,  list(var1d,var2d,var3d,var4d,var5d) )
    # Write some data to the file
    ncvar_put(nc, var1d, var1mon) # no start or count: write all values\
    ncvar_put(nc, var2d, var2mon) # no start or count: write all values\
    ncvar_put(nc, var3d, var3mon) # no start or count: write all values\
    ncvar_put(nc, var4d, var4mon) # no start or count: write all values\
    ncvar_put(nc, var5d, var5mon) # no start or count: write all values\
    # close ncdf
    nc_close(nc)
  }
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(files))
  message("Process time: ",ptmend[3]," secs")
  
}


