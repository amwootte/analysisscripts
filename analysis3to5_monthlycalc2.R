###############
# 3^5 Projected Change Analyses - monthly converter

library(ncdf4) # loading necessary libraries and extra functions
source("analysisfunctions.R")

##############
# User supplied inputs

varname = "pr" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname
dataunits = "mm"
appfunc = "sum" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/EDQM/*historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/EDQM/*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),9)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),3)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  yearmon = unique(substr(datesin,1,7))
  
  for(y in 1:length(yearmon)){
    yearidx = which(substr(datesin,1,7)==yearmon[y])
    
    test = nc_open(histfilelist[i])
    tempdata = ncvar_get(test,varname,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    if(varname=="pr") tempdata=tempdata*86400
    
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata)[1],dim(tempdata)[2],length(yearmon)))
      lat = ncvar_get(test,"lat")
      lon = ncvar_get(test,"lon")
      times = ncvar_get(test,"time")
      domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
      startdate = substr(test$dim[[4]]$units,12,21)
    }
    
    #if(any(yearidx>length(times))==TRUE){ yearidx = yearidx[1]:length(times) }
     nc_close(test)
    
    if(appfunc=="mean") yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
    if(appfunc=="sum") {
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    rm(tempdata)
    message("Finished Calcs for yearmon ",yearmon[y])
  }
  
  namesplit = do.call("c",strsplit(histfilelist[i],"_"))
  
  filename = paste(varname,"_mon_",namesplit[3],"_",namesplit[4],"_",namesplit[5],"_",namesplit[6],"_",namesplit[7],sep="")
    
  dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  dimT <- ncdim_def( "time", paste("months since ",substr(yearmon[1],1,4),"-01-15",sep=""),0:(length(yearmon)-1))
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  var1d <- ncvar_def(varname,dataunits,longname=paste("Monthly ",appfunc," of ",varname,sep=""), list(dimX,dimY,dimT), mv )
  nc <- nc_create(paste("/data2/3to5/I35/",varname,"/EDQM/monthly/",filename,sep="") ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, yearlydat) # no start or count: write all values\
  nc_close(nc)
  rm(yearlydat)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

#######
# monthly values on projected values

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  message("Starting work on file ",projfilelist[i])
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  yearmon = unique(substr(datesin,1,7))
  
  for(y in 1:length(yearmon)){
    yearidx = which(substr(datesin,1,7)==yearmon[y])
    
    test = nc_open(projfilelist[i])
    tempdata = ncvar_get(test,varname,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    if(varname=="pr") tempdata=tempdata*86400
    
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata)[1],dim(tempdata)[2],length(yearmon)))
      lat = ncvar_get(test,"lat")
      lon = ncvar_get(test,"lon")
      times = ncvar_get(test,"time")
      domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
      startdate = substr(test$dim[[4]]$units,12,21)
    }
    nc_close(test)
    
    if(appfunc=="mean") yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
    if(appfunc=="sum") {
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    rm(tempdata)
    message("Finished Calcs for yearmon ",yearmon[y])
  }
  
  namesplit = do.call("c",strsplit(projfilelist[i],"_"))
  
  filename = paste(varname,"_mon_",namesplit[3],"_",namesplit[4],"_",namesplit[5],"_",namesplit[6],"_",namesplit[7],sep="")
  
    
    dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  dimT <- ncdim_def( "time", paste("months since ",substr(yearmon[1],1,4),"-01-15",sep=""),0:(length(yearmon)-1))
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  var1d <- ncvar_def(varname,dataunits,longname=paste("Monthly ",appfunc," of ",varname,sep=""), list(dimX,dimY,dimT), mv )
  nc <- nc_create(paste("/data2/3to5/I35/",varname,"/EDQM/monthly/",filename,sep="") ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, yearlydat) # no start or count: write all values\
  nc_close(nc)
  rm(yearlydat)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}


