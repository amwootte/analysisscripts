###############
# 3^5 Projected Change Analyses - monthly drought metrics

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
library(SPEI)
source("analysisfunctions.R")

##############
# User supplied inputs

varname = "SPI" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
scale = 1
#noleap=TRUE # if your data has leap days change this to FALSE, otherwise leave it alone.  Should be false for 3^5.

difftype="absolute" # type of difference to take, can be either absolute or percent

applymask = NA # if you don't want to apply a state mask, leave this alone

colorchoicediff = "redtoblue" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 30 # maximum number of color bins allowed for plotting the projected changes

appfunc = "spei" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

dataunits=""

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/pr/EDQM/monthly/*historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/pr/EDQM/monthly/*rcp*.nc",sep=""),intern=T)

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

if(appfunc=="spei"){
  tmaxhistlist = system(paste("ls /data2/3to5/I35/tasmax/EDQM/monthly/*historical*.nc",sep=""),intern=T)
  tmaxprojlist = system(paste("ls /data2/3to5/I35/tasmax/EDQM/monthly/*rcp*.nc",sep=""),intern=T)
  
  tminhistlist = system(paste("ls /data2/3to5/I35/tasmin/EDQM/monthly/*historical*.nc",sep=""),intern=T)
  tminprojlist = system(paste("ls /data2/3to5/I35/tasmin/EDQM/monthly/*rcp*.nc",sep=""),intern=T)
}

dates = seq(as.Date("1981-01-15"),as.Date("2005-12-15"),by="month")

dohist = TRUE

if(dohist==TRUE){

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
 
  datesin = dates
  
  if(varname=="SPI"){
    filenames=histfilelist[i]
    varnamesin="pr"
  }
  
  if(varname=="SPEI"){
    filenames=paste(histfilelist[i],tmaxhistlist[i],tminhistlist[i],sep=",")
    varnamesin=c("pr","tasmax","tasmin")
  }
  
  droughtoutput=netcdfdroughtcalcs(filenames=filenames,varnames=varnamesin,dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appfunc=appfunc,scale=1)
  
  if(i==1){
    lon = droughtoutput[[1]]
    lat = droughtoutput[[2]]
    times = droughtoutput[[3]]
  }
  droughtoutput[[4]] = ifelse(droughtoutput[[4]]==Inf | droughtoutput[[4]]== -Inf,NA,droughtoutput[[4]])
  
  namesplit = do.call("c",strsplit(histfilelist[i],"_"))
  
  filename = paste(varname,scale,"_mon_",namesplit[3],"_",namesplit[4],"_",namesplit[5],"_",namesplit[6],"_",namesplit[7],sep="")
  
  dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  dimT <- ncdim_def( "time", paste("months since ",substr(datesin[1],1,4),"-01-15",sep=""),0:(length(datesin)-1))
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  var1d <- ncvar_def(varname,dataunits,longname=paste(varname ,"scale = ",scale," months",sep=""), list(dimX,dimY,dimT), mv )
  nc <- nc_create(paste("/data2/3to5/I35/",varname,"/EDQM/",filename,sep="") ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, droughtoutput[[4]]) # no start or count: write all values\
  nc_close(nc)
  rm(droughtoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}
}

########
# Future data grab

dates = seq(as.Date("2006-01-15"),as.Date("2099-12-15"),by="month")

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  
  datesin = dates
  
  if(varname=="SPI"){
    filenames=projfilelist[i]
    varnamesin="pr"
  }
  
  if(varname=="SPEI"){
    filenames=paste(projfilelist[i],tmaxprojlist[i],tminprojlist[i],sep=",")
    varnamesin=c("pr","tasmax","tasmin")
  }
  
  droughtoutput=netcdfdroughtcalcs(filenames=filenames,varnames=varnamesin,dimnames=c("lon","lat","time"),yearlydataperiod=c(2006,2099),datesdataperiod=datesin,appfunc=appfunc,scale=1)
  
  if(i==1){
    lon = droughtoutput[[1]]
    lat = droughtoutput[[2]]
    times = droughtoutput[[3]]
  }
  
  droughtoutput[[4]] = ifelse(droughtoutput[[4]]==Inf | droughtoutput[[4]]== -Inf,NA,droughtoutput[[4]])
  
  namesplit = do.call("c",strsplit(projfilelist[i],"_"))
  
  filename = paste(varname,scale,"_mon_",namesplit[3],"_",namesplit[4],"_",namesplit[5],"_",namesplit[6],"_",namesplit[7],sep="")
  
  dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  dimT <- ncdim_def( "time", paste("months since ",substr(datesin[1],1,4),"-01-15",sep=""),0:(length(datesin)-1))
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  var1d <- ncvar_def(varname,dataunits,longname=paste(varname ," scale = ",scale," months",sep=""), list(dimX,dimY,dimT), mv )
  nc <- nc_create(paste("/data2/3to5/I35/",varname,"/EDQM/",filename,sep="") ,  var1d )
  # Write some data to the file
  
  ncvar_put(nc, var1d, droughtoutput[[4]]) # no start or count: write all values\
  nc_close(nc)
  rm(droughtoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

