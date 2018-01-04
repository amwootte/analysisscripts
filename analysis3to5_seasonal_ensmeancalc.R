###############
# 3^5 Projected Change Analyses

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##############
# User supplied inputs

varname = "tasmin" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname # don't change this
#noleap=TRUE # if your data has leap days change this to FALSE, otherwise leave it alone.  Should be false for 3^5.

if(varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28") varin="tasmin"
if(varname=="pr25" | varname=="pr50" | varname=="mdrn") varin="pr"

if(varname=="tasmax" | varname=="tasmin") histunits = projunits = "degrees K" # set units
if(varname=="pr" | varname=="rx1day" | varname=="rx5day") histunits = projunits = "mm"
if(varname=="tmax95" | varname=="tmax100" | varname=="tmin32" | varname=="tmin28" | varname=="cdd" | varname=="cwd" | varname=="frd" | varname=="pr25" | varname=="pr50" | varname=="mdrn") histunits=projunits="days"
if(difftype=="absolute") diffunits=histunits
if(difftype=="percent") diffunits="%"

difftype="absolute" # type of difference to take, can be either absolute or percent
seasonin = "DJF" # which season should be calculated

applymask = NA # if you want to apply a state mask, leave this alone

appfunc = "mean" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

TC = FALSE   # Threshold calculator - should this be calculating a threshold? TRUE (calculate threshold) or FALSE(don't calculate threshold)
TH = 308.15  # Threshold value - what's the threshold the script should calculate for?
cond = "gte" # Threshold condition - "gte" = greater than or equal to, "lte" = less than or equal to, "gt" = greater than, "lt"= less than
# threshold value and threshold condition are ignored if TC=FALSE
# should be using these to calculate tmax95 and the others, temperature thresholds should be supplied in degrees K, precipitation thresholds in mm
# thresholds for each of the following are
# tmax95: 308.15 - tasmax greater than or equal to 95F
# tmax100: 310.928 - tasmax greater than or equal to 100F
# tmin32: 273.15 - tasmin less than or equal to 32F
# tmin28: 270.928 - tasmin less than or equal to 28F
# pr25: 25.4 - precipitation greater than or equal to 1 inch
# pr50: 50.8 - precipitation greater than or equal to 2 inches
# mdrn: 0.254 - precipitation greater than or equal to trace

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
  
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  yearlyoutput = netcdftoseasonalcalcs(histfilelist[i],varname=varin,season=seasonin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appliedfunction=appfunc)
  
  if(i==1){
    lon = yearlyoutput[[1]]
    lat = yearlyoutput[[2]]
    histlist = array(NA,dim=c(length(lon),length(lat),length(histfilelist)))
  }
  
  histlist[,,i] = climocalc(yearlyoutput[[3]],yearlydataperiod=c(1981,2005),climoperiod=c(1981,2005))
  rm(yearlyoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

########
# Future data grab

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
#if(noleap==TRUE) dates = dates[-which(substr(dates,6,10)=="02-29")]

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  #if(noleap==TRUE) dates = dates[-which(substr(dates,6,10)=="02-29")]
  
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  }else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  yearlyoutput = netcdftoseasonalcalcs(projfilelist[i],varname=varin,season=seasonin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(2041,2070),datesdataperiod=datesin,appliedfunction=appfunc)
  
  if(i==1){
    projlist = array(NA,dim=c(length(lon),length(lat),length(projfilelist)))
  }
  
  projlist[,,i] = climocalc(yearlyoutput[[3]],yearlydataperiod=c(2041,2070),climoperiod=c(2041,2070))
  rm(yearlyoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

######
# Difference Calcs

diffs = array(NA,dim=dim(projlist))
for(i in 1:length(projfilelist)){
  GCMin = projfilebreakdown$GCM[i]
  obsin = projfilebreakdown$obs[i]
  histidx = which(histfilebreakdown$GCM==GCMin & histfilebreakdown$obs==obsin)
  diffs[,,i]=diffcalc(projlist[,,i],histlist[,,histidx],type=difftype)
}

#####
# Group by Emissions Scenario

diffsg1 = array(NA,dim=c(length(lon),length(lat),length(unique(projfilebreakdown$scen))))
scens = unique(projfilebreakdown$scen)
for(s in 1:length(scens)){
  scenidx = which(projfilebreakdown$scen==scens[s])
  diffsg1[,,s] = apply(diffs[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

######
# Create Ensemble means netcdf

filename = paste(varname,"_ensmean_",seasonin,"_",difftype,"_",futureperiod[1],"-",futureperiod[2],".nc",sep="")

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat)

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

var1d <- ncvar_def("histmean",histunits,longname="Historical Mean", list(dimX,dimY), mv )

var2d <- ncvar_def("projmean_rcp26",projunits,longname="Projected Mean RCP 2.6", list(dimX,dimY), mv )
var3d <- ncvar_def("projmean_rcp45",projunits,longname="Projected Mean RCP 4.5", list(dimX,dimY), mv )
var4d <- ncvar_def("projmean_rcp85",projunits,longname="Projected Mean RCP 8.5", list(dimX,dimY), mv )

var5d <- ncvar_def("projmeandiff_rcp26",diffunits,longname="Mean Projected Change RCP 2.6", list(dimX,dimY), mv )
var6d <- ncvar_def("projmeandiff_rcp45",diffunits,longname="Mean Projected Change RCP 4.5", list(dimX,dimY), mv )
var7d <- ncvar_def("projmeandiff_rcp85",diffunits,longname="Mean Projected Change RCP 8.5", list(dimX,dimY), mv )

#######
# Create netcdf file

nc <- nc_create(paste("/data2/3to5/I35/ens_means/seasonal/",filename,sep="") ,  list(var1d,var2d,var3d,var4d,var5d,var6d,var7d) )

# Write some data to the file
ncvar_put(nc, var1d, histsg1) # no start or count: write all values\

ncvar_put(nc, var2d, projsg1[,,1])
ncvar_put(nc, var3d, projsg1[,,2])
ncvar_put(nc, var4d, projsg1[,,3])

ncvar_put(nc, var5d, diffsg1[,,1])
ncvar_put(nc, var6d, diffsg1[,,2])
ncvar_put(nc, var7d, diffsg1[,,3])

# close ncdf
nc_close(nc)

