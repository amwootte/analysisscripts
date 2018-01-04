###############
# 3^5 Projected Change Analyses - location information grab

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##############
# User supplied inputs

location_name = "Norman"
location_lat = 35.2226
location_lon = -97.4395

varname = "tasmax" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname # don't change this
#noleap=TRUE # if your data has leap days change this to FALSE, otherwise leave it alone.  Should be false for 3^5.

if(varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd") varin="pr"

difftype="absolute" # type of difference to take, can be either absolute or percent

appfunc = "mean" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

TC = FALSE  # Threshold calculator - should this be calculating a threshold? TRUE (calculate threshold) or FALSE(don't calculate threshold)
TH = 270.928  # Threshold value - what's the threshold the script should calculate for?
cond = "lte" # Threshold condition - "gte" = greater than or equal to, "lte" = less than or equal to, "gt" = greater than, "lt"= less than
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

####
# get model info.

test = nc_open(histfilelist[1])
lat = ncvar_get(test,"lat")
lon = ncvar_get(test,"lon")
times = ncvar_get(test,"time")
startdate = substr(test$dim[[4]]$units,12,21)
nc_close(test)

###
# create model grid

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
modelgrid$lon = modelgrid$lon-360

###
# get cells to use

pointarea = distfunc(location_lon,location_lat,modelgrid)
locstart = c(pointarea$R-1,pointarea$C-1)
loccount = c(3,3)

###

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

histlist = c()

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
  
  histlist[i] = netcdftopointtoclimocalcs(histfilelist[i],varname=varin,dimnames=c("lon","lat","time"),locstart=locstart,loccount=loccount,threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appliedfunction=appfunc)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

########
# Future data grab

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
#if(noleap==TRUE) dates = dates[-which(substr(dates,6,10)=="02-29")]

projlist = c()

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  }else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  projlist[i] = netcdftopointtoclimocalcs(projfilelist[i],varname=varin,dimnames=c("lon","lat","time"),locstart=locstart,loccount=loccount,threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(2041,2070),datesdataperiod=datesin,appliedfunction=appfunc)
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

######
# Difference Calcs

diffs = c()
for(i in 1:length(projfilelist)){
  GCMin = projfilebreakdown$GCM[i]
  obsin = projfilebreakdown$obs[i]
  histidx = which(histfilebreakdown$GCM==GCMin & histfilebreakdown$obs==obsin)
  if(difftype=="absolute") diffs[i]=projlist[i]-histlist[histidx]
  if(difftype=="percent") diffs[i]= ((projlist[i]-histlist[histidx])/histlist[histidx])*100
}

#####
# write out data files

modhist = c(rep(histlist[1:3],3),rep(histlist[4:6],3),rep(histlist[7:9],3))
projchangedat = projfilebreakdown

projchangedat$modhist = modhist
projchangedat$modfut = projlist
projchangedat$projchange = diffs

projchangedat = projchangedat[,c(3,5,7:11)]

filename = paste(varname,"_",location_name,"_",difftype,".csv",sep="")

write.table(projchangedat,file=filename,row.names=FALSE,sep=",")


