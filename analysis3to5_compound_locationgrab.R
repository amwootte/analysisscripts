###############
# 3^5 Projected Change Analyses - annual only

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##############
# User supplied inputs

location_name = "OKC"
location_lat = 35.4676
location_lon = -97.5164

varname = "gsl" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50

difftype="absolute" # type of difference to take, can be either absolute or percent

applymask = NA # if you don't want to apply a state mask, leave this alone

colorchoicediff = "bluetored" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 20 # maximum number of color bins allowed for plotting the projected changes

appfunc = "growing_season_length" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

###########
# 1. Data Gather and conversion

if(varname=="gsl" | varname=="heatwaves"){
histfilelist = system("ls /data2/3to5/I35/tasmax/EDQM/*historical*.nc",intern=T)
projfilelist = system("ls /data2/3to5/I35/tasmax/EDQM/*rcp*.nc",intern=T)

tmaxhistfilelist = histfilelist
tmaxprojfilelist = projfilelist

tminhistfilelist = system("ls /data2/3to5/I35/tasmin/EDQM/*historical*.nc",intern=T)
tminprojfilelist = system("ls /data2/3to5/I35/tasmin/EDQM/*rcp*.nc",intern=T)
} 

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

if(varname=="heatwaves"){
  tmaxq95filelist = system("ls tasmax*historical*_q95.nc",intern=T)
  tminq95filelist = system("ls tasmin*historical*_q95.nc",intern=T)
}

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

########

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

histlist=c()

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
  
  if(varname=="gsl"){
    filenames=paste(tmaxhistfilelist[i],tminhistfilelist[i],sep=",")
    varnamesin=c("tasmax","tasmin")
  }
  
  if(varname=="heatwaves"){
    filenames=paste(tmaxhistfilelist[i],tminhistfilelist[i],tmaxq95filelist[i],tminq95filelist[i],sep=",")
    varnamesin=c("tasmax","tasmin","tmaxq95","tminq95")
  }
  
  histlist[i]=netcdftopointcomboclimo(filenames=filenames,varnames=varnamesin,dimnames=c("lon","lat","time"),locstart=locstart,loccount=loccount,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,combofunction=appfunc)
  
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
  #if(noleap==TRUE) dates = dates[-which(substr(dates,6,10)=="02-29")]
  
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  }else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  
  if(varname=="gsl"){
    filenames=paste(tmaxprojfilelist[i],tminprojfilelist[i],sep=",")
    varnamesin=c("tasmax","tasmin")
  }
  
  if(varname=="heatwaves"){
    GCMin = projfilebreakdown$GCM[i]
    obsin = projfilebreakdown$obs[i]
    histidx = which(histfilebreakdown$GCM==GCMin & histfilebreakdown$obs==obsin)
    filenames=paste(tmaxprojfilelist[i],tminprojfilelist[i],tmaxq95filelist[histidx],tminq95filelist[histidx],sep=",")
    varnamesin=c("tasmax","tasmin","tmaxq95","tminq95")
  }
  
  projlist[i]=netcdftopointcomboclimo(filenames=filenames,varnames=varnamesin,dimnames=c("lon","lat","time"),locstart=locstart,loccount=loccount,yearlydataperiod=c(2041,2070),datesdataperiod=datesin,combofunction=appfunc)
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

#####
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

