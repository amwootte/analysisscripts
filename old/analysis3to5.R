###############
# 3^5 Projected Change Analyses - annual only

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##############
# User supplied inputs

varname = "pr50" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname # don't change this
#noleap=TRUE # if your data has leap days change this to FALSE, otherwise leave it alone.  Should be false for 3^5.

if(varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd") varin="pr"

difftype="absolute" # type of difference to take, can be either absolute or percent

applymask = NA # if you don't want to apply a state mask, leave this alone

colorchoicediff = "whitetored" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 30 # maximum number of color bins allowed for plotting the projected changes

appfunc = "sum" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

TC = TRUE   # Threshold calculator - should this be calculating a threshold? TRUE (calculate threshold) or FALSE(don't calculate threshold)
TH = 50.8  # Threshold value - what's the threshold the script should calculate for?
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
  message("Starting work on file ",histfilelist[i])
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  yearlyoutput = netcdftoyearlycalcs(histfilelist[i],varname=varin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appliedfunction=appfunc)
    
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
  
  yearlyoutput = netcdftoyearlycalcs(projfilelist[i],varname=varin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(2071,2099),datesdataperiod=datesin,appliedfunction=appfunc)
  
  if(i==1){
    projlist = array(NA,dim=c(length(lon),length(lat),length(projfilelist)))
  }
  
  projlist[,,i] = climocalc(yearlyoutput[[3]],yearlydataperiod=c(2071,2099),climoperiod=c(2071,2099))
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

################
lon=lon-360

if(is.na(applymask)==FALSE){
diffs = statemask(diffs,inputlat=lat,inputlon=lon,state=applymask)
diffsg1 = statemask(diffsg1,inputlat=lat,inputlon=lon,state=applymask)
projlist = statemask(projlist,inputlat=lat,inputlon=lon,state=applymask)
histlist = statemask(histlist,inputlat=lat,inputlon=lon,state=applymask)
} else {
  diffs = list(lon=lon,lat=lat,outputdata=diffs)
  diffsg1 = list(lon=lon,lat=lat,outputdata=diffsg1)
  projlist = list(lon=lon,lat=lat,outputdata=projlist)
  histlist = list(lon=lon,lat=lat,outputdata=histlist)
}

################
# Plotting

diffcolorbar = colorramp(diffs[[3]],colorchoice="browntogreen",Blimit=BINLIMIT,type="difference")

diffs_sort = diffs[[3]][,,order(projfilebreakdown$scen)]
projfilebreakdown = projfilebreakdown[order(projfilebreakdown$scen),]

pdf(paste("IndividualMembers_",varname,"_",difftype,".pdf",sep=""),onefile=TRUE,width=10,height=10)
par(mfrow=c(3,3))
for(i in 1:length(projfilelist)){
  GCM = projfilebreakdown$GCM[i]
  scen = projfilebreakdown$scen[i]
  obs = projfilebreakdown$obs[i]
  DS = projfilebreakdown$DS[i]
  
  testsfc1 = list(x=diffs[[1]],y=diffs[[2]],z=diffs_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()


scensin = scens[c(1,3)]
diffsg1_sort = diffsg1[[3]][,,c(1,3)]

pdf(paste("Group1_",varname,"_",difftype,".pdf",sep=""),onefile=TRUE,width=10,height=5)
diffcolorbar = colorramp(diffsg1[[3]],colorchoice="browntogreen",Blimit=BINLIMIT,type="difference")

par(mfrow=c(1,2))

for(i in 1:length(scensin)){
  testsfc1 = list(x=diffsg1[[1]],y=diffsg1[[2]],z=diffsg1_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scensin[i],sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

