###############
# 3^5 Projected Change Analyses - drought plotter

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##############
# User supplied inputs

varname = "SPEIlow" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname # don't change this
#noleap=TRUE # if your data has leap days change this to FALSE, otherwise leave it alone.  Should be false for 3^5.

if(varname=="SPIlow" | varname=="SPIhigh") varin="SPI" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="SPEIlow" | varname=="SPEIlow") varin="SPEI"

scale = 1 # monthly scale for SPI and SPEI

difftype="absolute" # type of difference to take, can be either absolute or percent

applymask = NA # if you don't want to apply a state mask, leave this alone

colorchoicediff = "greentobrown" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 30 # maximum number of color bins allowed for plotting the projected changes

appfunc = "sum" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

TC = TRUE   # Threshold calculator - should this be calculating a threshold? TRUE (calculate threshold) or FALSE(don't calculate threshold)
TH = 0  # Threshold value - what's the threshold the script should calculate for?
cond = "lte" # Threshold condition - "gte" = greater than or equal to, "lte" = less than or equal to, "gt" = greater than, "lt"= less than
# threshold value and threshold condition are ignored if TC=FALSE
# should be using these to calculate tmax95 and the others, temperature thresholds should be supplied in degrees K, precipitation thresholds in mm
# thresholds for each of the following are
# SPIlow, SPEIlow: -1.5 - SPEI, SPI values less than or equal to -1.5
# SPIhigh, SPEIhigh: 1.5 - SPEI, SPI values greater than or equal to 1.5

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,scale,"/EDQM/",varin,scale,"*historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,scale,"/EDQM/",varin,scale,"*rcp*.nc",sep=""),intern=T)

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

dates = seq(as.Date("1981-01-15"),as.Date("2005-12-15"),by="month")

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  datesin = dates
  
  test=nc_open(histfilelist[i])
  vardata = ncvar_get(test,varin)
  if(varname == "SPEIlow" | varname=="SPIlow") {
    vardata=ifelse(vardata<=TH,1,0)
  }
  if(varname == "SPEIhigh" | varname=="SPIhigh") {
    vardata=ifelse(vardata>=TH,1,0)
  }
  
  if(i==1){
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    histlist = array(NA,dim=c(length(lon),length(lat),length(histfilelist)))
  }
  nc_close(test)
  
  years = unique(as.numeric(substr(dates,1,4)))
  tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
  for(y in 1:length(years)){
    yearidx2 = which(as.numeric(substr(dates,1,4))==years[y])
    tmp[,,y]=apply(vardata[,,yearidx2],c(1,2),sum,na.rm=TRUE)
    tmp[,,y] = ifelse(is.na(vardata[,,1])==FALSE,tmp[,,y],NA)
  }
  
  histlist[,,i] = apply(tmp,c(1,2),mean,na.rm=TRUE)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

########
# Future data grab

dates = seq(as.Date("2006-01-15"),as.Date("2099-12-15"),by="month")
yearsused = 2044:2068
yearidx = which(as.numeric(substr(dates,1,4))>=2044 & as.numeric(substr(dates,1,4))<=2068)

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  #if(noleap==TRUE) dates = dates[-which(substr(dates,6,10)=="02-29")]
  
  datesin = dates
  
  test=nc_open(projfilelist[i])
  vardata = ncvar_get(test,varin,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
  if(varname == "SPEIlow" | varname=="SPIlow") {
    vardata=ifelse(vardata<=TH,1,0)
  }
  if(varname == "SPEIhigh" | varname=="SPIhigh") {
    vardata=ifelse(vardata>=TH,1,0)
  }
  
  if(i==1){
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    projlist = array(NA,dim=c(length(lon),length(lat),length(projfilelist)))
  }
  nc_close(test)
  
  years = unique(as.numeric(substr(dates[yearidx],1,4)))
  tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
  for(y in 1:length(years)){
    yearidx2 = which(as.numeric(substr(dates[yearidx],1,4))==years[y])
    tmp[,,y]=apply(vardata[,,yearidx2],c(1,2),sum,na.rm=TRUE)
    tmp[,,y] = ifelse(is.na(vardata[,,1])==FALSE,tmp[,,y],NA)
  }
  
  projlist[,,i] = apply(tmp,c(1,2),mean,na.rm=TRUE)
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

######
# Difference Calcs

source("analysisfunctions.R")

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

diffcolorbar = colorramp(diffs[[3]],colorchoice=colorchoicediff,Blimit=BINLIMIT,type="difference")

diffs_sort = diffs[[3]][,,order(projfilebreakdown$scen)]
projfilebreakdown = projfilebreakdown[order(projfilebreakdown$scen),]

pdf(paste("IndividualMembers_",varname,scale,"_",difftype,".pdf",sep=""),onefile=TRUE,width=10,height=10)
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

pdf(paste("Group1_",varname,scale,"_",difftype,".pdf",sep=""),onefile=TRUE,width=10,height=5)
diffcolorbar = colorramp(diffsg1[[3]],colorchoice=colorchoicediff,Blimit=BINLIMIT,type="difference")

par(mfrow=c(1,2))

for(i in 1:length(scensin)){
  testsfc1 = list(x=diffsg1[[1]],y=diffsg1[[2]],z=diffsg1_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scensin[i],sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

