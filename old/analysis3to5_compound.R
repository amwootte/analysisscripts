###############
# 3^5 Projected Change Analyses - annual only

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("/home/woot0002/analysisfunctions.R")

##############
# User supplied inputs

varname = "heatwaves" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50

difftype="absolute" # type of difference to take, can be either absolute or percent

applymask = NA # if you don't want to apply a state mask, leave this alone

colorchoicediff = "bluetored" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 20 # maximum number of color bins allowed for plotting the projected changes

appfunc = "heatwaves" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
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
  
  if(varname=="gsl"){
    filenames=paste(tmaxhistfilelist[i],tminhistfilelist[i],sep=",")
    varnamesin=c("tasmax","tasmin")
  }
  
  if(varname=="heatwaves"){
    filenames=paste(tmaxhistfilelist[i],tminhistfilelist[i],tmaxq95filelist[i],tminq95filelist[i],sep=",")
    varnamesin=c("tasmax","tasmin","tmaxq95","tminq95")
  }
  
  yearlyoutput=netcdftoyearlycombocalcs(filenames=filenames,varnames=varnamesin,dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=datesin,combofunction=appfunc)
  
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
  
  yearlyoutput=netcdftoyearlycombocalcs(filenames=filenames,varnames=varnamesin,dimnames=c("lon","lat","time"),yearlydataperiod=c(2041,2070),datesdataperiod=datesin,combofunction=appfunc)
  
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

diffcolorbar = diff_colorramp(diffs[[3]],colorchoice=colorchoicediff,Blimit=BINLIMIT)

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
diffcolorbar = diff_colorramp(diffsg1[[3]],colorchoice=colorchoicediff,Blimit=BINLIMIT)

par(mfrow=c(1,2))

for(i in 1:length(scensin)){
  testsfc1 = list(x=diffsg1[[1]],y=diffsg1[[2]],z=diffsg1_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scensin[i],sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

