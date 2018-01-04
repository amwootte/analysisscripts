###############
# Time Capsule Point grab for OKC

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

difftype="absolute" # type of difference to take, can be either absolute or percent

applymask = NA # if you want to apply a state mask, leave this alone

colorchoicediff = "bluetored" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 25 # maximum number of color bins allowed for plotting the projected changes

appfunc = "mean" # which functions do you apply for yearly calculations? "mean" is used for getting annual average temps for instance. "sum" would be used for annual total rainfall and thresholds
# for precipitation and all the threshold functions this should be "sum", otherwise use "mean"

TC = TRUE   # Threshold calculator - should this be calculating a threshold? TRUE (calculate threshold) or FALSE(don't calculate threshold)
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

distfunc = function(coords = c(lat,lon),dataset){
  
  lat = coords[1]
  lon = coords[2]
  
  dist = 3963 * acos(sin( lat / 57.2958 ) * sin( dataset$LATS / 57.2958 ) + cos( lat / 57.2958 ) * cos( dataset$LATS /57.2958 ) * cos( ( dataset$LONS / 57.2958) - ( lon / 57.2958 ) ))
  
  dataset[which(dist==min(dist)),]
  
}

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
  
  dateidx = which(substr(dates,6,10)=="02-29")
  
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-dateidx]
  
  test = nc_open(histfilelist[i])
  if(i==1){
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    histlist = matrix(NA,nrow=length(datesin),ncol=length(histfilelist))
    
    rows = 1:length(lon) # row number is the index for longitudes
    cols = 1:length(lat) # columns number is the index for latitude
    
    R = rep(rows,each=length(cols)) # This whole section gets lat, lon, rows, and cols, from a vector format to a table with the grids together
    C = rep(cols,length(rows))
    LONS = rep(lon,each=length(cols))
    LATS = rep(lat,length(rows))
    
    gridmeta = data.frame(R,C,LONS,LATS) 
    citylon = -95.8872 # this happens to be coordinates for Oklahoma City, OK
    citylat = 36.1994  # longitude values need to be in degrees W
    
    closest_point = distfunc(coords=c(citylat,citylon),gridmeta) # the result is a table with one row with the indices, lat, and lon for the closest grid point match. Included is the distance in miles.
    r=closest_point$R[1]
    c=closest_point$C[1]
  }
  
  temp = ncvar_get(test,varin,start=c(r,c,1),count=c(1,1,-1))
  nc_close(test)
  
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR") temp = temp[-dateidx]
  
  histlist[,i] = temp
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

######
# projected future

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  
  dateidx = which(substr(dates,6,10)=="02-29")
  
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-dateidx]
  
  test = nc_open(projfilelist[i])
  if(i==1){
    projlist = matrix(NA,nrow=length(datesin),ncol=length(projfilelist))
   }
  
  temp = ncvar_get(test,varin,start=c(r,c,1),count=c(1,1,-1))
  nc_close(test)
  
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR") temp = temp[-dateidx]
  
  projlist[,i] = temp
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

#################
# separate future by scenarios

rcp26_projlist = projlist[,which(as.character(projfilebreakdown$scen)=="rcp26")]
rcp45_projlist = projlist[,which(as.character(projfilebreakdown$scen)=="rcp45")]
rcp85_projlist = projlist[,which(as.character(projfilebreakdown$scen)=="rcp85")]

projfiles_rcp26 = projfilebreakdown[which(as.character(projfilebreakdown$scen)=="rcp26"),]
projfiles_rcp45 = projfilebreakdown[which(as.character(projfilebreakdown$scen)=="rcp45"),]
projfiles_rcp85 = projfilebreakdown[which(as.character(projfilebreakdown$scen)=="rcp85"),]

################
# Turn everything into data frames

dateshist = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
dateshist = dateshist[-which(substr(dateshist,6,10)=="02-29")]

datesfut = datesin
datesfut = datesfut[-which(substr(datesfut,6,10)=="02-29")]


histlist = data.frame(histlist)
names(histlist)= paste(histfilebreakdown$GCM,histfilebreakdown$obs,sep="_")
histlist = cbind(dateshist,histlist)

rcp26_projlist = data.frame(rcp26_projlist)
names(rcp26_projlist)= paste(projfiles_rcp26$GCM,projfiles_rcp26$obs,sep="_")
rcp26_projlist = cbind(datesfut,rcp26_projlist)

rcp45_projlist = data.frame(rcp45_projlist)
names(rcp45_projlist)= paste(projfiles_rcp45$GCM,projfiles_rcp45$obs,sep="_")
rcp45_projlist = cbind(datesfut,rcp45_projlist)

rcp85_projlist = data.frame(rcp85_projlist)
names(rcp85_projlist)= paste(projfiles_rcp85$GCM,projfiles_rcp85$obs,sep="_")
rcp85_projlist = cbind(datesfut,rcp85_projlist)

#######
# save output

save(list=c("varname","histfilebreakdown","projfiles_rcp26","projfiles_rcp45","projfiles_rcp85","histlist","rcp26_projlist","rcp45_projlist","rcp85_projlist"),file=paste(varname,"_tulsa_output.Rdata",sep=""))



