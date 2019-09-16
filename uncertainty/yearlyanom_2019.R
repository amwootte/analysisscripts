####################
#
# Yearly Anomaly calculator
#
# This script will do the following with each dataset
#
# 1) Calculate Yearly averages (or sums) - also doing scale factor corrections if needed
# 2) Calculate 1981-2000 climatological average
# 3) Calculate Yearly anomalies (absolute or percent difference)
# 4) Combine historical and future runs into one file
# 5) Make variable, scenario, GCM, and DS names match to common naming in master lists
# 6) Write new yearly netcdfs containing the yearly values, climatology, and averages
# 
# 2018 variation with the output from the 3^5 project
#####################

source("/data2/3to5/I35/scripts/analysisfunctions.R")

###
# Arguments set

datasetname = "I35" # common name of dataset to start with
varname = "rx5day" # common name of variable

if(varname=="tasmax" | varname=="tasmin" | varname=="pr"){
varin=varname
} else {
  if(varname=="tmax95"){
    varin = "tasmax"
  }
  if(varname=="tmin32"){
    varin = "tasmin"
  }
  if(varname=="pr25" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day"){
    varin = "pr"
  }
}


###
# load libraries and master list files
library(ncdf4)
library(maps) # these two just to check plotting issues if need be
library(fields)

#GCMnamelist = read.table("GCMnamelist_SE.csv",sep=",",header=TRUE,colClasses="character") # GCM namelist
#scennamelist = read.table("scennamelist_SE.csv",sep=",",header=TRUE,colClasses="character") # scenario namelist
#varnamelist = read.table("varnamelist_SE.csv",sep=",",header=TRUE,colClasses="character") # variable namelist

###
# determine files to use in calculation

files= c(system(paste("ls /data2/3to5/I35/",varin,"/EDQM/",varin,"_day*.nc",sep=""),intern=TRUE),
         system(paste("ls /data2/3to5/I35/",varin,"/DeltaSD/",varin,"_day*.nc",sep=""),intern=TRUE))

###########
# Start calculating yearly averages for each file

outputdata = list()
yearlist = list()

for(f in 1:length(files)){

test = nc_open(files[f])

if(f==1){
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
}

vardata = ncvar_get(test,varin)

time = ncvar_get(test,"time")
timeunits = test$var[[varin]]$dim[[3]]$units # actual data is index 6 in test$var for 3^5
nc_close(test)

if(datasetname=="regcm" | datasetname=="I35") startdate = as.Date(substr(timeunits,12,21))

if(datasetname=="cl10"){ 
startdate = as.Date(substr(timeunits,13,22))
times = seq(startdate,by="day",length.out=length(time))
} else {
if(datasetname!="I35"){
  times = startdate + time
} else{
  times = startdate + time
  if(length(grep("historical",files[f]))==1){
    dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
    if(length(times)!=length(dates)) times = dates[-grep("02-29",dates)]
  }
  if(length(grep("rcp",files[f]))==1){
    dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
    if(length(times)!=length(dates)) times = dates[-grep("02-29",dates)]
  }
}
}

years = unique(as.numeric(substr(times,1,4)))

yearlist[[f]] = years

yeararray = array(NA,dim=c(length(lon),length(lat),length(years)))

for(y in 1:length(years)){

if((varname=="tasmax" | varname=="tasmin" | varname=="pr")){

yearidx = which(as.numeric(substr(times,1,4))==years[y])

if(varname=="tasmax" | varname=="tasmin"){ 
  message("Calculating average tasmax or tasmin")
yeararray[,,y] = apply(vardata[,,yearidx],c(1,2),mean,na.rm=TRUE)
if(datasetname=="dcp"){
if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
}
}
}

if(varname=="pr"){
  message("calculating total precipitation")
temp = vardata[,,yearidx]*86400
yeararray[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
yeararray[,,y] = ifelse(is.na(vardata[,,yearidx[2]])==TRUE,NA,yeararray[,,y])
if(datasetname=="dcp"){
if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
}
}
}
#message("basic calcs ran")

} else {

yearidx = which(as.numeric(substr(times,1,4))==years[y])

if(length(yearidx)>0){

message("running for year ",y)

temp1 = vardata[,,yearidx]

if(varname=="tmax95"){
  message("calculating tmax95")
temp1 = temp1-273.15
temp = ifelse(temp1>=35,1,0)
yeararray[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
yeararray[,,y] = ifelse(is.na(temp1[,,2])==FALSE,yeararray[,,y],NA)
if(datasetname=="dcp"){
if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
}
}
} 

if(varname=="tmin32"){
  message("calculating tmin32")
  temp1 = temp1-273.15
temp = ifelse(temp1<=0,1,0)
yeararray[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
yeararray[,,y] = ifelse(is.na(temp1[,,2])==FALSE,yeararray[,,y],NA)
if(datasetname=="dcp"){
if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
}
}
}

if(varname=="pr25"){
  message("Calculating pr25")
temp1 = temp1*86400
temp = ifelse(temp1>=25.4,1,0)
yeararray[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
yeararray[,,y] = ifelse(is.na(temp1[,,2])==FALSE,yeararray[,,y],NA)
if(datasetname=="dcp"){
if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
}
}
}

if(varname=="r1mm"){
  message("calculating r1mm")
  temp1 = temp1*86400
  temp = ifelse(temp1>=1,1,0)
  yeararray[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
  yeararray[,,y] = ifelse(is.na(temp1[,,2])==FALSE,yeararray[,,y],NA)
  if(datasetname=="dcp"){
    if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
      yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
    }
  }
}

if(varname=="rx1day"){
  message("calculating rx1day")
  temp1 = temp1*86400
  yeararray[,,y] = apply(temp1,c(1,2),max,na.rm=TRUE)
  yeararray[,,y] = ifelse(is.na(temp1[,,2])==FALSE,yeararray[,,y],NA)
  if(datasetname=="dcp"){
    if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
      yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
    }
  }
}

if(varname=="rx5day"){
  message("calculating rx5day")
  temp1 = temp1*86400
  yeararray[,,y] = apply(temp1,c(1,2),calcrollsum,size=5)
  yeararray[,,y] = ifelse(is.na(temp1[,,2])==FALSE,yeararray[,,y],NA)
  if(datasetname=="dcp"){
    if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
      yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
    }
  }
}

}
}
message("Finished calculation for year ",years[y])
  }

outputdata[[f]] = yeararray
rm(yeararray)

message("Finished yearly averages for file ",f," / ",length(files))
rm(vardata)
gc()
}



if(datasetname=="cl10") lon = lon-360


#############
# Combine the historical and future files for each GCM / scenario combination

fileinfo = do.call(rbind,strsplit(files,"_",fixed=TRUE))
fileinfo2 = do.call(rbind,strsplit(fileinfo[,3],"-",fixed=TRUE))
fileinfo2 = as.data.frame(fileinfo2)
fileinfo2$obs = substr(fileinfo2[,3],4,4)
fileinfo2[,2] = as.character(fileinfo2[,2])
fileinfo2$scen = fileinfo[,4]
fileinfo2$GCM = NA
fileinfo2$GCM[grep("A1",fileinfo2[,3])] = "CCSM4"
fileinfo2$GCM[grep("A2",fileinfo2[,3])] = "MIROC5"
fileinfo2$GCM[grep("A3",fileinfo2[,3])] = "MPI-ESM-LR"

#GCMs = subset(GCMnamelist,is.na(GCMnamelist[,which(names(GCMnamelist)==datasetname)])==FALSE)[,which(names(GCMnamelist)==datasetname)]
GCMs = c("A1","A2","A3")
yearshist = 1981:2005
yearsfuture = 2006:2099
yearsfull = c(yearshist,yearsfuture)

outputdata2 = list()
fileslist = list()

for(g in 1:length(GCMs)){

GCMidx = grep(GCMs[g],files)

if(datasetname=="I35"){
  histidx = GCMidx[grep("historical",files[GCMidx])]
  futureidx = GCMidx[-grep("historical",files[GCMidx])]
}

filesused = files[GCMidx]
yearsused = yearlist[GCMidx]

if(g==1) counter=1
fllist = c()
for(i in 1:length(futureidx)){
  outarray = array(NA,dim=c(length(lon),length(lat),length(yearsfull))) # combine historical and future data
  splitname = strsplit(files[futureidx[i]],GCMs[g],fixed=TRUE)
  obsfind = fileinfo2$obs[futureidx[i]] 
  DSfind = fileinfo2[futureidx[i],2]
  histfileidx = which(fileinfo2[histidx,4]==obsfind & fileinfo2[histidx,2]==DSfind)
  outarray[,,which(yearsfull %in% yearlist[[histidx[histfileidx]]]) ] = outputdata[[histidx[histfileidx]]]
  outarray[,,which(yearsfull %in% yearlist[[futureidx[i]]]) ] = outputdata[[futureidx[i]]]
  outputdata2[[counter]] = outarray
  fileslist[[counter]] = paste(varname,"_",fileinfo2$GCM[futureidx[i]],"_",fileinfo2[futureidx[i],2],"_",fileinfo2$obs[futureidx[i]],"_",fileinfo2$scen[futureidx[i]],"_anom.nc",sep="")
  counter=counter+1
}
  
message("Finished combining files for GCM ",g," / ",length(GCMs))
}

rm(outputdata)

#############
# Calculate historical climatology

outputclimodat = list()
yearsidx = which(yearsfull %in% 1981:2005)

for(i in 1:length(outputdata2)){
outputclimodat[[i]] = apply(outputdata2[[i]][,,yearsidx],c(1,2),mean,na.rm=TRUE)
message("Finished climo calcs for file ",i," / ",length(outputdata2))
}

#############
# Calculate anomalies

outputanomdat = list()

for(i in 1:length(outputdata2)){
outputanomdat[[i]] = outputdata2[[i]]
for(t in 1:length(yearsfull)){

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25" | varname=="r1mm" | varname=="rx5day" | varname=="rx1day" | varname=="tasmax" | varname=="tasmin"){
outputanomdat[[i]][,,t] = outputdata2[[i]][,,t]-outputclimodat[[i]]
}

if(varname=="pr"){
outputanomdat[[i]][,,t] = ((outputdata2[[i]][,,t]-outputclimodat[[i]])/(outputclimodat[[i]]))*100
}

}
message("Finished anomaly calcs for file ",i," / ",length(outputdata2))
}

#############
# Write data with new file and variable names to uncertainty folder

if(varname=="tasmax" | varname=="tasmin"){
dataunits1 = dataunits2 = "K"
}

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25"){
dataunits1 = "days"
dataunits2 = "days"
}

if(varname=="pr"){
dataunits1 = "mm"
dataunits2 = "%"
}

if(varname=="rx1day" | varname=="rx5day"){
  dataunits1 = "mm"
  dataunits2 = "mm"
}

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )
dimT <- ncdim_def( "time","year",yearsfull)

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

for(i in 1:length(outputdata2)){

fileout = paste("/home/woot0002/uncertainty/anomdat/",fileslist[i],sep="")

var1d <- ncvar_def(varname, dataunits1, list(dimX,dimY,dimT), mv )
var2d <- ncvar_def(paste(varname,"_climo",sep=""), dataunits1, list(dimX,dimY), mv )
var3d <- ncvar_def(paste(varname,"_anom",sep=""), dataunits2, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , list(var1d,var2d,var3d) )
ncvar_put( nc1, var1d, outputdata2[[i]]) # no start or count: write all values\
ncvar_put( nc1, var2d, outputclimodat[[i]]) # no start or count: write all values\
ncvar_put( nc1, var3d, outputanomdat[[i]]) # no start or count: write all values\

nc_close(nc1)

message("Finished writing file ",i," / ",length(outputdata2))
}
 
