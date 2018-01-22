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

###
# Arguments set

datasetname = "I35" # common name of dataset to start with
varin = "tasmax" # common name of variable

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

files= system(paste("ls /data2/3to5/I35/",varin,"/EDQM/*.nc",sep=""),intern=T)

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
timeunits = test$var[[6]]$dim[[3]]$units # actual data is index 6 in test$var for 3^5
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

if((varin=="tasmax" | varin=="tasmin" | varin=="pr")){

yearidx = which(as.numeric(substr(times,1,4))==years[y])

if(varin=="tasmax" | varin=="tasmin"){ 
yeararray[,,y] = apply(vardata[,,yearidx],c(1,2),mean,na.rm=TRUE)
if(datasetname=="dcp"){
if(length(which(is.na(apply(vardata[,,yearidx],3,mean,na.rm=TRUE))==FALSE))<350){
yeararray[,,y] = matrix(NA,nrow=dim(yeararray[,,y])[1],ncol=dim(yeararray[,,y])[2])
}
}
}

if(varin=="pr"){
temp = vardata[,,yearidx]
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
temp = ifelse(temp1>=25.4,1,0)
yeararray[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
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

GCMs = subset(GCMnamelist,is.na(GCMnamelist[,which(names(GCMnamelist)==datasetname)])==FALSE)[,which(names(GCMnamelist)==datasetname)]
yearshist = 1981:2005
yearsfuture = 2006:2099
yearsfull = c(yearshist,yearsfuture)

outputdata2 = list()
fileslist = list()

for(g in 1:length(GCMs)){

GCMidx = grep(GCMs[g],files)

if(datasetname =="maca" & g==1) GCMidx = GCMidx[c(1,5,6)]
if(datasetname =="maca" & g==17) GCMidx = GCMidx[4:6]
if(datasetname =="cmip5.bcca" & g==16) GCMidx = GCMidx[6:10]
if(datasetname =="ccr" & GCMs[g]=="cccma_cgcm3_1") GCMidx = GCMidx[-grep("cccma_cgcm3_1_t63",files[grep(GCMs[g],files)])]

if(datasetname!="regcm"){
if(datasetname!="dcp" & datasetname!="ccr"){
histidx = GCMidx[1]
futureidx = GCMidx[-1]
} else {
if(datasetname=="dcp"){
histidx = GCMidx[grep("hist",files[GCMidx])]
futureidx = GCMidx[-grep("hist",files[GCMidx])]
}
if(datasetname=="ccr"){
histidx = GCMidx[grep("twentyc3m",files[GCMidx])]
futureidx = GCMidx[-grep("twentyc3m",files[GCMidx])]
}
}
}

filesused = files[GCMidx]
yearsused = yearlist[GCMidx]

if(g==1) counter=1

if(datasetname!="regcm"){

for(i in 1:length(futureidx)){

outarray = array(NA,dim=c(length(lon),length(lat),length(yearsfull))) # combine historical and future data
outarray[,,which(yearsfull %in% yearlist[[histidx]]) ] = outputdata[[histidx]]
if(datasetname!="cmip5.bcca"){
outarray[,,which(yearsfull %in% yearlist[[futureidx[i]]]) ] = outputdata[[futureidx[i]]]
} else {
outarray[,,which(yearsfull %in% yearsfuture) ] = outputdata[[futureidx[i]]][,,(which(yearsfull %in% yearsfuture)-56)] # cmip5.bcca has 2100, which the other datasets do not have.
}

outputdata2[[counter]] = outarray

if(datasetname=="maca"){ # makes common file names for later
filesplit = strsplit(files[futureidx[i]],"_",fixed=TRUE)[[1]]
ems = strsplit(filesplit[4],".",fixed=TRUE)[[1]][1]
emsidx = which(scennamelist[,which(names(scennamelist)==datasetname)]==ems)
GCMidx2 = which(GCMnamelist[,which(names(GCMnamelist)==datasetname)]==GCMs[g])
fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_",scennamelist[emsidx,1],"_anom.nc",sep="")
} 

if(datasetname=="cmip5.bcca"){ # makes common file names for later
filesplit = strsplit(files[futureidx[i]],"_",fixed=TRUE)[[1]]
ems = filesplit[7]
emsidx = which(scennamelist[,which(names(scennamelist)==datasetname)]==ems)
GCMidx2 = which(GCMnamelist[,which(names(GCMnamelist)==datasetname)]==GCMs[g])

fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_",scennamelist[emsidx,1],"_anom.nc",sep="")
} 

if(datasetname=="dcp"){ # makes common file names for later
filesplit = strsplit(files[futureidx[i]],"-",fixed=TRUE)[[1]]
ems = filesplit[2]
if(length(filesplit)==6) ems = filesplit[3]
emsidx = which(scennamelist[,which(names(scennamelist)==datasetname)]==ems)
GCMidx2 = which(GCMnamelist[,which(names(GCMnamelist)==datasetname)]==GCMs[g])
fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_",scennamelist[emsidx,1],"_anom.nc",sep="")
} 

if(datasetname=="ccr"){ # makes common file names for later
filesplit = strsplit(files[futureidx[i]],"-",fixed=TRUE)[[1]]
ems = strsplit(filesplit[1],"/",fixed=TRUE)[[1]][4]
emsidx = which(scennamelist[,which(names(scennamelist)==datasetname)]==ems)
GCMidx2 = which(GCMnamelist[,which(names(GCMnamelist)==datasetname)]==GCMs[g])
fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_",scennamelist[emsidx,1],"_anom.nc",sep="")
} 

if(datasetname=="cl10"){ # makes common file names for later
GCMidx2 = which(GCMnamelist[,which(names(GCMnamelist)==datasetname)]==GCMs[g])
#fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_a2_anom.nc",sep="")
fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_rcp85_anom.nc",sep="")
} 


counter=counter+1
}


} else { # regcm data is pre-combined, so I'm just doing file names for that dataset.

outputdata2[[counter]] = outputdata[[GCMidx]]

if(datasetname=="regcm"){ # makes common file names for later
GCMidx2 = which(GCMnamelist[,which(names(GCMnamelist)==datasetname)]==GCMs[g])
fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_a2_anom.nc",sep="")
#fileslist[[counter]] = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_rcp85_anom.nc",sep="")
}

counter=counter+1

}


message("Finished combining files for GCM ",g," / ",length(GCMs))
}

rm(outputdata)

#############
# Calculate historical climatology

outputclimodat = list()
yearsidx = which(yearsfull %in% 1981:2000)

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

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25" | varname=="tmax" | varname=="tmin"){
outputanomdat[[i]][,,t] = outputdata2[[i]][,,t]-outputclimodat[[i]]
}

if(varname=="pr"){
outputanomdat[[i]][,,t] = ((outputdata2[[i]][,,t]-outputclimodat[[i]])/(outputclimodat[[i]]))*100
}

}
message("Finished anomaly calcs for file ",i," / ",length(outputdata2))
}

#testsfc = list(x=lon,y=lat,z=outputanomdat[[1]][,,1])
#surface(testsfc,type="I")
#map("state",add=TRUE)

#############
# Write data with new file and variable names to uncertainty folder

if(varname=="tmax" | varname=="tmin"){
dataunits1 = dataunits2 = "C"
}

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25"){
dataunits1 = "days"
dataunits2 = "%"
}

if(varname=="pr"){
dataunits1 = "mm"
dataunits2 = "%"
}

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
dimT <- dim.def.ncdf( "time","year",yearsfull)

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

for(i in 1:length(outputdata2)){

fileout = paste("uncertainty/SE/anomdat/",fileslist[i],sep="")

var1d <- var.def.ncdf(varname, dataunits1, list(dimX,dimY,dimT), mv )
var2d <- var.def.ncdf(paste(varname,"_climo",sep=""), dataunits1, list(dimX,dimY), mv )
var3d <- var.def.ncdf(paste(varname,"_anom",sep=""), dataunits2, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , list(var1d,var2d,var3d) )
put.var.ncdf( nc1, var1d, outputdata2[[i]]) # no start or count: write all values\
put.var.ncdf( nc1, var2d, outputclimodat[[i]]) # no start or count: write all values\
put.var.ncdf( nc1, var3d, outputanomdat[[i]]) # no start or count: write all values\

close.ncdf(nc1)

message("Finished writing file ",i," / ",length(outputdata2))
}
 
