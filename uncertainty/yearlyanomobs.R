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
#####################

###
# Arguments set

datasetname = "PRISM" # common name of dataset to start with
varname = "pr25" # common name of variable

###
# load libraries and master list files
library(ncdf)
library(maps) # these two just to check plotting issues if need be
library(fields)

###
# determine files to use in calculation

# this is section is helping with the right file path on convection
if(datasetname=="PRISM") inname = "PRISM" 
if(datasetname=="cmip5.bcca") inname="cmip5_bcca"
if(datasetname=="regcm") inname="regcmdata"

if(varname=="tmax" | varname=="tmin" | varname=="pr"){
invar = varname
if(varname=="pr") invar = "ppt"
} else {
if(varname=="tmax95") invar="tmax"
if(varname=="tmin32") invar="tmin"
if(varname=="pr25") invar="ppt"
}

# find files to work with
files = system(paste("ls ",inname,"/res15km/SE/*",invar,"*.nc",sep=""),intern=TRUE)

#if(datasetname=="PRISM") files=files[-grep("climo",files)]


###########
# Start calculating yearly averages for each file

outputdata = list()
yearlist = list()

for(f in 1:length(files)){

test = open.ncdf(files[f])

if(f==1){
lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
}

vardata = get.var.ncdf(test,test$var[[1]]$name)

time = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
close.ncdf(test)

if(datasetname=="PRISM") startdate = as.Date(substr(timeunits,12,21))
times = seq(startdate,by="day",length.out=length(time))
years = unique(as.numeric(substr(times,1,4)))

yearlist[[f]] = years

yeararray = array(NA,dim=c(length(lon),length(lat),length(years)))

for(y in 1:length(years)){

if(varname=="tmax" | varname=="tmin"){
yearidx = which(as.numeric(substr(times,1,4))==years[y])
yeararray[,,y] = apply(vardata[,,yearidx],c(1,2),mean,na.rm=TRUE)
} else {
yearidx = which(as.numeric(substr(times,1,4))==years[y])
if(varname=="tmax95") temp = ifelse(vardata[,,yearidx]>=35,1,0)
if(varname=="tmin32") temp = ifelse(vardata[,,yearidx]<=0,1,0)
#if(varname=="pr25") temp = ifelse(vardata[,,yearidx]/10>=25.4,1,0)
#if(varname=="pr") temp = vardata[,,yearidx]/10
if(varname=="pr25") temp = ifelse(vardata[,,yearidx]>=25.4,1,0)
if(varname=="pr") temp = vardata[,,yearidx]

yeararray[,,y] = apply(temp,c(1,2),sum,na.rm=TRUE)
if(varname=="pr") yeararray[,,y] = ifelse(yeararray[,,y]==0,NA,yeararray[,,y])
yeararray[,,y] = ifelse(is.na(vardata[,,1])==FALSE,yeararray[,,y],NA)
}

}

outputdata[[f]] = yeararray
rm(yeararray)

message("Finished yearly averages for file ",f," / ",length(files))

}

if(varname=="pr") outputdata[[1]] = ifelse(outputdata[[1]]==0,NA,outputdata[[1]]/10)

rm(vardata)

#############
# Calculate historical climatology

outputclimodat = list()
yearsidx = which(yearlist[[1]] %in% 1981:2000)

for(i in 1:length(outputdata)){
outputclimodat[[i]] = apply(outputdata[[i]][,,yearsidx],c(1,2),mean,na.rm=TRUE)
message("Finished climo calcs for file ",i," / ",length(outputdata))
}

#############
# Calculate anomalies

outputanomdat = list()

for(i in 1:length(outputdata)){
outputanomdat[[i]] = outputdata[[i]]
for(t in 1:length(yearlist[[1]])){

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25" | varname=="tmax" | varname=="tmin"){
outputanomdat[[i]][,,t] = outputdata[[i]][,,t]-outputclimodat[[i]]
}

if(varname=="pr"){
outputanomdat[[i]][,,t] = ((outputdata[[i]][,,t]-outputclimodat[[i]])/(outputclimodat[[i]]))*100
}
}
message("Finished anomaly calcs for file ",i," / ",length(outputdata))
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
dimT <- dim.def.ncdf( "time","year",yearlist[[1]])

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

for(i in 1:length(outputdata)){

fileout = paste("uncertainty/SE/anomdat/",varname,"_",datasetname,"_anom.nc",sep="")

var1d <- var.def.ncdf(varname, dataunits1, list(dimX,dimY,dimT), mv )
var2d <- var.def.ncdf(paste(varname,"_climo",sep=""), dataunits1, list(dimX,dimY), mv )
var3d <- var.def.ncdf(paste(varname,"_anom",sep=""), dataunits2, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , list(var1d,var2d,var3d) )
put.var.ncdf( nc1, var1d, outputdata[[i]]) # no start or count: write all values\
put.var.ncdf( nc1, var2d, outputclimodat[[i]]) # no start or count: write all values\
put.var.ncdf( nc1, var3d, outputanomdat[[i]]) # no start or count: write all values\

close.ncdf(nc1)

message("Finished writing file ",i," / ",length(outputdata))
}




















 
