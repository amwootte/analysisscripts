#######################
#
# Bias correction script

###
# Arguments set

datasetname = "cl10" # common name of dataset to start with
obsdatasetname = "PRISM" # common name of dataset to start with
varname = "tmin" # common name of variable
method ="qqmap"

###
# load libraries, master list files, and functions
library(ncdf)
library(maps) # these two just to check plotting issues if need be
library(fields)

GCMnamelist = read.table("GCMnamelist.csv",sep=",",header=TRUE,colClasses="character") # GCM namelist
scennamelist = read.table("scennamelist.csv",sep=",",header=TRUE,colClasses="character") # scenario namelist
varnamelist = read.table("varnamelist.csv",sep=",",header=TRUE,colClasses="character") # variable namelist

source("BCfunctions.R")

###
# determine files to use in calculation

# this is section is helping with the right file path on convection
if(datasetname=="ccr") inname = "wicci" 
if(datasetname=="cmip5.bcca") inname="cmip5_bcca"
if(datasetname=="regcm") inname="regcmdata"
if(datasetname=="cl10") inname="CLAREnCE10"
if(datasetname=="dcp") inname="dcp"
if(datasetname=="maca") inname="maca"

# find right variable name for the dataset given the common name
colidx = which(names(varnamelist)==datasetname)

if(varname=="tmax" | varname=="tmin" | varname=="pr"){
rowidx = which(varnamelist[,1]==varname)
} else {
if(varname=="tmax95") rowidx = which(varnamelist[,1]=="tmax")
if(varname=="tmin32") rowidx = which(varnamelist[,1]=="tmin")
if(varname=="pr25") rowidx = which(varnamelist[,1]=="pr")
}

invar = varnamelist[rowidx,colidx]

# find files to work with

files = c("CLAREnCE10/res12km/veg_current_20th_HadCM3_TAMIN.nc","CLAREnCE10/res12km/veg_current_21st_HadCM3_TMIN.nc")

#if(datasetname!="maca"){
#files = system(paste("ls ",inname,"/SE/res15km/*",invar,"*.nc",sep=""),intern=TRUE)
#system(paste("ls ",inname,"/res15km/*",invar,"*.nc >> filelist.txt",sep=""))
#} else {
#files = system(paste("ls ",inname,"/res15km/SE/*",invar,"*.nc",sep=""),intern=TRUE)
#system(paste("ls ",inname,"/res15km/qqmapBC/*",invar,"*.nc >> filelist.txt",sep=""))
#}

#if(datasetname == "cmip5.bcca" | datasetname=="dcp"){
#files = files[-grep("climo",files)]
#}

#####################
# Gather model data

outputdata = timesmod = timeunitsmod = timemod = list()
projfiles = files

for(f in 1:length(projfiles)){

test = open.ncdf(projfiles[f])

if(f==1){
lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
}

vardata = get.var.ncdf(test,test$var[[1]]$name)

if(datasetname=="maca" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=vardata-273.15
if(datasetname=="cl10" & varname=="pr") {
vardata = ifelse(vardata<0,0,vardata*3600)
vardata = ifelse(vardata>=600,NA,vardata)
}

if(datasetname=="cl10" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=ifelse(vardata>1000,NA,vardata)
if(datasetname=="cl10" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=ifelse(vardata== -999,NA,vardata)
if(datasetname=="cl10" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=vardata/10 + 273.15
if(datasetname=="regcm" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=vardata+273.15

if(datasetname=="cmip5.bcca" & (varname=="tmax" | varname=="tmin")) vardata = ifelse(vardata< -50,NA,vardata)
if(datasetname=="cmip5.bcca" & varname=="pr") vardata = vardata/10
if(datasetname=="dcp" & (varname=="tmax" | varname=="tmin")){
vardata=ifelse(vardata== -999,NA,vardata)
vardata=(vardata/100)
}

if(datasetname=="dcp" & varname=="pr"){
vardata=ifelse(vardata== -999,NA,vardata)
vardata=(vardata/10)
vardata = ifelse(vardata<0,NA,vardata)
}

time = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
if(f==1){
missingval = test$var[[1]]$missval
varunits = test$var[[1]]$units
}
close.ncdf(test)

if(datasetname=="regcm") startdate = as.Date(substr(timeunits,12,21))

if(datasetname=="cl10"){ 
startdate = as.Date(substr(timeunits,13,22))
times = seq(startdate,by="day",length.out=length(time))
} else {
times = startdate + time
}

timemod[[f]] = time
timeunitsmod[[f]] = timeunits
timesmod[[f]] = times

outputdata[[f]] = vardata

message("Finished getting data for file ",f)
}

############
# Get observed mean

###
# determine files to use in calculation

# this is section is helping with the right file path on convection
if(obsdatasetname=="PRISM") inname = "PRISM" 

if(varname=="tmax" | varname=="tmin" | varname=="pr"){
invar = varname
if(varname=="pr") invar = "ppt"
} else {
if(varname=="tmax95") invar="tmax"
if(varname=="tmin32") invar="tmin"
if(varname=="pr25") invar="ppt"
}

# find files to work with
files = system(paste("ls ",inname,"/res12km/*",invar,"*.nc",sep=""),intern=TRUE)

#if(obsdatasetname=="PRISM") files=files[-grep("climo",files)]

###########
# Get obs data

test = open.ncdf(files)
vardata = get.var.ncdf(test,test$var[[1]]$name)
timeobs = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
close.ncdf(test)

if(obsdatasetname=="PRISM") startdate = as.Date(substr(timeunits,12,21))
timesobs = seq(startdate,by="day",length.out=length(timeobs))

timeidx = which(as.numeric(substr(timesobs,1,4))>=1981 & as.numeric(substr(timesobs,1,4))<=1999)
vardata = vardata[,,timeidx]
timesobs = timesobs[timeidx]

############
# mod monthly climo

outputdata[[1]] = outputdata[[1]]-273.15
outputdata[[2]] = outputdata[[2]]-273.15 

outputclimo = list() 

for(f in 1:2){

years = unique(as.numeric(substr(timesmod[[f]],1,4)))
timeidx = which(as.numeric(substr(timesmod[[f]],1,4))>=1981 & as.numeric(substr(timesmod[[f]],1,4))<=1999)

timesin = timesmod[[f]][timeidx]
years = 1981:1999

montharray = array(NA,dim=c(length(lon),length(lat),12))

for(m in 1:12){
monidx = which(as.numeric(substr(timesin,6,7))==m)
montharray[,,m] = apply(outputdata[[f]][,,monidx],c(1,2),mean,na.rm=TRUE)
}

outputclimo[[f]] = montharray
rm(montharray)
}

############
# obs monthly climo

timeobsfull = timesobs

obsclimo = list() 

years = unique(as.numeric(substr(timesobs,1,4)))
timeidx = which(as.numeric(substr(timesobs,1,4))>=1981 & as.numeric(substr(timesobs,1,4))<=1999)

timesobs = timesobs[timeidx]
years = 1981:1999

montharray = array(NA,dim=c(length(lon),length(lat),12))

for(m in 1:12){
monidx = which(as.numeric(substr(timesobs,6,7))==m)
montharray[,,m] = apply(vardata[,,monidx],c(1,2),mean,na.rm=TRUE)
}

obsclimo = montharray

############
# obs daily climo

obsdayclimo = list() 

years = unique(as.numeric(substr(timesobs,1,4)))
timeidx = which(as.numeric(substr(timesobs,1,4))>=1981 & as.numeric(substr(timesobs,1,4))<=1999)

timesobs = timesobs[timeidx]
years = 1981:1999

days = unique(substr(timesobs,6,10))
dayarray = array(NA,dim=c(length(lon),length(lat),length(days)))

for(d in 1:length(days)){
monidx = which(substr(timesobs,6,10)==days[d])
dayarray[,,d] = apply(vardata[,,monidx],c(1,2),mean,na.rm=TRUE)
}

obsdayclimo = dayarray


############
# BC methods

# Delta

if(method=="delta"){
outputBC = list()
for(f in 1:length(outputdata)){
c=ifelse(f>3,f-3,f)
outputBC[[f]] = deltaBC(modeldata=outputdata[[f]],modclimo=outputclimo[[c]],obsdailyclimo=obsdayclimo,timesmod[[f]])
}
}


plot=FALSE

if(plot==TRUE){

timecheck1 = which(as.numeric(substr(timesmod[[3]],6,7))>=6 & as.numeric(substr(timesmod[[3]],6,7))<=8 & as.numeric(substr(timesmod[[3]],1,4))<2000)
timecheck2 = which(as.numeric(substr(timesmod[[3]],6,7))>=6 & as.numeric(substr(timesmod[[3]],6,7))<=8 & as.numeric(substr(timesmod[[3]],1,4))>=2000)

testsfc1 = list(x=lon,y=lat,z=apply(outputBC[[3]][,,timecheck1],c(1,2),mean,na.rm=TRUE))
testsfc2 = list(x=lon,y=lat,z=apply(outputBC[[3]][,,timecheck2],c(1,2),mean,na.rm=TRUE))
testsfc3 = list(x=lon,y=lat,z=((apply(outputBC[[3]][,,timecheck2],c(1,2),mean,na.rm=TRUE)-apply(outputBC[[3]][,,timecheck1],c(1,2),mean,na.rm=TRUE))))
testsfc4 = list(x=lon,y=lat,z=((apply(outputdata[[3]][,,timecheck2],c(1,2),mean,na.rm=TRUE)-apply(outputdata[[3]][,,timecheck1],c(1,2),mean,na.rm=TRUE))))

par(mfrow=c(1,4))
surface(testsfc1,zlim=c(27,35))
map("state",add=TRUE)
surface(testsfc2,zlim=c(27,35))
map("state",add=TRUE)
surface(testsfc3,zlim=c(1.5,4))
map("state",add=TRUE)
surface(testsfc4,zlim=c(1.5,4))
map("state",add=TRUE)

}

# Unbiasing

if(method=="unbiasing"){
outputBC = list()
for(f in 1:length(outputdata)){
c=ifelse(f>3,f-3,f)
outputBC[[f]] = unbiasingBC(modeldata=outputdata[[f]],modclimo=outputclimo[[c]],obsclimo,timesmod[[f]])
}
}

# Scaling

if(method=="scaling"){
outputBC = list()
for(f in 1:length(outputdata)){
c=ifelse(f>3,f-3,f)
outputBC[[f]] = scalingBC(modeldata=outputdata[[f]],modclimo=outputclimo[[c]],obsclimo,timesmod[[f]])
if(varname=="pr"){
outputBC[[f]] = ifelse(outputBC[[f]]<0,0,outputBC[[f]])
}
}
}

# QQmap

if(method=="qqmap"){

timeidx = which(as.numeric(substr(timesmod[[1]],1,4))>=1981 & as.numeric(substr(timesmod[[1]],1,4))<=1999)
timesmodhist = timesmod[[1]][timeidx]

outputBC = list()

#for(f in 1:length(outputdata)){
#c=ifelse(f>3,f-3,f)

outputBC[[1]] = qqmapBC(modeldata=outputdata[[1]],modelhist=outputdata[[1]][,,timeidx],timesmodhist,timesmod[[1]],obsdata=vardata)

outputBC[[2]] = qqmapBC(modeldata=outputdata[[2]],modelhist=outputdata[[1]][,,timeidx],timesmodhist,timesmod[[2]],obsdata=vardata)


#}

}



#####
# Write out new files

filesplit = strsplit(projfiles,"/",fixed=TRUE)

for(f in 1:length(projfiles)){

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
dimT <- dim.def.ncdf( "time",timeunitsmod[[f]],timemod[[f]])

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use
#fileout = paste(filesplit[[f]][1],"/res15km/",method,"BC/",filesplit[[f]][3],sep="")
fileout = paste(filesplit[[f]][1],"/res12km/",method,"BC/",filesplit[[f]][3],sep="")

var1d <- var.def.ncdf(varname, varunits, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , var1d )
put.var.ncdf( nc1, var1d, outputBC[[f]]) # no start or count: write all values\

close.ncdf(nc1)

message("Finished writing file ",f," / ",length(outputBC))
}





