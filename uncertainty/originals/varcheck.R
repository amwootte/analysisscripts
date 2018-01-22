####################
#
# Variability check vs. PRISM
# 
#####################

###
# Arguments set

datasetname = "cl10" # common name of dataset to start with
varname = "tmax" # common name of variable

###
# load libraries and master list files
library(ncdf)
library(maps) # these two just to check plotting issues if need be
library(fields)

GCMnamelist = read.table("GCMnamelist.csv",sep=",",header=TRUE,colClasses="character") # GCM namelist
scennamelist = read.table("scennamelist.csv",sep=",",header=TRUE,colClasses="character") # scenario namelist
varnamelist = read.table("varnamelist.csv",sep=",",header=TRUE,colClasses="character") # variable namelist

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
if(datasetname!="cl10"){
files = system(paste("ls ",inname,"/res15km/*",invar,"*.nc",sep=""),intern=TRUE)
} else {
files = system(paste("ls ",inname,"/res15kmG/*",invar,"*.nc",sep=""),intern=TRUE)
}

if(datasetname == "cmip5.bcca" | datasetname=="dcp"){
files = files[-grep("climo",files)]
}

if(datasetname == "cl10") files = files[-grep("21st",files)]

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

if(datasetname=="maca" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=vardata-273.15
if(datasetname=="cl10" & varname=="pr") {
vardata = ifelse(vardata>=10E10,NA,vardata)
vardata = ifelse(vardata<0,0,vardata)
}
if(datasetname=="cl10" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=ifelse(vardata>1000,NA,vardata)
if(datasetname=="cl10" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=vardata/10
if(datasetname=="cmip5.bcca" & (varname=="tmax" | varname=="tmin")) vardata = ifelse(vardata< -50,NA,vardata)
if(datasetname=="dcp" & (varname=="tmax" | varname=="tmin")){
vardata=ifelse(vardata== -999,NA,vardata)
vardata=(vardata/100)
}

if(datasetname=="dcp" & varname=="pr"){
vardata = ifelse(vardata<0,NA,vardata)
}

time = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
close.ncdf(test)

if(datasetname=="maca") startdate = as.Date(substr(timeunits,12,21))
if(datasetname=="cmip5.bcca") startdate = as.Date(substr(timeunits,11,20))
if(datasetname=="dcp") startdate = as.Date(substr(timeunits,12,21))
if(datasetname=="ccr") startdate = as.Date(substr(timeunits,12,19))
if(datasetname=="regcm") startdate = as.Date(substr(timeunits,12,21))

if(datasetname=="cl10"){ 
startdate = as.Date(substr(timeunits,13,22))
times = seq(startdate,by="day",length.out=length(time))
} else {
times = startdate + time
}

years = unique(as.numeric(substr(times,1,4)))

vararray = apply(vardata,c(1,2),var,na.rm=TRUE)
yearlist[[f]] = years
outputdata[[f]] = vararray

message("Finished yearly averages for file ",f," / ",length(files))

}

