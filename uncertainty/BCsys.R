######################
#
# Systematic Bias Finder
#
######################

###
# Arguments set

datasetname = "cl10" # common name of dataset to start with
obsdatasetname = "PRISM" # common name of dataset to start with
varname = "pr" # common name of variable

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
files = system(paste("ls ",inname,"/res15kmG/*20th*",invar,"*.nc",sep=""),intern=TRUE)
}

if(datasetname == "cmip5.bcca" | datasetname=="dcp"){
files = files[-grep("climo",files)]
}

fileslist = files

###########
# Start calculating yearly averages for each file

outputdata = list()

for(f in 1:length(files)){

test = open.ncdf(files[f])

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
if(datasetname=="cl10" & (varname=="tmax" | varname=="tmin" | varname=="tmax95" | varname=="tmin32")) vardata=vardata/10
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
timeidx = which(as.numeric(substr(times,1,4))>=1981 & as.numeric(substr(times,1,4))<=1999)

vardata= vardata[,,timeidx]
times = times[timeidx]
years = 1981:1999

montharray = array(NA,dim=c(length(lon),length(lat),12))

for(m in 1:12){

if(varname=="tmax" | varname=="tmin"){
monidx = which(as.numeric(substr(times,6,7))==m)
montharray[,,m] = apply(vardata[,,monidx],c(1,2),mean,na.rm=TRUE)
} 

if(varname=="pr25" | varname=="pr"){
monidx = which(as.numeric(substr(times,6,7))==m)
montharray[,,m] = apply(vardata[,,monidx],c(1,2),mean,na.rm=TRUE)
}

}

outputdata[[f]] = montharray
rm(montharray)

message("Finished yearly averages for file ",f," / ",length(files))

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
files = system(paste("ls ",inname,"/res15km/*",invar,"*.nc",sep=""),intern=TRUE)

if(obsdatasetname=="PRISM") files=files[-grep("climo",files)]


###########
# Start calculating yearly averages for each file

montharray = array(NA,dim=c(length(lon),length(lat),12))

test = open.ncdf(files)

vardata = get.var.ncdf(test,test$var[[1]]$name)

time = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
close.ncdf(test)

if(obsdatasetname=="PRISM") startdate = as.Date(substr(timeunits,21,30))
times = seq(startdate,by="day",length.out=length(time))

timeidx = which(as.numeric(substr(times,1,4))>=1981 & as.numeric(substr(times,1,4))<=1999)

vardata= vardata[,,timeidx]
times = times[timeidx]
years = 1981:1999

for(m in 1:12){

if(varname=="tmax" | varname=="tmin"){
monidx = which(as.numeric(substr(times,6,7))==m)
montharray[,,m] = apply(vardata[,,monidx],c(1,2),mean,na.rm=TRUE)
} 

if(varname=="pr25" | varname=="pr"){
monidx = which(as.numeric(substr(times,6,7))==m)
montharray[,,m] = apply(vardata[,,monidx],c(1,2),mean,na.rm=TRUE)
}

}

#############
# mean difference

meandifference = lapply(outputdata,"-",montharray)

#testsfc = list(x=lon,y=lat,z=meandifference[[1]][,,6])
#surface(testsfc,type="I")

#############
# Write data with new file and variable names to uncertainty folder

if(varname=="tmax" | varname=="tmin"){
dataunits1 = dataunits2 = "C"
}

if(varname=="pr"){
dataunits1 = "mm"
dataunits2 = "mm"
}

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
dimT <- dim.def.ncdf( "time","month",1:12)

# Make varables of various dimensionality, for illustration purposes
mv <- -999 # missing value to use

colidx = which(names(varnamelist)==datasetname)

for(i in 1:length(outputdata)){

if(datasetname=="regcm") GCMidx2 = which(GCMnamelist[,colidx]==strsplit(strsplit(fileslist[i],"/")[[1]][3],"_")[[1]][2])
if(datasetname=="cl10") GCMidx2 = which(GCMnamelist[,colidx]==strsplit(strsplit(fileslist[i],"/")[[1]][3],"_")[[1]][4])

outfilename = paste(varname,"_",GCMnamelist[GCMidx2,1],"_",datasetname,"_a2_sysbias.nc",sep="")
fileout = paste("uncertainty/biascors/",outfilename,sep="")

var1d <- var.def.ncdf("mod_mean", dataunits1, list(dimX,dimY,dimT), mv )
var2d <- var.def.ncdf("obs_mean", dataunits1, list(dimX,dimY,dimT), mv )
var3d <- var.def.ncdf("diff_mean", dataunits2, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , list(var1d,var2d,var3d) )
put.var.ncdf( nc1, var1d, outputdata[[i]]) # no start or count: write all values\
put.var.ncdf( nc1, var2d, montharray) # no start or count: write all values\
put.var.ncdf( nc1, var3d, meandifference[[i]]) # no start or count: write all values\

close.ncdf(nc1)

message("Finished writing file ",i," / ",length(outputdata))
}
 












