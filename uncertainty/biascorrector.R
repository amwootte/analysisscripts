################
#
# BC correct daily dynamic files

###
# Arguments set

datasetname = "cl10" # common name of dataset to start with
obsdatasetname = "PRISM" # common name of dataset to start with
varname = "tmin" # common name of variable
method ="qqmap"

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

if(datasetname=="cl10"){
files = files[1:3]
}


outputdata = list()
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

if(f==1){
time = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
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

outputdata[[f]] = vardata

}

timemod = time
timeunitsmod = timeunits

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
# Get obs data

test = open.ncdf(files)
vardata = get.var.ncdf(test,test$var[[1]]$name)
timeobs = get.var.ncdf(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
close.ncdf(test)

if(obsdatasetname=="PRISM") startdate = as.Date(substr(timeunits,21,30))
timesobs = seq(startdate,by="day",length.out=length(timeobs))

timeidx = which(as.numeric(substr(timesobs,1,4))>=1981 & as.numeric(substr(timesobs,1,4))<=1999)
vardata = vardata[,,timeidx]
timesobs = timesobs[timeidx]

############
# mod monthly averages

timesmod = times

outputmonavg = list() 

for(f in 1:length(outputdata)){

#years = unique(as.numeric(substr(times,1,4)))
#timeidx = which(as.numeric(substr(times,1,4))>=1981 & as.numeric(substr(times,1,4))<=1999)
#times = times[timeidx]
#years = 1981:1999

yearmon = unique(substr(times,1,7))
montharray = array(NA,dim=c(length(lon),length(lat),length(yearmon)))

for(m in 1:length(yearmon)){
monidx = which(substr(times,1,7)==yearmon[m])
montharray[,,m] = apply(outputdata[[f]][,,monidx],c(1,2),mean,na.rm=TRUE)
}

outputmonavg[[f]] = montharray
rm(montharray)
}

############
# mod monthly climo

outputclimo = list() 

for(f in 1:length(outputdata)){

years = unique(as.numeric(substr(times,1,4)))
timeidx = which(as.numeric(substr(times,1,4))>=1981 & as.numeric(substr(times,1,4))<=1999)

times = times[timeidx]
years = 1981:1999

montharray = array(NA,dim=c(length(lon),length(lat),12))

for(m in 1:12){
monidx = which(as.numeric(substr(times,6,7))==m)
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
monidx = which(as.numeric(substr(times,6,7))==m)
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
outputbc = array(NA,dim=dim(outputdata[[f]]))
for(i in 1:dim(outputbc)[3]){
monidx = as.numeric(substr(timesmod[i],6,7))
dayidx = which(days==substr(timesmod[i],6,10))
outputbc[,,i] = obsdayclimo[,,dayidx]-outputclimo[[f]][,,monidx]+outputdata[[f]][,,i]
}
outputBC[[f]] = outputbc
}
}

# Unbiasing

if(method=="unbiasing"){
outputBC = list()
for(f in 1:length(outputdata)){
outputbc = array(NA,dim=dim(outputdata[[f]]))
for(i in 1:dim(outputbc)[3]){
monidx = as.numeric(substr(timesmod[i],6,7))
yearmonidx = which(yearmon==substr(timesmod[i],1,7))
dayidx = which(days==substr(timesmod[i],6,10))
outputbc[,,i] = outputdata[[f]][,,i]-outputclimo[[f]][,,monidx]+obsclimo[,,monidx]#+obsdayclimo[,,dayidx]outputmonavg[[f]][,,yearmonidx]
}
outputBC[[f]] = outputbc
}
}

# Scaling

if(method=="scaling"){
outputBC = list()
for(f in 1:length(outputdata)){
outputbc = array(NA,dim=dim(outputdata[[f]]))
for(i in 1:dim(outputbc)[3]){
monidx = as.numeric(substr(timesmod[i],6,7))
yearmonidx = which(yearmon==substr(timesmod[i],1,7))
outputbc[,,i] = (outputdata[[f]][,,i]/outputmonavg[[f]][,,yearmonidx])*obsclimo[,,monidx]
}
outputBC[[f]] = outputbc
}
}

# QQmap

if(method=="qqmap"){

timeidx = which(as.numeric(substr(timesmod,1,4))>=1981 & as.numeric(substr(timesmod,1,4))<=1999)
timetemp = timesmod[timeidx]
outputBC = list()
#daysoff = unique(substr(timesmod,6,10))
#days = c()
#days[1:59] = daysoff[1:59]
#days[60] = daysoff[366]
#days[61:366] = daysoff[60:365]

for(f in 1:length(outputdata)){

outputtemp = outputdata[[f]][,,timeidx]
outputbc = array(NA,dim=dim(outputdata[[f]]))

for(m in 1:12){

monidx = which(as.numeric(substr(timesmod,6,7))==m)
monidxobs = which(as.numeric(substr(timesobs,6,7))==m)

for(r in 1:dim(outputbc)[1]){
for(c in 1:dim(outputbc)[2]){

if(any(!is.na(outputtemp[r,c,monidxobs])) & any(!is.na(vardata[r,c,monidxobs]))){
pred = outputtemp[r,c,monidxobs]
ePrd = ecdf(pred)
sim = quantile(vardata[r,c,monidxobs],probs=ePrd(outputdata[[f]][r,c,monidx]),na.rm=TRUE,type=4)
outputbc[r,c,monidx] = sim
}

message("Finished calc for row ",r," and col ",c," and month ",m)

}
}

}

outputBC[[f]] = outputbc
}

}

#testsfc1 = list(x=lon-360,y=lat,z=vardata[,,6759])
#testsfc2 = list(x=lon-360,y=lat,z=outputdata[[f]][,,11142])
#testsfc3 = list(x=lon-360,y=lat,z=outputbc[,,11142])

#par(mfrow=c(1,3))
#surface(testsfc1,main = paste("PRISM ",timesmod[11142],sep=""))
#map("state",add=TRUE)
#surface(testsfc2,main = paste("RAW ",timesmod[11142],sep=""))
#map("state",add=TRUE)
#surface(testsfc3,main = paste("BC ",timesmod[11142],sep=""))
#map("state",add=TRUE)



#####
# Write out new files

filesplit = strsplit(projfiles,"/",fixed=TRUE)

for(f in 1:length(projfiles)){

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
dimT <- dim.def.ncdf( "time",timeunitsmod,timemod)

# Make varables of various dimensionality, for illustration purposes
mv <- -999 # missing value to use
#fileout = paste(filesplit[[f]][1],"/res15km/",method,"BC/",filesplit[[f]][3],sep="")
fileout = paste(filesplit[[f]][1],"/",method,"BC/",filesplit[[f]][3],sep="")

var1d <- var.def.ncdf(varname, varunits, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , var1d )
put.var.ncdf( nc1, var1d, outputBC[[f]]) # no start or count: write all values\

close.ncdf(nc1)

message("Finished writing file ",f," / ",length(outputBC))
}


