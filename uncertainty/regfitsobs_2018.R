####################
#
# Regression Fit calculator - observations edition
#
# This script will do the following with each dataset
#
# 1) Pulls in yearly anomalies
# 2) Calculates fourth order polynomial regression to each file
# 3) Writes out fits and residuals in one file
# 
#####################

###
# Arguments set

datasetname = "METDATA" # common name of dataset to start with
varname = "tasmax" # common name of variable

###
# load libraries and master list files
library(ncdf4)
library(maps) # these two just to check plotting issues if need be
library(fields)

###
# determine files to use in calculation

# find files to work with

files = system(paste("ls /home/woot0002/uncertainty/anomdat/",varname,"*",datasetname,"*.nc",sep=""),intern=TRUE)

if(varname=="pr"){
files = files[-grep("pr25",files)]
}

###########
# Line fit array function for later

allNAcheck<- function(vecvals){
all(is.na(vecvals)==TRUE)
}


linefitarray <- function (time, y){
  stopifnot(length(dim(y)) == 3, dim(y)[3] == length(time))
  yMatrix <- matrix(aperm(y, c(3, 1, 2)), dim(y)[3])
  allNAs = apply(yMatrix,2,allNAcheck)
  allNAsidx = which(allNAs==TRUE)

  rowswdata = which(apply(yMatrix,1,allNAcheck)==FALSE)
  
  if(length(rowswdata)<length(time)){
  NAidx = 1:length(time)
  NAidx = NAidx[-rowswdata]
  } else {
  NAidx = NULL
  }
  
  yMatrix2 = yMatrix
  if(length(allNAsidx)>0) yMatrix2[,allNAs] = 0

  #test1 = aperm(array(yMatrix2, dim(y)[c(3, 1, 2)]), c(2, 3, 1))

  #fit <- lm(yMatrix2 ~  I(time^4)+I(time^3)+I(time^2)+time)
  #tempfitvals = yMatrix2
  #tempfitvals[rowswdata,] = fitted.values(fit)
  #tempresvals = yMatrix2
  #tempresvals[rowswdata,] = residuals(fit)

  tempfitvals = yMatrix2
  tempresvals = yMatrix2

   for(i in 1:ncol(yMatrix2)){
 	test1 = glm(yMatrix2[,i] ~  I(time^4)+I(time^3)+I(time^2)+time)
	tempfitvals[as.numeric(names(fitted.values(test1))),i] = fitted.values(test1)   
	tempresvals[as.numeric(names(fitted.values(test1))),i] = residuals(test1)
   }
  
  tempfitvals[,allNAsidx]=NA
  tempresvals[,allNAsidx]=NA

  fitvals = aperm(array(tempfitvals, dim(y)[c(3, 1, 2)]), c(2, 3, 1))
  resvals = aperm(array(tempresvals, dim(y)[c(3, 1, 2)]), c(2, 3, 1))
  
  # NA check 
  if(length(NAidx)>0){
  fitvals[,,NAidx]=NA
  resvals[,,NAidx]=NA
  }

  results = list(fitvals,resvals)
  results 
}

###########
# Start calculating yearly averages for each file

fittedvals = list()
resvals = list()

for(f in 1:length(files)){

test = nc_open(files[f])

if(f==1){
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
}

vardata = ncvar_get(test,paste(varname,"_anom",sep=""))

spunits = test$var[[3]]$units

time = ncvar_get(test,"time")
timeunits = test$var[[1]]$dim[[3]]$units
nc_close(test)

fitresults = linefitarray(time,vardata)

fittedvals[[f]] = fitresults[[1]]
resvals[[f]] = fitresults[[2]]
 
message("Finished regression fits for file ",f," / ",length(files))

}

#############
# Write fitted values and 

if(varname=="tasmax" | varname=="tasmin"){
dataunits1 = "C"
}

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25"){
dataunits1 = "days"
}

if(varname=="pr"){
dataunits1 = "%"
}

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )
dimT <- ncdim_def( "time","year",time)

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

for(f in 1:length(files)){

filepart = strsplit(files[f],"/",fixed=TRUE)[[1]]
filepart2 = substr(filepart[length(filepart)],1,nchar(filepart[length(filepart)])-7)

fileout = paste("/home/woot0002/uncertainty/regfits/",filepart2,"regfit.nc",sep="")

var1d <- ncvar_def("fittedvalues", dataunits1, list(dimX,dimY,dimT), mv )
var2d <- ncvar_def("residuals", dataunits1, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , list(var1d,var2d) )
ncvar_put( nc1, var1d, fittedvals[[f]]) # no start or count: write all values\
ncvar_put( nc1, var2d, resvals[[f]]) # no start or count: write all values\

nc_close(nc1)

message("Finished writing file ",f," / ",length(files))
}




















 
