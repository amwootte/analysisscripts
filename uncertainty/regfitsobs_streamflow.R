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

datasetname = "08144500" # common name of dataset to start with
varname = "TAVG_12M" # common name of variable

###
# load libraries and master list files
library(ncdf4)
library(maps) # these two just to check plotting issues if need be
library(fields)

###
# determine files to use in calculation

# find files to work with

load(paste("/home/woot0002/uncertainty/anomdat/",varname,"_",datasetname,"_anom.Rdata",sep=""))
year = 1981:2005
#if(varname=="pr"){
#files = files[-grep("pr25",files)]
#}

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
# Start calculating regressions

test1 = glm(outputanomdat ~  I(year^4)+I(year^3)+I(year^2)+year)
fittedvals = as.vector(fitted.values(test1))
resvals = as.vector(residuals(test1))

obstable =  data.frame(year)
obstable$fittedvals = NA
obstable$fittedvals[which(is.na(outputanomdat)==FALSE)] = fittedvals

obstable$resvals = NA
obstable$resvals[which(is.na(outputanomdat)==FALSE)] = resvals

#############
# Write fitted values and 

fileout = paste("/home/woot0002/uncertainty/regfits/",varname,"_",datasetname,"_regfit.Rdata",sep="")
save(list=c("obstable"),file=fileout)




















 
