####################
#
# Regression Fit calculator
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
load(paste("/home/woot0002/uncertainty/anomdat/CPREP_",varname,"_",datasetname,"_anom.Rdata",sep=""))

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

  #fit <- glm(yMatrix2 ~  I(time^4)+I(time^3)+I(time^2)+time)
  tempfitvals = yMatrix2
  tempresvals = yMatrix2

   for(i in 1:ncol(yMatrix2)){
 	fit = glm(yMatrix2[,i] ~  I(time^4)+I(time^3)+I(time^2)+time)
	#tempfitvals[as.numeric(names(fitted.values(test1))),i] = fitted.values(test1)   
	#tempresvals[as.numeric(names(fitted.values(test1))),i] = residuals(test1)
  
  if(length(rowswdata)==length(fitted.values(fit))){
  tempfitvals[rowswdata,i] = fitted.values(fit)
  tempresvals[rowswdata,i] = residuals(fit)
  } else {
  tempfitvals[as.numeric(rownames(fitted.values(fit))),i] = fitted.values(fit)
  tempresvals[as.numeric(rownames(fitted.values(fit))),i] = residuals(fit)
  tempfitvals[rowswdata[-which(rowswdata %in% as.numeric(rownames(fitted.values(fit))))],i] = NA
  tempresvals[rowswdata[-which(rowswdata %in% as.numeric(rownames(fitted.values(fit))))],i] = NA
  }

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

fittedvals = outputanomdat
resvals = outputanomdat

for(f in 1:81){

  Fin = f+1

  test1 = glm(outputanomdat[,Fin] ~  I(outputanomdat[,1]^4)+I(outputanomdat[,1]^3)+I(outputanomdat[,1]^2)+outputanomdat[,1])
  
  idxin = which(is.na(outputanomdat[,Fin])==FALSE)
  fittedvals[idxin,Fin] = as.vector(fitted.values(test1))
  resvals[idxin,Fin] = as.vector(residuals(test1))
  
message("Finished regression fits for file ",f," / 81")

}

#############
# Write fitted values and 

fileout = paste("/home/woot0002/uncertainty/regfits/CPREP_",varname,"_",datasetname,"_regfit.Rdata",sep="")
save(list=c("fittedvals","resvals"),file=fileout)


