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

datasetname = "08151500" # common name of dataset to start with
varname = "streamflowmod" # common name of variable

###
# load libraries and master list files
library(ncdf4)
library(maps) # these two just to check plotting issues if need be
library(fields)


###
# determine files to use in calculation

load(paste("/home/woot0002/CPREP_",datasetname,"_2006-2099_fixed.Rdata",sep=""))

###########
# Start calculating yearly averages for ensemble member

histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
projdates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

for(f in 1:length(projresults)){

  histidx = which(histfilebreakdown$DS==projfilebreakdown$DS[f] & histfilebreakdown$GCM==projfilebreakdown$GCM[f] & histfilebreakdown$obs==projfilebreakdown$obs[f])
  
  histtmp = histresults[[histidx]][[1]]
  projtmp = projresults[[f]][[1]]
  
  if(length(histtmp)<length(histdates)){
    hdin = histdates[-which(substr(histdates,6,10)=="02-29")]
  } else {
    hdin = histdates
  }
  
  if(length(projtmp)<length(projdates)){
    pdin = projdates[-which(substr(projdates,6,10)=="02-29")]
  } else {
    pdin = projdates
  }
  
  hyears = as.numeric(substr(hdin,1,4))
  pyears = as.numeric(substr(pdin,1,4))
  
  histframe = data.frame(hyears,histtmp)
  projframe = data.frame(pyears,projtmp)
  
  names(histframe) = c("year",paste(projfilebreakdown$GCM[f],projfilebreakdown$DS[f],projfilebreakdown$obs[f],projfilebreakdown$scen[f],sep="_"))
  names(projframe) = c("year",paste(projfilebreakdown$GCM[f],projfilebreakdown$DS[f],projfilebreakdown$obs[f],projfilebreakdown$scen[f],sep="_"))
  
  tmpdataframe = rbind(histframe,projframe)
  tmpyearly = aggregate(tmpdataframe[,2],by=list(year=tmpdataframe$year),mean,na.rm=TRUE)
  names(tmpyearly)[2] = names(tmpdataframe)[2]
  
  if(f==1){
    outputdata = tmpyearly
  } else {
    outputdata = merge(outputdata,tmpyearly,by="year")
  }
  
message("Finished yearly averages for file ",f," / ",length(projresults))
gc()
}

#############
# Calculate historical climatology

yearsidx = which(outputdata$year %in% 1981:2005)

outputclimodat = apply(outputdata[yearsidx,2:ncol(outputdata)],2,mean,na.rm=TRUE)

#############
# Calculate anomalies

outputanomdat = outputdata

for(i in 2:ncol(outputdata)){

  outputanomdat[,i] = outputdata[,i]-outputclimodat[(i-1)]
  
message("Finished anomaly calcs for file ",i-1," / ",ncol(outputdata)-1)
}

#############
# Write data with new file and variable names to uncertainty folder

fileout = paste("/home/woot0002/uncertainty/anomdat/CPREP_",varname,"_",datasetname,"_anom.Rdata",sep="")

save(list=c("outputanomdat"),file=fileout)


