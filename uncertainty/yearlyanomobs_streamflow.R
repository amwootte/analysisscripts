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
source("/data2/3to5/I35/scripts/analysisfunctions.R")

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
# load data file

# find files to work with
load(paste("/home/woot0002/streamflowreg_3daymovingmean_",datasetname,".Rdata",sep=""))

#if(datasetname=="PRISM") files=files[-grep("climo",files)]
coluse = grep(varname,names(fulldata2))

###########
# Start calculating yearly averages

if(varname=="TMAX" | varname=="TMIN" | varname=="streamflowmod" | varname == "MAXPRCP1W_2W" | varname == "MAXPRCP1W_12M" | varname == "MAXTMAX1W_2W" | varname == "TMAXAVG_1W" | varname == "TAVG_3M" | varname == "MAXTMAX1W_12M" | varname == "TAVG_12M"){
  outputdata = aggregate(fulldata2[,coluse],by=list(year=fulldata2$YEAR),mean,na.rm=TRUE)
}
if(varname=="PRCP" | varname == "PRCPAVG_3D" | varname == "MDRN_4W" | varname == "HOT90_3M" | varname == "HOT90_12M" | varname == "PRCPAVG_6M" | varname == "PRCPAVG_12M"){
  outputdata = aggregate(fulldata2[,coluse],by=list(year=fulldata2$YEAR),sum,na.rm=TRUE)
}

names(outputdata)[2] = varname

#############
# Calculate historical climatology

outputclimodat = mean(outputdata[,2],na.rm=TRUE)

#############
# Calculate anomalies

outputanomdat = outputdata
outputanomdat = outputanomdat[,2]-outputclimodat


#############
# Write data with new file and variable names to uncertainty folder

fileout = paste("/home/woot0002/uncertainty/anomdat/",varname,"_",datasetname,"_anom.Rdata",sep="")

save(list=c("outputanomdat"),file=fileout)





















 
