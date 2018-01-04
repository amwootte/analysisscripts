######################
#
# 3^5 analysis and imagery wrapper script

###
# Analysis Steps
# 
# 1) Individual Member Calculation of Historical and Projected Values and Projected Change
#   a) Historical Value Calculation
#   b) Projected Value Calculation
#   c) Projected Change Calculation
#   d) Write out files
# 2) Convert to Ensemble Means by Emissions Scenario
#   a) Ensemble Mean Calculation
#   b) Write out files
# 3) Produce Imagery for Individual Members - inputs needed: individual members calculation
# 4) Produce Imagery for Ensemble Means - inputs needed: ensemble means
# 5) File format conversion - inputs needed: ensemble mean calculation
# 6) Area Range Calculation - inputs needed: individual members calculation
# 7) Location Point Calculation - inputs needed: individual members calculation
#
# Steps 3-7 require inputs which are the result of 1 and 2
#
###
#
# Required Abilities
# 1) Modularity - separate parts can be run individually
# 2) Flexibility - Can take new variables and abilities
# 3) Need to account for new directories w/ new sets of projections
#
###
#
# Arguments needed
# 1) short variable name -- tasmax,tasmin,pr,tmax100,tmax95,tmin32,tmin28,frd,gsl,heatwaves,pr25,pr50,mdrn,rx1day,rx5day,cdd,cwd
# 2) DS technique folders in 3^5
# 3) steps to run -- c(1:7), or c(1,2,5,6,7)
# 4) applied function -- sum, mean, heatwaves, growing_season_length,max,rx5day,lastfreeze,maxdryspell,maxwetspell
# 5) difference type -- absolute,relative, or percent
# 6) Temporal period (either say annual, or say a specific season or month) -- annual, JJA, July
# 7) future time period -- c(2041,2070), c(2071,2099)
# 8) units for historical and projected values -- "degrees_K", "mm"
# 9) units for projected change -- "degrees_K", "%"
# 10) Color choices for plotting
#   a) raw values colorbar
#   b) difference values colorbar
# 11) Limit of color bar bins
# 12) legend type
#   a) observed color bar type - raw
#   b) difference color bar type - difference
#
####

library(ncdf4)
library(maps)
library(fields)
library(sp)

source("analysisfunctions.R")

# set arguments
varname = "tmax95"
DStechs = "EDQM"
steps = c(1:7)
appfunc = "sum"
difftype = "absolute"
tempperiod = "annual"
futureperiod = c(2071,2099)
varunits = "Number_of_days"
changeunits = "Number_of_days"
BINLIMIT=30
colorchoicedata = "yellowtored"
colorchoicediff = "bluetored"
obsbartype = "raw"
diffbartype = "difference"
step1_filename = NA
step2_filename = NA

if(6 %in% steps | 3 %in% steps | 7 %in% steps | 2 %in% steps){
  if(1 %in% steps == FALSE & is.na(step1_filename)==TRUE){
    stop("To Calculate steps 2,3,6, and 7, you must calculate step2 or provide the file name and path to the Ensemble means file",.call=TRUE)
  }
}

if(4 %in% steps | 5 %in% steps){
  if(2 %in% steps == FALSE & is.na(step2_filename)==TRUE){
    stop("To Calculate steps 4 and 5, you must calculate step2 or provide the file name and path to the Ensemble means file",.call=TRUE)
  }
}

varin = varname
if(varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd") varin="pr"

TC = FALSE

if(varname == "heatwaves"){
  TC=FALSE
  usecompound=TRUE
  appfunc = "heatwaves"
}

if(varname == "gsl"){
  TC=FALSE
  usecompound=TRUE
  appfunc = "growing_season_length"
}

if(varname == "tmax95"){
  TC=TRUE
  TH = 308.15
  cond = "gte"
  usecompound=FALSE
  appfunc = "sum"
}

if(varname == "tmax100"){
  TC=TRUE
  TH = 310.98
  cond = "gte"
  usecompound=FALSE
  appfunc = "sum"
}

if(varname == "tmin32"){
  TC=TRUE
  TH = 273.15
  cond = "lte"
  usecompound=FALSE
  appfunc = "sum"
}

if(varname == "tmin28"){
  TC=TRUE
  TH = 270.928
  cond = "lte"
  appfunc = "sum"
  usecompound=FALSE
}

if(varname == "frd"){
  TC=FALSE
  appfunc = "lastfreeze"
  usecompound=FALSE
}

if(varname == "mdrn"){
  TC=TRUE
  TH = 0.254
  cond = "gte"
  appfunc = "sum"
  usecompound=FALSE
}
if(varname == "pr25"){
  TC=TRUE
  TH = 25.4
  cond = "gte"
  appfunc = "sum"
  usecompound=FALSE
}
if(varname == "pr50"){
  TC=TRUE
  TH = 50.8
  cond = "gte"
  appfunc = "sum"
  usecompound=FALSE
}
if(varname == "rx1day"){
  TC=FALSE
  appfunc = "max"
  usecompound=FALSE
}
if(varname == "rx5day"){
  TC=FALSE
  appfunc = "rx5day"
  usecompound=FALSE
}

########
# Find all file names

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),9)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),3)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

histlist = paste(histfilelist,collapse=",")
projlist = paste(projfilelist,collapse=",")

########

if(1 %in% steps){
  # run individual calc
  command = paste("Rscript step1.R -v ",varname," -i ",histlist," -p ",projlist," -a ",appfunc," -d ",difftype," -u ",varunits," -x ",changeunits," -f ",futureperiod[1],",",futureperiod[2]," -T ",TC," -H ",TH," -c ",cond,sep="")
  system(command,intern=TRUE)
  step1_filename = paste("/data2/3to5/I35/all_mems/",varname,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],".nc",sep="")
}

########

if(2 %in% steps){
  # run ensemble mean calcs
  command = paste("Rscript step2.R -i ",step1_filename," -p ",paste(projnotes,collapse=","),sep="")
  system(command,intern=TRUE)
  step2_filename = paste("/data2/3to5/I35/ensmeans/",varname,"_ensmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],".nc",sep="")
}

########

if(3 %in% steps){
  # run individual member imagery
  command = paste("Rscript step3.R",varname,step1_filename,difftype,colorbarobs,colorbardiffs,BINLIMIT,obsbartype,diffbartype,futureperiod,tempperiod,sep=" ")
  system(command,intern=TRUE)
}

########

if(4 %in% steps){
  # run ensemble mean imagery
  
}

########

if(5 %in% steps){
  # run ensemble mean file format conversion
  
}

########

if(6 %in% steps){
  # area range calculation
  
}

########

if(7 %in% steps){
  # single location calculation
  
}


