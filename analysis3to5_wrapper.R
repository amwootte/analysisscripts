######################
#
# 3^5 analysis and imagery wrapper script

###
# Analysis Steps
# 
# 1) Individual Member Calculation of Historical and Projected Values and Projected Change (var_calc.R)
#   a) Historical Value Calculation
#   b) Projected Value Calculation
#   c) Projected Change Calculation
#   d) Write out files
# 2) Convert to Ensemble Means by Emissions Scenario (ens_mean.R)
#   a) Ensemble Mean Calculation
#   b) Write out files
# 3) Produce Imagery for Individual Members - inputs needed: individual members calculation (plot_indiv.R)
# 4) Produce Imagery for Ensemble Means - inputs needed: ensemble means (plot_ens.R)
# 5) File format conversion - inputs needed: ensemble mean calculation (format_change.R)
# 6) Area Range Calculation - inputs needed: individual members calculation (area_range.R)
# 7) Location Point Calculation - inputs needed: individual members calculation (point_calcs.R)
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
# 6) seasonin (either say annual, or say a specific season) -- annual, JJA
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

source("/data2/3to5/I35/scripts/analysisfunctions.R")

# set arguments
varname = "heatwaves" # these are all required arguments for step 1
steps = c(1,2,3,4) # others can be set based upon what varname is.
difftype = "absolute"
#tempperiod = "annual"
futureperiod = c(2071,2099)
varunits = "events_per_year"
changeunits = "events_per_year"
BINLIMIT=30
colorchoicediff = "bluetored"
diffbartype = "difference"
seasonin = "ann"
useobs=  FALSE

step1_filename = NA # if these are NA and you are not running step1 or step2, then other options that rely on these will break
step2_filename = NA
  
outfileformat = "GTiff" # file format for step 5

lon = c(-101,-94)+360 # information needed for step 6 if regiontype = "box"
lat = c(33,35)
regiontype = "box"
regionname = "RedRiver"

locationname = "BatonRouge" # step 7 required variables
loc_lon = 268.8129
loc_lat = 30.4515

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

############

if(6 %in% steps | 3 %in% steps | 7 %in% steps | 2 %in% steps){
  if(1 %in% steps == FALSE & is.na(step1_filename)==TRUE){
    stop("To Calculate steps 2,3,6, and 7, you must calculate step1 or provide the file name and path to the Individual members file",.call=TRUE)
  }
}

if(4 %in% steps | 5 %in% steps){
  if(2 %in% steps == FALSE & is.na(step2_filename)==TRUE){
    stop("To Calculate steps 4 and 5, you must calculate step2 or provide the file name and path to the Ensemble means file",.call=TRUE)
  }
}

varin = varname
if(varname=="tmax90" | varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd") varin="pr"

TC = FALSE 
TH = NA
cond=NA

if(varname == "tasmax"){
  TC=FALSE
  appfunc = "mean"
}

if(varname == "tasmin"){
  TC=FALSE
  appfunc = "mean"
}

if(varname == "pr"){
  TC=FALSE
  appfunc = "sum"
}

if(varname == "heatwaves"){
  TC=FALSE
  appfunc = "heatwaves"
  varin = "tasmax"
}

if(varname == "gsl"){
  TC=FALSE
  appfunc = "growing_season_length"
  varin = "tasmax"
}

if(varname == "tmax90"){
  TC=TRUE
  TH = 305.372
  cond = "gte"
  appfunc = "sum"
}


if(varname == "tmax95"){
  TC=TRUE
  TH = 308.15
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmax100"){
  TC=TRUE
  TH = 310.928
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmin32"){
  TC=TRUE
  TH = 273.15
  cond = "lte"
  appfunc = "sum"
}

if(varname == "tmin28"){
  TC=TRUE
  TH = 270.928
  cond = "lte"
  appfunc = "sum"
}

if(varname == "frd"){
  TC=FALSE
  appfunc = "lastfreeze"
}

if(varname == "mdrn"){
  TC=TRUE
  TH = 0.254
  cond = "gte"
  appfunc = "sum"
}

if(varname == "pr25"){
  TC=TRUE
  TH = 25.4
  cond = "gte"
  appfunc = "sum"
}

if(varname == "pr50"){
  TC=TRUE
  TH = 50.8
  cond = "gte"
  appfunc = "sum"
}

if(varname == "rx1day"){
  TC=FALSE
  appfunc = "max"
}

if(varname == "rx5day"){
  TC=FALSE
  appfunc = "rx5day"
}

if(varname == "cdd"){
  TC=FALSE
  appfunc = "maxdryspell"
}

if(varname == "cwd"){
  TC=FALSE
  appfunc = "maxwetspell"
}

########
# Find all file names

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

projnotes = paste(projnotes,collapse=",")
histnotes = paste(histnotes,collapse=",")

histlist = paste(histfilelist,collapse=",")
projlist = paste(projfilelist,collapse=",")

########

if(1 %in% steps){
  # run individual calc
  if(TC==FALSE) command = paste("Rscript /data2/3to5/I35/scripts/var_calc.R -v ",varname," -i ",histlist," -p ",projlist," -a ",appfunc," -d ",difftype," -u ",varunits," -x ",changeunits," -f ",futureperiod[1],",",futureperiod[2]," -S ",seasonin,sep="")
  if(TC==TRUE) command = paste("Rscript /data2/3to5/I35/scripts/var_calc.R -v ",varname," -i ",histlist," -p ",projlist," -a ",appfunc," -d ",difftype," -u ",varunits," -x ",changeunits," -f ",futureperiod[1],",",futureperiod[2]," -T ",TC," -H ",TH," -c ",cond," -S ",seasonin,sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand1.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
  step1_filename = paste("/data2/3to5/I35/all_mems/",varname,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_",seasonin,".nc",sep="")
}

########

if(2 %in% steps){
  # run ensemble mean calcs
  command = paste("Rscript ens_mean.R -i ",step1_filename," -p ",projnotes," -d ",histnotes," -g DS", sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand2.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
  step2_filename = paste("/data2/3to5/I35/ens_means/",varname,"_ensmean_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_",seasonin,".nc",sep="")
}

########

if(3 %in% steps){
  # run individual member imagery
  #step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2041-2070_ann.nc"
  command = paste("Rscript plot_indiv.R -i ",step1_filename," -p ",projnotes," -c ",colorchoicediff, " -d ", diffbartype," -b ",BINLIMIT,sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand3.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
}

########

if(4 %in% steps){
  # run ensemble mean imagery
  command = paste("Rscript plot_ens.R -i ",step2_filename," -c ",colorchoicediff, " -d ", diffbartype," -b ",BINLIMIT,sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand4.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
}

########

if(5 %in% steps){
  # run ensemble mean file format conversion
  command = paste("Rscript format_change.R -i ",step2_filename," -o ",outfileformat,sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand5.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
}

########

if(6 %in% steps){
  # area range calculation
  #step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2041-2070_ann.nc"
  if(regiontype=="box") command = paste("Rscript area_range.R -i ",step1_filename," -s ",histnotes," -p ",projnotes," -t ",regiontype," -n ",regionname," -x ",paste(lon[1]+360,lon[2]+360,sep=",")," -y ",paste(lat[1],lat[2],sep=","),sep="")
  if(regiontype=="shape") command = paste("Rscript area_range.R -i ",step1_filename," -s ",histnotes," -p ",projnotes," -t ",regiontype," -n ",regionname," -f ",shapefile," -d ",shapedimension," -a ",areaname,sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand6.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
}

########

if(7 %in% steps){
  # single location calculation
  #step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2041-2070_ann.nc"
  command = paste("Rscript point_calcs.R -i ",step1_filename," -s ",histnotes," -p ",projnotes," -n ",locationname," -x ",loc_lon," -y ",loc_lat,sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand7.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
}

########

if(8 %in% steps){
  # run individual calc
  if(TC==FALSE) command = paste("Rscript ts_pull.R -v ",varname," -i ",histlist," -p ",projlist," -a ",appfunc," -S ",seasonin," -O ",useobs," -n ",locationname," -x ",loc_lon," -y ",loc_lat,sep="")
  if(TC==TRUE) command = paste("Rscript ts_pull.R -v ",varname," -i ",histlist," -p ",projlist," -a ",appfunc," -T ",TC," -H ",TH," -c ",cond," -S ",seasonin," -O ",useobs," -n ",locationname," -x ",loc_lon," -y ",loc_lat,sep="")
  #write.table(command,"/data2/3to5/I35/scripts/testcommand8.txt",sep=",",row.names=FALSE)
  system(command,intern=TRUE)
}


