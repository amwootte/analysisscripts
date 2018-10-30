########
# Step 1 analysis 3^5
#

step1 = function(varname,histfilelist,projfilelist,appfunc,difftype,varunits,changeunits,futureperiod,TC,TH,cond,seasonin){

# must pass in the following: histfilelist,histfilebreakdown, projfilelist,projfilebreakdown,TC,TH,cond, appfunc, varname, futureperiod,tempperiod,varunits,changeunits,usecompound

# requirements
# varname - short varname, like tasmax, tasmin, pr, tmax95,tmin32, pr25,mdrn
# histlist - comma separated list of model historical files
# projlist - comma separated list of model future files
# appfunc - function to apply, like mean, sum, or any of those in the analysisfunctions.R file
# difftype - type of differencing function, "absolute" or "percent"
# varunits - units of the variables themselves, used for netcdf writing
# changeunits - units of the projected change, used for netcdf writing
# futureperiod - comma delimited start and end year, 1981,2005
# TC - is this a threshold calculation? TRUE or FALSE (default)
# TH - threshold value, only used if TC is TRUE. Default is NA
# cond - condition connected with the threshold. Used if TC is TRUE. gt,gte,lt,lte are options
  
  futureperiod = do.call("c",strsplit(futureperiod,",",fixed=TRUE))
  histfilelist = do.call("c",strsplit(histfilelist,",",fixed=TRUE))
  projfilelist = do.call("c",strsplit(projfilelist,",",fixed=TRUE))
 
  varin = varname
  if(varname=="tmax90" | varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
  if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
  if(varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd") varin="pr"
  
  if(varname=="heatwaves" | varname=="gsl"){
    if(length(grep("tasmax",projfilelist))==length(projfilelist)) varin="tasmax"
    if(length(grep("tasmin",projfilelist))==length(projfilelist)) varin="tasmin"
  }
  
   
  # Creating file breakdown tables
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

# 1a- historical calcs

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  if(appfunc!="heatwaves" & appfunc!="growing_season_length"){
  if(seasonin=="ann"){
    yearlyoutput = netcdftoyearlycalcs(histfilelist[i],varname=varin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appliedfunction=appfunc)
  } else {
    yearlyoutput = netcdftoseasonalcalcs(histfilelist[i],varname=varin,season=seasonin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appliedfunction=appfunc)
  }
  }
  
  if(appfunc=="growing_season_length" | appfunc=="heatwaves"){
    message("Running compound calculation, please be patient")
    yearlyoutput = netcdfcombocalcs(histfilelist[i],varname=varin,dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=datesin,combofunction=appfunc)
  }
  
  #print(range(yearlyoutput[[3]],na.rm=TRUE))
  
  if(i==1){
    lon = yearlyoutput[[1]]
    lat = yearlyoutput[[2]]
    histlist = array(NA,dim=c(length(lon),length(lat),length(histfilelist)))
  }
  
  histlist[,,i] = climocalc(yearlyoutput[[3]],yearlydataperiod=c(1981,2005),climoperiod=c(1981,2005))
  rm(yearlyoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}
#print(range(histlist,na.rm=TRUE))
# 1b- Future calcs

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  message("Starting work on file ",projfilelist[i])
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  }else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  if(appfunc!="heatwaves" & appfunc!="growing_season_length"){
  if(seasonin=="ann"){
    yearlyoutput = netcdftoyearlycalcs(projfilelist[i],varname=varin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=futureperiod,datesdataperiod=datesin,appliedfunction=appfunc)
  } else {
    yearlyoutput = netcdftoseasonalcalcs(projfilelist[i],varname=varin,season=seasonin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=futureperiod,datesdataperiod=datesin,appliedfunction=appfunc)
  }
  }
  
  if(appfunc=="growing_season_length" | appfunc=="heatwaves"){
    yearlyoutput = netcdfcombocalcs(projfilelist[i],varname=varin,dimnames=c("lon","lat","time"),yearlydataperiod=futureperiod,datesdataperiod=datesin,combofunction=appfunc)
  }
  
  #print(range(yearlyoutput[[3]],na.rm=TRUE))
  
  if(i==1){
    projlist = array(NA,dim=c(length(lon),length(lat),length(projfilelist)))
  }
  
  projlist[,,i] = climocalc(yearlyoutput[[3]],yearlydataperiod=futureperiod,climoperiod=futureperiod)
  
  rm(yearlyoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}
#print(range(projlist,na.rm=TRUE))
# 1-c Projected Difference Calculation

diffs = array(NA,dim=dim(projlist))
for(i in 1:length(projfilelist)){
  GCMin = projfilebreakdown$GCM[i]
  obsin = projfilebreakdown$obs[i]
  DSin = projfilebreakdown$DS[i]
  histidx = which(histfilebreakdown$GCM==GCMin & histfilebreakdown$obs==obsin & histfilebreakdown$DS==DSin)
  diffs[,,i]=diffcalc(projlist[,,i],histlist[,,histidx],type=difftype)
}

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

# notes on the order of the ensemble members for the output netcdf file.

# 1-d Write out netcdf files

step1_filename = paste(varname,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_",seasonin,".nc",sep="")

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat)
dimF <- ncdim_def("ens_mem_proj",projnotes,1:nrow(projfilebreakdown))
dimH <- ncdim_def("ens_mem_hist",histnotes,1:nrow(histfilebreakdown))

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

var1d <- ncvar_def("histmean",varunits,longname="Historical Mean All Members", list(dimX,dimY,dimH), mv )

var2d <- ncvar_def("projmean",varunits,longname="Projected Mean All Members", list(dimX,dimY,dimF), mv )
var3d <- ncvar_def("projmeandiff",changeunits,longname="Projected Change all Members", list(dimX,dimY,dimF), mv )

#######
# Create netcdf file

nc <- nc_create(paste("/data2/3to5/I35/all_mems/",step1_filename,sep="") ,  list(var1d,var2d,var3d) )

# Write some data to the file
ncvar_put(nc, var1d, histlist) # no start or count: write all values\
ncvar_put(nc, var2d, projlist)
ncvar_put(nc, var3d, diffs)

# close ncdf
nc_close(nc)

returnlist = list(step1_filename,projnotes,histnotes)
return(returnlist)
}

###
# Argument parser

library(optparse)
source("/home/woot0002/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

ParseArgs <- function(arg.list){
  option_list <- list(
    #Input and output options first: input, output, timestamp, overlay
    make_option(c('-v', "--var"), action="store", default='na',
                dest='varname',
                help=paste("The variable of interest (tasmax, tasmin, pr,tmax95,tmin32,pr25) in the input file.", 
                           "If not present, an error will be thrown.")),
    make_option(c("-i", "--histlist"), action="store", default="",
                dest='histlist',
                help=paste("Comma delimited list of filenames for the model historical input files;", 
                           "an error will be thrown if this argument is not present")),
    make_option(c("-p", "--projlist"), action="store", default="",
                dest='projlist',
                help=paste("Comma delimited list of filenames for the model future input files;", 
                           "an error will be thrown if this argument is not present")),
    make_option(c("-a","--appfunc"), action="store", default='mean',dest='appfunc',
                help=paste("The function to be applied prior to determining climatology. ", 
                           "Includes the following options: mean (default), sum,",
                           " max, rx5day (5 day maximum), lastfreeze (Date of last freeze),",
                           "maxdryspell (longest dry spell), and maxwetspell (longest wet spell).")),
    
    make_option(c("-d","--difftype"), action="store", default='absolute',dest='difftype',
                help=paste("The difference calculation to be applied to determine projected change. ", 
                           "Includes the following options: absolute (default) ",
                           "and percent.")),
    
    make_option(c("-u","--varunits"), action="store", default='',dest='varunits',
                help=paste("The units of the variables, used for writing netcdfs.")),
    
    make_option(c("-x","--changeunits"), action="store", default='',dest='changeunits',
                help=paste("The units of the projected change, used for writing netcdfs. ")),
    
    make_option(c("-f","--futureperiod"), action="store", default='2041,2070',dest='futureperiod',
                help=paste("The future period used for the calculation of climatology and projected change. ", 
                           "This should be two values separated by commas for example ",
                           " year1,year2. Will throw an error if the years don't exist in the data.")),
    
    make_option(c("-T","--thres_calc"), action="store", default=FALSE,dest='TC',
                help=paste("Option to calculate certain variables that are threshold based. ", 
                           "Defaults to FALSE, but must be TRUE to calculate variables like ",
                           "tmax95, tmin32, mdrn, pr25. ",
                           "Should also be paired with appfunc option sum. ")),
    make_option(c("-H","--thres_val"), action="store", default='na',dest='TH',
                help=paste("The value of the threshold to be applied if -T is TRUE. ", 
                           "Defaults to a missing value otherwise. Will throw an error ",
                           "if TC is TRUE and -H is not specified")),
    make_option(c("-c","--condition"), action="store", default='na',dest='cond',
                help=paste("The condition that determines how the threshold will", 
                           "be applied to the input data. May be one of LT, LE,", 
                           "GT, GE")),
    make_option(c("-S","--season"), action="store", default='ann',dest='seasonin',
                help=paste("The season to use for calculations. Defaults to 'ann', ", 
                           "but can also be DJF, MAM, JJA, SON."))
  )
  
  description = paste('Given the desired variable name, the list of files (with paths) ', 
                      "and associated options for thresholds, future period, and applied function ",
                      "determine historical climatology, future climatology, and projected change ",
                      "for all available members from a collection of downscaled projections.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --var varname --histlist histlist --projlist projlist ", 
                "[-a appfunc] [-d difftype] [-u varunits] [-x changeunits] [-f futureperiod]", 
                "[-T calcthres] [-H thresvalue] [-c condition] [-S seasonin]", 
                " [-h --help]")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

step1(varname=parsed.args$varname,histfilelist=parsed.args$histlist,projfilelist=parsed.args$projlist,
      appfunc=parsed.args$appfunc,difftype=parsed.args$difftype,varunits=parsed.args$varunits,
      changeunits=parsed.args$changeunits,futureperiod=parsed.args$futureperiod,TC=parsed.args$TC,
      TH=parsed.args$TH,cond=parsed.args$cond,seasonin=parsed.args$seasonin)

