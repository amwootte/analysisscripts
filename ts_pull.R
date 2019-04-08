########
# Step 8 analysis 3^5
# pull time series for location

ts_pull = function(varname,histfilelist,projfilelist,loc_name,loc_lon,loc_lat,useobs=FALSE){

# must pass in the following: histfilelist, projfilelist,varname, varunits,changeunits,loc_name,loc_lon,loc_lat,useobs

# requirements
# varname - short varname, like tasmax, tasmin, pr - should not be used for anything aside from the daily variables. 
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
  #varname = varname
  #loc_name = locationname
  #loc_lon = loc_lon
  #loc_lat = loc_lat
  #histfilelist = histfilelist
  #projfilelist = projfilelist
  
  histfilelist = do.call("c",strsplit(histfilelist,",",fixed=TRUE))
  projfilelist = do.call("c",strsplit(projfilelist,",",fixed=TRUE))
 
  varin = varname
   
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

#############
# 1a- historical calcs

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
years = 1981:2005

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR" & histfilebreakdown$DS[i]!="DeltaSD"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  #### dates to use
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  test = nc_open(histfilelist[i])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  if(i==1){
    histdates = datesin
    ###
    # create model grid
    
    LON = rep(lon,each=length(lat))
    LAT = rep(lat,length(lon))
    R = rep(1:length(lon),each=length(lat))
    C = rep(1:length(lat),length(lon))
    modelgrid = data.frame(R,C,LON,LAT)
    names(modelgrid) = c("R","C","lon","lat")
    if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
    
    ###
    # get cells to use
    loc_lon = as.numeric(loc_lon)
    loc_lat = as.numeric(loc_lat)
    if(loc_lon>0) loc_lon=loc_lon-360
    
    pointarea = distfunc(loc_lon,loc_lat,modelgrid)
    locstart = c(pointarea$R-1,pointarea$C-1)
    locend = c(pointarea$R+1,pointarea$C+1)
  }
  
  tmp=c()
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(datesin,1,4))==years[y])
    test = nc_open(histfilelist[i])
    tempdata = ncvar_get(test,varin,start=c(locstart[1],locstart[2],yearidx[1]),count=c(3,3,length(yearidx)))
    nc_close(test)
    tmp = c(tmp,apply(tempdata,3,mean,na.rm=TRUE))
    rm(tempdata)
    message("Finished time series grab for year ",years[y])
  }
  
  tmpframe=data.frame(datesin,tmp)
  names(tmpframe)[2] = paste(histfilebreakdown$GCM[i],histfilebreakdown$DS[i],histfilebreakdown$obs[i],sep="_")
  filename = paste(loc_name,"_",varname,"_",names(tmpframe)[2],"_",histfilebreakdown$scen[i],".csv",sep="")
  write.table(tmpframe,paste("/data2/3to5/I35/ts_output/",filename,sep=""),sep=",",row.names=FALSE)
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

#############
# 1b- Future calcs

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
years = 2006:2099

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
  
  test = nc_open(projfilelist[i])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  if(i==1){
    projdates = datesin
    ###
    # create model grid
    
    LON = rep(lon,each=length(lat))
    LAT = rep(lat,length(lon))
    R = rep(1:length(lon),each=length(lat))
    C = rep(1:length(lat),length(lon))
    modelgrid = data.frame(R,C,LON,LAT)
    names(modelgrid) = c("R","C","lon","lat")
    if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
    
    ###
    # get cells to use
    loc_lon = as.numeric(loc_lon)
    loc_lat = as.numeric(loc_lat)
    if(loc_lon>0) loc_lon=loc_lon-360
    
    pointarea = distfunc(loc_lon,loc_lat,modelgrid)
    locstart = c(pointarea$R-1,pointarea$C-1)
    locend = c(pointarea$R+1,pointarea$C+1)
  }
  
  tmp=c()
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(datesin,1,4))==years[y])
    test = nc_open(projfilelist[i])
    tempdata = ncvar_get(test,varin,start=c(locstart[1],locstart[2],yearidx[1]),count=c(3,3,length(yearidx)))
    nc_close(test)
    tmp = c(tmp,apply(tempdata,3,mean,na.rm=TRUE))
    rm(tempdata)
    message("Finished time series grab for year ",years[y])
  }
  
  tmpframe=data.frame(datesin,tmp)
  names(tmpframe)[2] = paste(projfilebreakdown$GCM[i],projfilebreakdown$DS[i],projfilebreakdown$obs[i],sep="_")
  filename = paste(loc_name,"_",varname,"_",names(tmpframe)[2],"_",projfilebreakdown$scen[i],".csv",sep="")
  write.table(tmpframe,paste("/data2/3to5/I35/ts_output/",filename,sep=""),sep=",",row.names=FALSE)
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}
#print(range(projlist,na.rm=TRUE))

}

###
# Argument parser

library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
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
    make_option(c("-n", "--loc_name"), action="store",
                dest='loc_name',
                help=paste("Name of the Location you are grabbing information for. ", 
                           "Anything can be used, but it's recommended to use something short and simple. ",
                           "For example, using OKC instead of Oklahoma City. If not provided, the script will throw an error. ")),
    make_option(c("-x", "--loc_lon"), action="store",
                dest='loc_lon',
                help=paste("Longitude of the desired location. ", 
                           "If not provided, the script will throw an error. ")),
    make_option(c("-y", "--loc_lat"), action="store",
                dest='loc_lat',
                help=paste("Latitude of the desired location. If not provided, the script will throw an error. ")),
    make_option(c("-O", "--useobs"), action="store",default=FALSE,
                dest='useobs',
                help=paste("Include the historical observations from the METDATA dataset for the same location. ",
                           "Given METDATA has a finer resolution, a 5 by 5 block around the location is taken compared to a 3x3 block in 3^5."))
  )
  
  description = paste('Given the desired variable name, the list of files (with paths) ', 
                      "and associated options for thresholds, future period, applied function, temporal resolution, and location,  ",
                      "pull a time series for a location of interest ",
                      "from all available members from a collection of downscaled projections.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --var varname --histlist histlist --projlist projlist ", 
                "-n loc_name -x loc_lon -y loc_lat [-O useobs]", 
                " [-h --help]")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

ts_pull(varname=parsed.args$varname,histfilelist=parsed.args$histlist,projfilelist=parsed.args$projlist,
      useobs=parsed.args$useobs,loc_name=parsed.args$loc_name,loc_lon=parsed.args$loc_lon,loc_lat=parsed.args$loc_lat)

