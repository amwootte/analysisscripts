##################
#
# Step 7: Location grab and calculation

point_calcs = function(step1_filename,histnotes,projnotes,loc_name,loc_lon,loc_lat){
  # must pass in the following: loc_lon,loc_lat,loc_name,step1_filename,histnotes,projnotes
  print(class(histnotes))
  print(class(projnotes))
  histnotes = do.call("c",strsplit(as.character(histnotes),",",fixed=TRUE))
  projnotes = do.call("c",strsplit(as.character(projnotes),",",fixed=TRUE))
  
  split1 = strsplit(step1_filename,"/",fixed=TRUE)[[1]]
  split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]
  
  varname = split2[1]
  futureperiod = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
  difftype = split2[3]
  seasonin = substr(split2[5],1,3)
  
test = nc_open(step1_filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
histlist = ncvar_get(test,"histmean")
projlist = ncvar_get(test,"projmean")
diffs = ncvar_get(test,"projmeandiff")

varunits = test$var[[1]]$units
changeunits = test$var[[3]]$units

nc_close(test)

histfilebreakdown = do.call(rbind,strsplit(histnotes,",",fixed=TRUE))
histfilebreakdown = do.call(rbind,strsplit(histfilebreakdown,"_",fixed=TRUE))
histfilebreakdown = data.frame(histfilebreakdown)
names(histfilebreakdown) = c("GCM","DS","obs")

projfilebreakdown = do.call(rbind,strsplit(projnotes,",",fixed=TRUE))
projfilebreakdown = do.call(rbind,strsplit(projfilebreakdown,"_",fixed=TRUE))
projfilebreakdown = data.frame(projfilebreakdown)
names(projfilebreakdown) = c("GCM","DS","obs","scen")

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

###
# get point values

histvals = apply(histlist[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
projvals = apply(projlist[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
diffvals = apply(diffs[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)

####
# Create output table

projchangedat = projfilebreakdown
projchangedat$modfut = projvals
projchangedat$projchange = diffvals

projchangedat$modhist = NA

for(i in 1:nrow(projchangedat)){
  GCMin = projfilebreakdown$GCM[i]
  obsin = projfilebreakdown$obs[i]
  DSin = projfilebreakdown$DS[i]
  histidx = which(histfilebreakdown$GCM==GCMin & histfilebreakdown$obs==obsin & histfilebreakdown$DS==DSin)
  projchangedat$modhist[i] = histvals[histidx]
}

#projchangedat = projchangedat[,c(3,5,7:11)]

####
# Write out file

filename = paste("/data2/3to5/I35/point_output/",varname,"_",loc_name,"_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_",seasonin,".csv",sep="")
write.table(projchangedat,file=filename,row.names=FALSE,sep=",")

}

###############
###
# Argument parser

library(optparse)
source("/home/woot0002/scripts/analysisfunctions.R")
source("colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

ParseArgs <- function(arg.list){
  option_list <- list(
    #Input and output options first: input, output, timestamp, overlay
    make_option(c('-i', "--input"), action="store",
                dest='filename',
                help=paste("Input file from the first step with the climatology and change calcs from step1. ", 
                           "Must include the path to the file. If not present, an error will be thrown.")),
    make_option(c("-s", "--histnotes"), action="store",
                dest='histnotes',
                help=paste("Listing of available information from each member. ", 
                           "Must have the same number of members as the number of historical members in the input file.",
                           " Will throw an error if the list is not given. Must take this structure GCM_DS_obs.",
                           " Comma delimited for all members")),
    make_option(c("-p", "--projnotes"), action="store",
                dest='projnotes',
                help=paste("Listing of available information from each member. ", 
                           "Must have the same number of members as the projected change in the input file.",
                           " Will throw an error if the list is not given. Must take this structure GCM_DS_obs_scen.",
                           " Comma delimited for all members")),
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
                help=paste("Latitude of the desired location. If not provided, the script will throw an error. "))
    
  )
  
  description = paste('Given the filename from step1.R, the metadata list for individual members,', 
                      " and location information this script provides the projection information for all members for the desired location.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --input filename --projnotes projnotes --histnotes histnotes -n loc_name -x loc_lon -y loc_lat")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

point_calcs(step1_filename=parsed.args$filename,histnotes=parsed.args$histnotes,projnotes = parsed.args$projnotes,loc_name=parsed.args$loc_name,loc_lon=parsed.args$loc_lon,loc_lat=parsed.args$loc_lat)




