##################
#
# Climatology location grab

locgrab_climo = function(filename,loc_name,loc_lon,loc_lat){
  # must pass in the following: loc_lon,loc_lat,loc_name,step1_filename,histnotes,projnotes
  
  split1 = strsplit(filename,"/",fixed=TRUE)[[1]]
  split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]
  
  varname = split2[1]
  
test = nc_open(filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
climodat = ncvar_get(test,varname)

varunits = test$var[[1]]$units

nc_close(test)

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
if(loc_lon>0 & all(lon<0)==TRUE){
if(loc_lon>0) loc_lon=loc_lon-360
}

pointarea = distfunc(loc_lon,loc_lat,modelgrid)
locstart = c(pointarea$R-1,pointarea$C-1)
locend = c(pointarea$R+1,pointarea$C+1)

###
# get point values

climovals = mean(climodat[locstart[1]:locend[1],locstart[2]:locend[2]],na.rm=TRUE)

####
# Create output table

output = data.frame(loc_name,loc_lon,loc_lat,climovals)

####
# Write out file

outfilename = paste("/data2/3to5/I35/point_output/",varname,"_",loc_name,"_histclimo.csv",sep="")
write.table(output,file=outfilename,row.names=FALSE,sep=",")

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
  usage = paste("usage: %prog --input filename -n loc_name -x loc_lon -y loc_lat")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

locgrab_climo(filename=parsed.args$filename,loc_name=parsed.args$loc_name,loc_lon=parsed.args$loc_lon,loc_lat=parsed.args$loc_lat)




