########
# Step 2 analysis 3^5
#

step2 = function(step1_filename,projnotes){

# must pass in the following: step1_filename,projnotes,histnotes
#step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2071-2099.nc"
#projnotes = projnotes

  projnotes = do.call("c",strsplit(projnotes,",",fixed=TRUE))
  
split1 = strsplit(step1_filename,"/",fixed=TRUE)[[1]]
split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]

varname = split2[1]
futureperiod = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
difftype = split2[3]
seasonin = substr(split2[5],1,3)
###

test = nc_open(step1_filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
histlist = ncvar_get(test,"histmean")
projlist = ncvar_get(test,"projmean")
diffs = ncvar_get(test,"projmeandiff")

varunits = test$var[[1]]$units
changeunits = test$var[[3]]$units

nc_close(test)

projfilebreakdown = do.call(rbind,strsplit(projnotes,",",fixed=TRUE))
projfilebreakdown = do.call(rbind,strsplit(projfilebreakdown,"_",fixed=TRUE))
projfilebreakdown = data.frame(projfilebreakdown)
names(projfilebreakdown) = c("GCM","DS","obs","scen")

#####
# Group by Emissions Scenario

histsg1 = apply(histlist,c(1,2),mean,na.rm=TRUE)

diffsg1 = projsg1 = array(NA,dim=c(length(lon),length(lat),length(unique(projfilebreakdown$scen))))
scens = unique(projfilebreakdown$scen)
for(s in 1:length(scens)){
  scenidx = which(projfilebreakdown$scen==scens[s])
  projsg1[,,s] = apply(projlist[,,scenidx],c(1,2),mean,na.rm=TRUE)
  diffsg1[,,s] = apply(diffs[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

######
# Create Ensemble means netcdf

step2_filename = paste("/data2/3to5/I35/ens_means/",varname,"_ensmean_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_",seasonin,".nc",sep="")

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat)

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

var1d <- ncvar_def("histmean",varunits,longname="Historical Mean", list(dimX,dimY), mv )

var2d <- ncvar_def("projmean_rcp26",varunits,longname="Projected Mean RCP 2.6", list(dimX,dimY), mv )
var3d <- ncvar_def("projmean_rcp45",varunits,longname="Projected Mean RCP 4.5", list(dimX,dimY), mv )
var4d <- ncvar_def("projmean_rcp85",varunits,longname="Projected Mean RCP 8.5", list(dimX,dimY), mv )

var5d <- ncvar_def("projmeandiff_rcp26",changeunits,longname="Mean Projected Change RCP 2.6", list(dimX,dimY), mv )
var6d <- ncvar_def("projmeandiff_rcp45",changeunits,longname="Mean Projected Change RCP 4.5", list(dimX,dimY), mv )
var7d <- ncvar_def("projmeandiff_rcp85",changeunits,longname="Mean Projected Change RCP 8.5", list(dimX,dimY), mv )

#######
# Create netcdf file

nc <- nc_create(step2_filename,list(var1d,var2d,var3d,var4d,var5d,var6d,var7d))

# Write some data to the file
ncvar_put(nc, var1d, histsg1) # no start or count: write all values\

ncvar_put(nc, var2d, projsg1[,,1])
ncvar_put(nc, var3d, projsg1[,,2])
ncvar_put(nc, var4d, projsg1[,,3])

ncvar_put(nc, var5d, diffsg1[,,1])
ncvar_put(nc, var6d, diffsg1[,,2])
ncvar_put(nc, var7d, diffsg1[,,3])

# close ncdf
nc_close(nc)

return(step2_filename)
}

###
# Argument parser

library(optparse)
source("analysisfunctions.R")
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
    make_option(c("-p", "--projnotes"), action="store",
                dest='projnotes',
                help=paste("Listing of available information from each member. ", 
                           "Must have the same number of members as the projected change in the input file.",
                           " Will throw an error if the list is not given. Must take this structure GCM_DS_obs_scen.",
                           " Comma delimited for all members"))
  )
  
  description = paste('Given the filename from step1.R and the metadata list for individual members ', 
                      "this will calculate ensemble mean change and ensemble mean value by emissions scenario.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --input filename --projnotes projnotes")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

step2(step1_filename=parsed.args$filename,projnotes = parsed.args$projnotes)


