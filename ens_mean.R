########
# Step 2 analysis 3^5
# Calculate Ensemble Means

ens_mean = function(step1_filename,histnotes,projnotes,groups="DS"){

# must pass in the following: step1_filename,projnotes,histnotes
#step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2071-2099.nc"
#projnotes = projnotes

  #step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2041-2070_ann.nc"
  #groups = "DS"
  #histnotes = do.call("c",strsplit(histnotes,",",fixed=TRUE))
  #projnotes = do.call("c",strsplit(projnotes,",",fixed=TRUE))

  if(length(grep(",",groups))>0){
    groups = strsplit(groups,",",fixed=TRUE)[[1]]
  }
  
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

histfilebreakdown = do.call(rbind,strsplit(histnotes,",",fixed=TRUE))
histfilebreakdown = do.call(rbind,strsplit(histfilebreakdown,"_",fixed=TRUE))
histfilebreakdown = data.frame(histfilebreakdown)
names(histfilebreakdown) = c("GCM","DS","obs")

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

#####
# Group by Emissions Scenario and alternate groupings

groupcount = length(groups)

if(groupcount==1){
  
  groupnameidx = which(names(projfilebreakdown)==groups[1])
  classes = unique(as.character(projfilebreakdown[,groupnameidx]))
  
  for(g in 1:length(classes)){
    histidx = which(eval(parse(text=paste("histfilebreakdown$",groups[1],sep="")))==classes[g])
    assign(paste("histsg1",classes[g],sep="_"),apply(histlist[,,histidx],c(1,2),mean,na.rm=TRUE))
    scens = unique(projfilebreakdown$scen)
    for(s in 1:length(scens)){
      scenidx = which(projfilebreakdown$scen==scens[s] & eval(parse(text=paste("projfilebreakdown$",groups[1],sep="")))==classes[g])
      assign(paste("projsg1",scens[s],classes[g],sep="_"),apply(projlist[,,scenidx],c(1,2),mean,na.rm=TRUE))
      assign(paste("diffsg1",scens[s],classes[g],sep="_"),apply(diffs[,,scenidx],c(1,2),mean,na.rm=TRUE))
    }
  }
  
}

if(groupcount==2){
  
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

itemlist = "list(var1d,var2d,var3d,var4d,var5d,var6d,var7d"

if(groupcount<1){
  itemlist = paste(itemlist,")",sep="")
}

if(groupcount==1){
  for(g in 1:length(classes)){
    assign(paste("varhd",classes[g],sep="_"),ncvar_def(paste("histmean",classes[g],sep="_"),varunits,longname=paste("Historical Mean - ",classes[g],sep=""), list(dimX,dimY), mv ))
    itemlist = paste(itemlist,",",paste("varhd",classes[g],sep="_"),sep="")
    for(s in 1:length(scens)){
      assign(paste("varpd",scens[s],classes[g],sep="_"),ncvar_def(paste("projmean",scens[s],classes[g],sep="_"),varunits,longname=paste("Projected Mean - ",scens[s],", ",classes[g],sep=""), list(dimX,dimY), mv ))
      assign(paste("vardd",scens[s],classes[g],sep="_"),ncvar_def(paste("projmeandiff",scens[s],classes[g],sep="_"),varunits,longname=paste("Projected Mean Change - ",scens[s],", ",classes[g],sep=""), list(dimX,dimY), mv ))
      itemlist = paste(itemlist,",",paste("varpd",scens[s],classes[g],sep="_"),",",paste("vardd",scens[s],classes[g],sep="_"),sep="")
    }
  }
  itemlist = paste(itemlist,")",sep="")
}

#######
# Create netcdf file

nc <- nc_create(step2_filename,eval(parse(text=itemlist)))

# Write some data to the file
ncvar_put(nc, var1d, histsg1) # no start or count: write all values\
ncvar_put(nc, var2d, projsg1[,,1])
ncvar_put(nc, var3d, projsg1[,,2])
ncvar_put(nc, var4d, projsg1[,,3])
ncvar_put(nc, var5d, diffsg1[,,1])
ncvar_put(nc, var6d, diffsg1[,,2])
ncvar_put(nc, var7d, diffsg1[,,3])

if(groupcount==1){
  for(g in 1:length(classes)){
    ncvar_put(nc,eval(parse(text=paste("varhd",classes[g],sep="_"))),eval(parse(text=paste("histsg1",classes[g],sep="_"))))
    for(s in 1:length(scens)){
      ncvar_put(nc,eval(parse(text=paste("varpd",scens[s],classes[g],sep="_"))),eval(parse(text=paste("projsg1",scens[s],classes[g],sep="_"))))
      ncvar_put(nc,eval(parse(text=paste("vardd",scens[s],classes[g],sep="_"))),eval(parse(text=paste("diffsg1",scens[s],classes[g],sep="_"))))
    }
  }
}

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
                help=paste("Listing of available projection information from each member. ", 
                           "Must have the same number of members as the projected change in the input file.",
                           " Will throw an error if the list is not given. Must take this structure GCM_DS_obs_scen.",
                           " Comma delimited for all members")),
    make_option(c("-d", "--histnotes"), action="store",
                dest='histnotes',
                help=paste("Listing of available historical information from each member. ", 
                           "Must have the same number of members as the historical means in the input file.",
                           " Will throw an error if the list is not given. Must take this structure GCM_DS_obs.",
                           " Comma delimited for all members")),
    make_option(c("-g", "--groups"), action="store",
                dest='groups',default="DS",
                help=paste("Groupings for initial analysis. ", 
                           "Aside from the ensemble mean, pick one to two other groups by which you want ensemble means.",
                           " This can be one of three options - DS, GCM, or obs. This can also be two groups if comma delimited (e.g. DS,obs).",
                           " The default for this is DS alone"))
  )
  
  description = paste('Given the filename from var_calc.R and the metadata list for individual members ', 
                      "this will calculate ensemble mean change and ensemble mean value by emissions scenario, and any additional groupings desired.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --input filename --projnotes projnotes --histnotes histnotes --groups groups")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

ens_mean(step1_filename=parsed.args$filename,projnotes = parsed.args$projnotes,histnotes= parsed.args$histnotes,groups = parsed.args$groups)


