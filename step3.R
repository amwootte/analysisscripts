########
# Step 3 analysis 3^5
#

step3 = function(step1_filename,projnotes,colorchoicediff,BINLIMIT,diffbartype){

# must pass in the following: step1_filename ,projnotes
# colorchoicediff, BINLIMIT, diffbartype
# step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2071-2099.nc"
# projnotes = projnotes
# colorchoicediff = "bluetored"
# BINLIMIT = 30
# diffbartype = "difference"

split1 = strsplit(step1_filename,"/",fixed=TRUE)[[1]]
split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]

varname = split2[1]
futureperiod = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
difftype = split2[3]
###

test = nc_open(step1_filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
diffs = ncvar_get(test,"projmeandiff")

varunits = test$var[[1]]$units
changeunits = test$var[[3]]$units

nc_close(test)

projfilebreakdown = do.call(rbind,strsplit(projnotes,",",fixed=TRUE))
projfilebreakdown = do.call(rbind,strsplit(projfilebreakdown,"_",fixed=TRUE))
projfilebreakdown = data.frame(projfilebreakdown)
names(projfilebreakdown) = c("GCM","DS","obs","scen")

########
# Plotting Differences

if(all(lon>0)) lon =lon-360

diffcolorbar = colorramp(diffs,colorchoice=colorchoicediff,Blimit=BINLIMIT,type=diffbartype,use_fixed_scale = FALSE)

diffs_sort = diffs[,,order(projfilebreakdown$scen)]
projfilebreakdown = projfilebreakdown[order(projfilebreakdown$scen),]

plotfilename = paste("/data2/3to5/I35/plots/all_mems/IndividualMembers_",varname,"_",futureperiod[1],"-",futureperiod[2],"_",difftype,"change.pdf",sep="")

pdf(plotfilename,onefile=TRUE,width=10,height=10)
par(mfrow=c(3,3))
for(i in 1:length(projfilelist)){
  GCM = projfilebreakdown$GCM[i]
  scen = projfilebreakdown$scen[i]
  obs = projfilebreakdown$obs[i]
  DS = projfilebreakdown$DS[i]
  
  testsfc1 = list(x=lon,y=lat,z=diffs_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()
}

###
# Argument parser

library(optparse)
source("analysisfunctions.R")
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
    make_option(c("-p", "--projnotes"), action="store",
                dest='projnotes',
                help=paste("Listing of available information from each member. ", 
                           "Must have the same number of members as the projected change in the input file.",
                           " Will throw an error if the list is not given. Must take this structure GCM_DS_obs_scen.",
                           " Comma delimited for all members")),
    make_option(c("-c", "--color_choice"), action="store", default="bluetored",
                dest='colorchoicediff',
                help=paste("Listing of available information from each member. ", 
                           "Must have the same number of members as the projected change in the input file.",
                           " Will throw an error if the list is not given. Must take this structure GCM_DS_obs_scen.",
                           " Comma delimited for all members")),
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


