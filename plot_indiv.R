########
# Step 3 analysis 3^5
# Plot individual members

plot_indiv = function(step1_filename,projnotes,colorchoicediff,BINLIMIT,diffbartype,use_fixed_scale,fixed_scale){

# must pass in the following: step1_filename ,projnotes
# colorchoicediff, BINLIMIT, diffbartype
# step1_filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2041-2070_ann.nc"
# projnotes = projnotes
# colorchoicediff = "bluetored"
# BINLIMIT = 30
# diffbartype = "difference"
# use_fixed_scale = FALSE

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

if(use_fixed_scale==TRUE){
  print(fixed_scale)
  fixed_scale = do.call("c",strsplit(fixed_scale,",",fixed=TRUE))
  print(fixed_scale)
  if(1 %in% grep("n",fixed_scale)) fixed_scale[1] =  paste("-",substr(fixed_scale[1],2,nchar(fixed_scale[1])),sep="")
  if(2 %in% grep("n",fixed_scale)) fixed_scale[2] =  paste("-",substr(fixed_scale[2],2,nchar(fixed_scale[2])),sep="")
  print(fixed_scale)
  fixed_scale=as.numeric(fixed_scale)
  print(fixed_scale)
}

diffcolorbar = colorramp(diffs,colorchoice=colorchoicediff,Blimit=BINLIMIT,type=diffbartype,use_fixed_scale = use_fixed_scale,fixed_scale=fixed_scale)

diffs_sort = diffs[,,order(projfilebreakdown$scen)]
projfilebreakdown = projfilebreakdown[order(projfilebreakdown$scen),]

plotfilename = paste("/data2/3to5/I35/plots/all_mems/IndividualMembers_",varname,"_",futureperiod[1],"-",futureperiod[2],"_",difftype,"_",seasonin,"_change.pdf",sep="")

pdf(plotfilename,onefile=TRUE,width=10,height=10)
par(mfrow=c(3,3))
for(i in 1:length(projnotes)){
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
                help=paste("Color bar choices for plotting. ", 
                           "This includes the following options: bluetored (default), redtoblue, yellowtored,",
                           " whitetored, browntogreen, greentobrown, whitetogreen, and whitetobrown.")),
    make_option(c("-d", "--diffbartype"), action="store", default="difference",
                dest='diffbartype',
                help=paste("The type of plot being produced, important for getting the colorbar right in the internal functions. ", 
                           "This includes the following options: difference (default), raw, or ratio.")),
    make_option(c("-b", "--bin_limit"), action="store", default=30,
                dest='BINLIMIT',
                help=paste("Maximum number of bins allowed in the colorbar. Defaults to 30.")),
    make_option(c("-x","--use_fixed_scale"), action="store", default=FALSE,dest='use_fixed_scale',
                help=paste("Should there be a fixed scale used?", 
                           "Defaults to FALSE (dynamic scale based on input data). ",
                           "Will throw an error if --use_fixed_scale is TRUE, ",
                           "and no fixed scale is provided.")),
    make_option(c("-g","--fixed_scale"), action="store", default='-100,100',dest='fixed_scale',
                help=paste("Fixed scale for plotting, listed as two comma separated values min and max (-g -100,100). ", 
                           "Will throw an error if --use_fixed_scale is TRUE, ",
                           "and no fixed scale is provided."))
  )
  
  description = paste('Given the filename from step1.R and the metadata list for individual members ', 
                      "this will produce the multi-panel plot of the individual members by emissions scenarios.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --input filename --projnotes projnotes -c color_choice -d difftype -b BINLIMIT  -x use_fixed_scale -g fixed_scale")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

plot_indiv(step1_filename=parsed.args$filename,projnotes = parsed.args$projnotes,colorchoicediff=parsed.args$colorchoicediff,BINLIMIT=parsed.args$BINLIMIT,diffbartype=parsed.args$diffbartype,use_fixed_scale=parsed.args$use_fixed_scale,fixed_scale=parsed.args$fixed_scale)



