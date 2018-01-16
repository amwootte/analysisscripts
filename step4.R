########
# Step 4 analysis 3^5
#

step4 = function(step2_filename,colorchoicediff,BINLIMIT,diffbartype,use_fixed_scale,fixed_scale){

# must pass in the following: step2_filename, colorchoicediff, BINLIMIT, diffbartype
# file format recommended /data2/3to5/I35/ens_means/tasmax_ensmean_absolute_2071-2099.nc
  
split1 = strsplit(step2_filename,"/",fixed=TRUE)[[1]]
split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]

varname = split2[1]
futureperiod = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
difftype = split2[3]
seasonin = substr(split2[5],1,3)

test = nc_open(step2_filename)
rcp26dat = ncvar_get(test,"projmeandiff_rcp26")
rcp45dat = ncvar_get(test,"projmeandiff_rcp45")
rcp85dat = ncvar_get(test,"projmeandiff_rcp85")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
changeunits = test$var[[3]]$units

nc_close(test)

####
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
plotfilename = paste("/data2/3to5/I35/plots/ens_means/EnsMean_",varname,"_",futureperiod[1],"-",futureperiod[2],"_",difftype,"_",seasonin,"_change.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=10,height=5)
diffcolorbar = colorramp(c(rcp26dat,rcp45dat,rcp85dat),colorchoice=colorchoicediff,Blimit=BINLIMIT,use_fixed_scale = use_fixed_scale,fixed_scale=fixed_scale)

par(mfrow=c(1,2))

  testsfc1 = list(x=lon,y=lat,z=rcp26dat)
  surface(testsfc1,type="I",main="Projected Difference from Historical Climate\nScen: RCP2.6",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  #testsfc1 = list(x=lon,y=lat,z=rcp45dat)
  #surface(testsfc1,type="I",main="Projected Difference from Historical Climate\nScen: RCP4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  #map("state",add=TRUE)
  
  testsfc1 = list(x=lon,y=lat,z=rcp85dat)
  surface(testsfc1,type="I",main="Projected Difference from Historical Climate\nScen: RCP8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)

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
                help=paste("Input file is the output file of something created with step2.R. ", 
                           "Must include the path to the file. If file and path not present, an error will be thrown.")),
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
  
  description = paste('Given the filename from step2.R this will produce panel plots  ', 
                      "of ensemble mean projected change by emissions scenario")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --input filename -c color_choice -d difftype -b BINLIMIT -x use_fixed_scale -g fixed_scale")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

step4(step2_filename=parsed.args$filename,colorchoicediff=parsed.args$colorchoicediff,BINLIMIT=parsed.args$BINLIMIT,diffbartype=parsed.args$diffbartype,use_fixed_scale=parsed.args$use_fixed_scale,fixed_scale = parsed.args$fixed_scale)
