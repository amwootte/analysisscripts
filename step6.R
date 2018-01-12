##################
#
# Step 6: Area Range calculation

step6 = function(step1_filename,histnotes,projnotes,regionname,regiontype,boxXextent,boxYextent,shapefile,shapedimension,areaname){
  
  # must pass in the following: step1_filename,histnotes,projnotes
  # regionname
  # regiontype (box,shapefile)
  # boxXextent = Xmin,Xmax
  # boxYextend = Ymin,Ymax
  # shapefile - shapefile name
  
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
  
  print(summary(modelgrid$lon))
  
  if(regiontype == "box"){
  ###
  # get cells to use
  Xextent = as.numeric(strsplit(boxXextent,",",fixed=TRUE)[[1]])
  Yextent = as.numeric(strsplit(boxYextent,",",fixed=TRUE)[[1]])
  
  print(Xextent)
  print(Yextent)
  
  if(all(Xextent>0)) Xextent=Xextent-360
  
  # get closest points in the grid which match the box
  point1 = distfunc(Xextent[1],Yextent[1],modelgrid) #Xmin, Ymin
  point2 = distfunc(Xextent[2],Yextent[1],modelgrid) #Xmax, Ymin
  point3 = distfunc(Xextent[2],Yextent[2],modelgrid) #Xmax, Ymax
  point4 = distfunc(Xextent[1],Yextent[2],modelgrid) #Xmin, Ymax
  
  print(point1)
  print(point3)
  
  locstart = c(point1$R[1],point1$C[1])
  locend = c(point3$R[1],point3$C[1])
  
  print(locstart)
  print(locend)
  
  ###
  # get point values
  
  histvals = apply(histlist[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
  projvals = apply(projlist[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
  diffvals = apply(diffs[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
  
  }
  
  #######
  # for shapefile based calculations
  
  if(regiontype=="shape"){
    test = readShapePoly(shapefile) # appropriate filename format
    projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
    
    test.sub <- test[as.character(eval(parse(text=paste("test@data$",shapedimension,sep="")))) %in% areaname, ] # in here are two important details. 1) column name in the data array, 2) item in the data array to subset by
    tmpV = histlist[,,1] # need the 
    
    modrasV = raster(t(tmpV)[length(lat):1,])
    if(all(lon>0)){
      extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
    } else {
      extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
    }
    mod.subV <- crop(modrasV, extent(test.sub))
    mod.subV <- mask(modrasV, test.sub)
    
    testin = matrix(getValues(mod.subV),nrow=length(lon),ncol=length(lat))
    testin = testin[,length(lat):1]
    
    maskuse1 = array(NA,dim=c(length(lon),length(lat),dim(histlist)[3]))
    maskuse2 = maskuse3 = array(NA,dim=c(length(lon),length(lat),dim(projlist)[3]))
    for(i in 1:dim(histlist)[3]) maskuse1[,,i] = ifelse(is.na(testin)==FALSE,histlist[,,i],NA)
    for(i in 1:dim(projlist)[3]){
      maskuse2[,,i] = ifelse(is.na(testin)==FALSE,projlist[,,i],NA) 
      maskuse3[,,i] = ifelse(is.na(testin)==FALSE,diffs[,,i],NA)
    } 
    
    histvals = apply(maskuse1,3,mean,na.rm=TRUE)
    projvals = apply(maskuse2,3,mean,na.rm=TRUE)
    diffvals = apply(maskuse3,3,mean,na.rm=TRUE)
    
    #print(diffs[119,78,])
    #print(testin[119,78])
    #print(maskuse3[119,78,])
    #print(diffvals)
    
  }
  
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
  
  filename = paste("/data2/3to5/I35/area_output/",varname,"_",regionname,"_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_",seasonin,".csv",sep="")
  write.table(projchangedat,file=filename,row.names=FALSE,sep=",")
  
}

###############
###
# Argument parser

library(optparse)
source("analysisfunctions.R")
source("colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maptools)

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
    make_option(c("-n", "--regionname"), action="store",
                dest='regionname',
                help=paste("Name of the region you are grabbing information for. ", 
                           "Anything can be used, but it's recommended to use something short and simple. ",
                           "If not provided, the script will throw an error. ")),
    make_option(c("-t", "--regiontype"), action="store", default="box",
                dest='regiontype',
                help=paste("Type of region you want. There are two possible examples you can use. ", 
                           "The first is box (default) and the other is shape (for shapefiles). ",
                           "If box is used then -x and -y must be provided, if shapefile is used then -f, -d, and -a must be provided. ",
                           "Currently only box is functional. ")),
    make_option(c("-x", "--boxXextent"), action="store",
                dest='boxXextent',
                help=paste("Extent of the box desired in the x direction. Should be longitude in degrees East ", 
                           "and supplied in the follow form, Xmin,Xmax. If not supplied an error will be thrown. ",
                           "An error will also be throw if -t is box and this is not provided or ",
                           "If -t is shapefile and this is provided.")),
    make_option(c("-y", "--boxYextent"), action="store",
                dest='boxYextent',
                help=paste("Extent of the box desired in the y direction. Should be latitude in degrees North ", 
                           "and supplied in the follow form, Ymin,Ymax. ",
                           "An error will be thrown if -t is box and this is not provided or ",
                           "If -t is shapefile and this is provided.")),
    make_option(c("-f", "--shapefile"), action="store",
                dest='shapefile',
                help=paste("Shapefile to define the region of interest. Note, do not use the extension (.shp), it will cause an error. ", 
                           "This must include the path to the shapefile desired for use. ",
                           "An error will be thrown if -t is shapefile and this is not provided or ",
                           "If -t is box and this is provided.")),
    make_option(c("-d", "--shapedimension"), action="store",
                dest='shapedimension',
                help=paste("Name of the column in the .dbf of the shapefile from which you are defining the polygon. ", 
                           "An error will be thrown if -t is shape and this is not provided or ",
                           "If -t is box and this is provided.")),
    make_option(c("-a", "--areaname"), action="store",
                dest='areaname',
                help=paste("The name of the area of interest in the column specified by -d. ", 
                           "For example if the STATE_NAME column is defined in -d, then Oklahoma could be provided to -a as the shape. ",
                           "An error will be thrown if -t is shape and this is not provided or ",
                           "If -t is box and this is provided. In addition, this will fail if both -d and -a are not provided if -t is shape."))
  )
  
  description = paste('Given the filename from step1.R, the metadata list for individual members,', 
                      " and region information this script provides the projection information for all members for the desired region.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --input filename --projnotes projnotes --histnotes histnotes -n regionname -t regiontype -x Xmin,Xmax -y Ymin,Ymax -f shapefile")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

step6(step1_filename=parsed.args$filename,histnotes=parsed.args$histnotes,projnotes = parsed.args$projnotes,regiontype=parsed.args$regiontype,regionname=parsed.args$regionname,boxXextent=parsed.args$boxXextent,boxYextent=parsed.args$boxYextent,shapefile=parsed.args$shapefile,shapedimension=parsed.args$shapedimension,areaname=parsed.args$areaname)
