##################
#
# Step 5: File Format Change 

step5 = function(step2_filename,outfiletype){

# must pass in the following: step2_filename, outfiletype
# configured currently only for annual tempperiod, GeoTIFFs, and the ensemble mean calcs
  split1 = strsplit(step2_filename,"/",fixed=TRUE)[[1]]
  split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]
  
  varname = split2[1]
  futureperiod = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
  difftype = split2[3]
  
test = nc_open(step2_filename)
projdiff_rcp26= ncvar_get(test,"projmeandiff_rcp26")
projdiff_rcp45= ncvar_get(test,"projmeandiff_rcp45")
projdiff_rcp85= ncvar_get(test,"projmeandiff_rcp85")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

if(outfiletype="GTiff"){

dataras = raster(t(projdiff_rcp26[,length(lat):1]))
if(all(lon>0)) {
  extent(dataras) = c(lon[1]-360,lon[length(lon)]-360,lat[1],lat[length(lat)])
} else {
  extent(dataras) = c(lon[1],lon[length(lon)],lat[1],lat[length(lat)])
}
crs(dataras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- writeRaster(dataras, filename=paste("/data2/3to5/I35/GeoTIFFs/",varname,"_rcp26_meanchange_",difftype,"_",futureperiod[1],"-",futureperiod[2],".tif",sep=""), format="GTiff", overwrite=TRUE)

dataras = raster(t(projdiff_rcp45[,length(lat):1]))
if(all(lon>0)) {
  extent(dataras) = c(lon[1]-360,lon[length(lon)]-360,lat[1],lat[length(lat)])
} else {
  extent(dataras) = c(lon[1],lon[length(lon)],lat[1],lat[length(lat)])
}
crs(dataras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- writeRaster(dataras, filename=paste("/data2/3to5/I35/GeoTIFFs/",varname,"_rcp45_meanchange_",difftype,"_",futureperiod[1],"-",futureperiod[2],".tif",sep=""), format="GTiff", overwrite=TRUE)

dataras = raster(t(projdiff_rcp85[,length(lat):1]))
if(all(lon>0)) {
  extent(dataras) = c(lon[1]-360,lon[length(lon)]-360,lat[1],lat[length(lat)])
} else {
  extent(dataras) = c(lon[1],lon[length(lon)],lat[1],lat[length(lat)])
}
crs(dataras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- writeRaster(dataras, filename=paste("/data2/3to5/I35/GeoTIFFs/",varname,"_rcp85_meanchange_",difftype,"_",futureperiod[1],"-",futureperiod[2],".tif",sep=""), format="GTiff", overwrite=TRUE)

}

}

#####
# library loading and option parsing

library(optparse)
source("analysisfunctions.R")
source("colorramp.R")
library(ncdf4)
library(rgdal)
library(raster)
library(rasterVis)
library(maps)
library(maptools)

ParseArgs <- function(arg.list){
  option_list <- list(
    #Input and output options first: input, output, timestamp, overlay
    make_option(c('-i', "--input"), action="store",
                dest='filename',
                help=paste("Input file is the output file of something created with step2.R. ", 
                           "Must include the path to the file. If file and path not present, an error will be thrown.")),
    make_option(c("-o", "--outfiletype"), action="store", default="GTiff",
                dest='colorchoicediff',
                help=paste("Outfile type options.", 
                           "Currently GTiff is the only one functioning and the default for GeoTIFF output."))
  )
  
  description = paste('Given the filename from step2.R this will convert the format of the projected change. ', 
                      "Currently only configured to GeoTiffs")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --input filename -o outfiletype")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

step5(step2_filename=parsed.args$filename,outfiletype=parsed.args$outfiletype)
