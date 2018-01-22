######################
#
# Systematic Bias Finder
#
######################

###
# Arguments set

datasetname = "regcm" # common name of dataset to start with
varname = "tmax" # common name of variable

###
# load libraries and master list files
library(ncdf)
library(maps) # these two just to check plotting issues if need be
library(fields)


#######
# Pull in PRISM climatology

climofile = system(paste("ls uncertainty/anomdat/",varname,"_*PRISM*.nc",sep=""),intern=TRUE)

test = open.ncdf(climofile)

lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
obsclimo = get.var.ncdf(test,paste(varname,"_climo",sep=""))

close.ncdf(test)

#######
# Pull in dataset climatology

climomodfiles = system(paste("ls uncertainty/anomdat/",varname,"_*",datasetname,"*.nc",sep=""),intern=TRUE)

modclimo = array(NA,dim=c(length(lon),length(lat),length(climomodfiles)))

for(f in 1:length(climomodfiles)){

test = open.ncdf(climomodfiles[f])

modclimo[,,f] = get.var.ncdf(test,paste(varname,"_climo",sep=""))

close.ncdf(test)

}

#######
# climo diffs

if(varname=="tmax" | varname=="tmin"){
dataunits1 = dataunits2 = "C"
}

if(varname=="pr"){
dataunits1 = "mm"
dataunits2 = "%"
}

for(f in 1:length(climomodfiles)) {

moddiffs = modclimo[,,f] - obsclimo

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )

# Make varables of various dimensionality, for illustration purposes
mv <- -999 # missing value to use

filesplit = strsplit(strsplit(climomodfiles[f],"/",fixed=TRUE)[[1]][3],"_",fixed=TRUE)[[1]]

fileout = paste("uncertainty/biascors/",filesplit[1],"_",filesplit[2],"_",filesplit[3],"_",filesplit[4],"_sysbias.nc",sep="")

var2d <- var.def.ncdf(paste(varname,"_sysbias",sep=""), dataunits1, list(dimX,dimY), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , var2d )
put.var.ncdf( nc1, var2d, moddiffs) # no start or count: write all values\

close.ncdf(nc1)

}








