################
#
# Gridded weight calculator

####
# library and function load

# load libraries

library(ncdf)
library(fields)
library(maps)

# load weight calc function

weight_calc = function(xvals,filelist=projfiles,subs = GCMs){
 
  xmdf = xvals[-length(xvals)]
  xobs = xvals[length(xvals)]
  w = c()
   
  for(i in 1:length(subs)){
    
    idx = grep(subs[i],projfiles)    
    temp = xmdf[idx]
    
    w[i] = abs(1/(round(xobs,5)+abs(round(mean(temp,na.rm=TRUE),5)-round(xobs,5)))) # for inequal weights
     #w[i] = 1/length(subs) # for equal weights	

  }
  
  # Add in correction of the occasional Inf. value in w
  #if(length(which(is.na(w)==TRUE))>0){
  #w[which(is.na(w)==TRUE)]=min(w[-which(is.na(w)==TRUE)],na.rm=TRUE)
  #}

  if(length(which(w==Inf))>0){
  w[which(w==Inf)]=min(w[-which(w==Inf)],na.rm=TRUE)
  }

  W = w/sum(w,na.rm=TRUE)
  
  if(any(is.na(w)==TRUE)) W=rep(NA,length(w))
  
  W  
}

###
# set parameters

varname = "pr25"
weighttype = "GCM"
obsdatasetname = "PRISM"
regfilepath = "regfits"
anfilepath = "analysis"

refyear = 1999 # reference year for weighting calculation
#refperiod = 1981:2000 #reference climo period for testing with the weights
###
# Gather Projection Data - fitted anomalies for reference year only

# determine files
files = system(paste("ls uncertainty/SE/",regfilepath,"/",varname,"*.nc",sep=""),intern=TRUE)

obsidx = grep(obsdatasetname,files)
projfiles = files[-obsidx] # files with projection data
obsfile = files[obsidx] # file with observations

if(varname=="pr" & length(grep("pr25",files))>0){
projfiles = projfiles[-grep("pr25",projfiles)]
#obsfile = obsfile[-grep("pr25",obsfile)]
}

if(varname=="tmax" & length(grep("tmax95",files))>0){
projfiles = projfiles[-grep("tmax95",projfiles)]
obsfile = obsfile[-grep("tmax95",obsfile)]
}

if(varname=="tmin" & length(grep("tmin32",files))>0){
projfiles = projfiles[-grep("tmin32",projfiles)]
obsfile = obsfile[-grep("tmin32",obsfile)]
}

# determine year index for the projdatafiles and reference year

years = 1950:2099
yearidx = which(years==refyear) # reference period
#yearidx = which(years %in% refperiod) # climo period idx

# Gather projection data for reference year

for(f in 1:length(projfiles)){

test = open.ncdf(projfiles[f])

if(f==1){
lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
}

vardata = get.var.ncdf(test,"fittedvalues")
time = get.var.ncdf(test,"time")

close.ncdf(test)

if(f==1) projrefyeararray = array(NA,dim=c(length(lon),length(lat),length(projfiles)))

projrefyeararray[,,f] = vardata[,,yearidx] # single year test
#projrefyeararray[,,f] = apply(vardata[,,yearidx],c(1,2),mean,na.rm=TRUE) # climo test

}



###
# Gather observed data - fitted anomaly for reference year only

test = open.ncdf(obsfile)
obsdata = get.var.ncdf(test,"fittedvalues")
timeobs = get.var.ncdf(test,"time")

close.ncdf(test)

timeobsidx = which(timeobs==refyear) # single year test
#timeobsidx = which(timeobs %in% refperiod) # climo test

obsdatain = obsdata[,,timeobsidx] # single year test
#obsdatain = apply(obsdata[,,timeobsidx],c(1,2),mean,na.rm=TRUE) # climo test

###
# Run Weighting calculation function for GCM Weights

splitnames = do.call(rbind,strsplit(projfiles,"_",fixed=TRUE))

GCMs = unique(splitnames[,2])
DSs = unique(splitnames[,3])

#xmdf = projrefyeararray[18,12,]
#xobs = obsdatain[18,12]

xvals = array(NA,dim=c(dim(projrefyeararray)[1],dim(projrefyeararray)[2],dim(projrefyeararray)[3]+1))

for(j in 1:(dim(projrefyeararray)[3]+1)){
if(j<(dim(projrefyeararray)[3]+1)){
xvals[,,j]=projrefyeararray[,,j]
} else {
xvals[,,j]=obsdatain
}
}

#weight_calc(xvals,filelist=projfiles,subs = DSs)

test = apply(xvals,c(1,2),weight_calc,filelist=projfiles,subs=GCMs) # Calculate GCM weights
GCMweights = aperm(test, c(2, 3, 1))

test =apply(xvals,c(1,2),weight_calc,filelist=projfiles,subs=DSs) # calculate DS weights
DSweights = aperm(test, c(2, 3, 1))

##
# write weights to netcdf file

dataunits = ""

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
dimG <- dim.def.ncdf( "GCMs", "GCMs", 1:length(GCMs) )
dimD <- dim.def.ncdf( "DSs", "DSs", 1:length(DSs) )

# Make varables of various dimensionality, for illustration purposes
mv <- -999 # missing value to use

fileout = paste("uncertainty/SE/",anfilepath,"/",varname,"_GCMweights.nc",sep="") # Write GCM weights file
var1d <- var.def.ncdf("GCMweights", dataunits, list(dimX,dimY,dimG), mv )
nc1 <- create.ncdf(fileout , var1d )
put.var.ncdf( nc1, var1d, GCMweights) # no start or count: write all values\
close.ncdf(nc1)


fileout = paste("uncertainty/SE/",anfilepath,"/",varname,"_DSweights.nc",sep="") # Write GCM weights file
var1d <- var.def.ncdf("DSweights", dataunits, list(dimX,dimY,dimD), mv )
nc1 <- create.ncdf(fileout , var1d )
put.var.ncdf( nc1, var1d, DSweights) # no start or count: write all values\
close.ncdf(nc1)

