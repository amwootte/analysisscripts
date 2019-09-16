################
#
# Gridded weight calculator

####
# library and function load

# load libraries

library(ncdf4)
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
  
  if(all(is.na(xvals)==FALSE)==TRUE & all(w==Inf)==TRUE){
    message("Applying equal weight in situation where all values match obs")
    w = rep(1/length(subs),length(subs))
  }
  
  W = w/sum(w,na.rm=TRUE)
  
  if(any(is.na(w)==TRUE)) W=rep(NA,length(w))
  
  W  
}

###
# set parameters

varname = "tmax95"
weighttype = "GCM"
obsdatasetname = "METDATA"
regfilepath = "regfits"
anfilepath = "analysis"

refyear = 2005 # reference year for weighting calculation
#refperiod = 1981:2000 #reference climo period for testing with the weights
###
# Gather Projection Data - fitted anomalies for reference year only

# determine files
files = system(paste("ls /home/woot0002/uncertainty/",regfilepath,"/",varname,"_*.nc",sep=""),intern=TRUE)

obsidx = grep(obsdatasetname,files)
projfiles = files[-obsidx] # files with projection data
obsfile = files[obsidx] # file with observations

#if(varname=="pr" & length(grep("pr25",files))>0){
#projfiles = projfiles[-grep("pr25",projfiles)]
#obsfile = obsfile[-grep("pr25",obsfile)]
#}

#if(varname=="tmax" & length(grep("tmax95",files))>0){
#projfiles = projfiles[-grep("tmax95",projfiles)]
#obsfile = obsfile[-grep("tmax95",obsfile)]
#}

#if(varname=="tmin" & length(grep("tmin32",files))>0){
#projfiles = projfiles[-grep("tmin32",projfiles)]
#obsfile = obsfile[-grep("tmin32",obsfile)]
#}

# determine year index for the projdatafiles and reference year

years = 1981:2099
yearidx = which(years==refyear) # reference period
#yearidx = which(years %in% refperiod) # climo period idx

# Gather projection data for reference year

for(f in 1:length(projfiles)){

test = nc_open(projfiles[f])

if(f==1){
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
}

vardata = ncvar_get(test,"fittedvalues")
time = ncvar_get(test,"time")

nc_close(test)

if(f==1) projrefyeararray = array(NA,dim=c(length(lon),length(lat),length(projfiles)))

projrefyeararray[,,f] = vardata[,,yearidx] # single year test
#projrefyeararray[,,f] = apply(vardata[,,yearidx],c(1,2),mean,na.rm=TRUE) # climo test

}



###
# Gather observed data - fitted anomaly for reference year only

test = nc_open(obsfile)
obsdata = ncvar_get(test,"fittedvalues")
timeobs = ncvar_get(test,"time")

nc_close(test)

timeobsidx = which(timeobs==refyear) # single year test
#timeobsidx = which(timeobs %in% refperiod) # climo test

obsdatain = obsdata[,,timeobsidx] # single year test
#obsdatain = apply(obsdata[,,timeobsidx],c(1,2),mean,na.rm=TRUE) # climo test

###
# Run Weighting calculation function for GCM Weights

splitnames = do.call(rbind,strsplit(projfiles,"_",fixed=TRUE))

GCMs = unique(splitnames[,2])
DSs = unique(paste(splitnames[,3],splitnames[,4],sep="_"))

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

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )
dimG <- ncdim_def( "GCMs", "GCMs", 1:length(GCMs) )
dimD <- ncdim_def( "DSs", "DSs", 1:length(DSs) )

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

fileout = paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_GCMweights.nc",sep="") # Write GCM weights file
var1d <- ncvar_def("GCMweights", dataunits, list(dimX,dimY,dimG), mv )
nc1 <- nc_create(fileout , var1d )
ncvar_put( nc1, var1d, GCMweights) # no start or count: write all values\
nc_close(nc1)

fileout = paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_DSweights.nc",sep="") # Write GCM weights file
var1d <- ncvar_def("DSweights", dataunits, list(dimX,dimY,dimD), mv )
nc1 <- nc_create(fileout , var1d )
ncvar_put( nc1, var1d, DSweights) # no start or count: write all values\
nc_close(nc1)

save(list=c("GCMs","DSs"),file=paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_fileorder.Rdata",sep=""))

