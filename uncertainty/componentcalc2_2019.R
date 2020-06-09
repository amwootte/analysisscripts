################
#
# Gridded component calculator

############
# library and function load

###########
# Just in case random things got saved, clean the workspace first

rm(list=ls(all=TRUE))

############
# load libraries

library(ncdf4)
library(fields)
library(maps)

############
# load calculator functions

source("/home/woot0002/scripts/uncertainty/2018/componentfunctions_2018.R")

###
# set variable name

varname = "rx1day"
regfilepath = "regfits"
anfilepath = "analysis"

############
# Load fitted values and residuals from all GCMs

# determine files
files = system(paste("ls /home/woot0002/uncertainty/",regfilepath,"/",varname,"_*.nc",sep=""),intern=TRUE)
obsidx = grep("METDATA",files)
files = files[-obsidx] # files with projection data

#if(varname=="pr" & length(grep("pr25",files))>0){
#files = files[-grep("pr25",files)]
#}

#if(varname=="tmax" & length(grep("tmax95",files))>0){
#files = files[-grep("tmax95",files)]
#}

#if(varname=="tmin" & length(grep("tmin32",files))>0){
#files = files[-grep("tmin32",files)]
#}

# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file

fittedvals = list()
resvals = list()

for(f in 1:length(files)){

test = nc_open(files[f])

if(f==1){
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
time = ncvar_get(test,"time")
}

fitv = ncvar_get(test,"fittedvalues")
resv = ncvar_get(test,"residuals")

fittedvals[[f]] = aperm(apply(fitv,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))
resvals[[f]] = aperm(apply(resv,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

nc_close(test)

}

fittedvalarray = resvalarray = array(NA,dim=c(length(lon),length(lat),length(time),length(files)))

for(f in 1:length(files)){
fittedvalarray[,,,f] = fittedvals[[f]]
resvalarray[,,,f] = resvals[[f]]
}

rm(fittedvals)
rm(resvals)

############
# Load GCM and DS weights

test = nc_open(paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_GCMweights.nc",sep=""))
GCMweights = ncvar_get(test,"GCMweights")
nc_close(test)

test = nc_open(paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_DSweights.nc",sep=""))
DSweights = ncvar_get(test,"DSweights")
nc_close(test)

#############
# Get unique GCM, DS, and scenario names from filenames

filesplit = do.call(rbind,strsplit(files,"_",fixed=TRUE))
GCMs = unique(filesplit[,2])
DSs = unique(paste(filesplit[,3],filesplit[,4],sep="_"))
scens = unique(filesplit[,5])

##############
# Calculate V - Internal Variability

V = intvar_sp(resvalarray,weights=list(GCMweights,DSweights),GCMs,DSs,files=files)

##############
# Calculate M - GCM Uncertainty

M = GCMunc(fittedvalarray,weights=GCMweights,scens,DSs,GCMs,files=files)

##############
# Calculate D - DS Uncertainty

D = DSunc(fittedvalarray,weights=DSweights,scens,DSs,GCMs,files=files)

##############
# Calculate S - Scenario Uncertainty

S = scenunc(fittedvalarray,scens,DSs,GCMs,weightsGCM=GCMweights,weightsDS=DSweights,files=files)

##############
# Smoothed values of M, D, and S

M_smooth = M #aperm(apply(M,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))
D_smooth = D #aperm(apply(D,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))
S_smooth = S #aperm(apply(S,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

#####
# map check

#testsfc1=list(x=lon,y=lat,z=V)
#testsfc2=list(x=lon,y=lat,z=M_smooth[,,100])
#testsfc3=list(x=lon,y=lat,z=D_smooth[,,100])
#testsfc4=list(x=lon,y=lat,z=S_smooth[,,100])

#par(mfrow=c(2,2))
#surface(testsfc1,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)
#surface(testsfc2,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)
#surface(testsfc3,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)
#surface(testsfc4,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)

###############
# Calculate mean change

G=meanchange(fittedvalarray,scens,GCMs,DSs,GCMweights,DSweights,files)
G_smooth = G #aperm(apply(G,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

################
# Calculate Total Variance

T = array(NA,dim=dim(M))
for(t in 1:dim(M)[3]) T[,,t] = M[,,t]+D[,,t]+S[,,t]+V
T_smooth = T #aperm(apply(T,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

######
# map check 2

#testsfc1=list(x=lon,y=lat,z=G_smooth[,,100])
#testsfc2=list(x=lon,y=lat,z=T_smooth[,,100])

#par(mfrow=c(2,1))
#surface(testsfc1,type="I")
#map("state",add=TRUE)
#surface(testsfc2,type="I")
#map("state",add=TRUE)

##################
# Write netcdf with components

if(varname=="tasmax" | varname=="tasmin"){
dataunits1 = "K^2"
dataunits2 = "K"
}

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25" | varname=="r1mm"){
dataunits1 = "days^2"
dataunits2 = "days"
}

if(varname=="pr"){
dataunits1 = "%^2"
dataunits2 = "%"
}

if(varname=="rx1day"){
  dataunits1 = "mm^2"
  dataunits2 = "mm"
}

if(varname=="rx5day"){
  dataunits1 = "mm^2"
  dataunits2 = "mm"
}

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )
dimT <- ncdim_def( "time","year",time)

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

fileout = paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_components_fix2.nc",sep="")

var1d <- ncvar_def("V", dataunits1, list(dimX,dimY), mv )
var2d <- ncvar_def("M", dataunits1, list(dimX,dimY,dimT), mv )
var3d <- ncvar_def("M_smooth", dataunits1, list(dimX,dimY,dimT), mv )
var4d <- ncvar_def("D", dataunits1, list(dimX,dimY,dimT), mv )
var5d <- ncvar_def("D_smooth", dataunits1, list(dimX,dimY,dimT), mv )
var6d <- ncvar_def("S", dataunits1, list(dimX,dimY,dimT), mv )
var7d <- ncvar_def("S_smooth", dataunits1, list(dimX,dimY,dimT), mv )
var8d <- ncvar_def("G", dataunits2, list(dimX,dimY,dimT), mv )
var9d <- ncvar_def("G_smooth", dataunits2, list(dimX,dimY,dimT), mv )
var10d <- ncvar_def("T", dataunits1, list(dimX,dimY,dimT), mv )
var11d <- ncvar_def("T_smooth", dataunits1, list(dimX,dimY,dimT), mv )


# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d,var11d) )
ncvar_put( nc1, var1d, V) # no start or count: write all values\
ncvar_put( nc1, var2d, M) # no start or count: write all values\
ncvar_put( nc1, var3d, M_smooth) # no start or count: write all values\
ncvar_put( nc1, var4d, D) # no start or count: write all values\
ncvar_put( nc1, var5d, D_smooth) # no start or count: write all values\
ncvar_put( nc1, var6d, S) # no start or count: write all values\
ncvar_put( nc1, var7d, S_smooth) # no start or count: write all values\
ncvar_put( nc1, var8d, G) # no start or count: write all values\
ncvar_put( nc1, var9d, G_smooth) # no start or count: write all values\
ncvar_put( nc1, var10d, T) # no start or count: write all values\
ncvar_put( nc1, var11d, T_smooth) # no start or count: write all values\

nc_close(nc1)

