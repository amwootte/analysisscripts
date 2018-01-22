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

library(ncdf)
library(fields)
library(maps)

############
# load calculator functions

source("componentfunctions.R")

###
# set variable name

varname = "pr25"
regfilepath = "regfits"
anfilepath = "analysis"

############
# Load fitted values and residuals from all GCMs

# determine files
files = system(paste("ls uncertainty/SE/",regfilepath,"/",varname,"*.nc",sep=""),intern=TRUE)
obsidx = grep("PRISM",files)
files = files[-obsidx] # files with projection data

if(varname=="pr" & length(grep("pr25",files))>0){
files = files[-grep("pr25",files)]
}

if(varname=="tmax" & length(grep("tmax95",files))>0){
files = files[-grep("tmax95",files)]
}

if(varname=="tmin" & length(grep("tmin32",files))>0){
files = files[-grep("tmin32",files)]
}

# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file

fittedvals = list()
resvals = list()

for(f in 1:length(files)){

test = open.ncdf(files[f])

if(f==1){
lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
time = get.var.ncdf(test,"time")
}

fitv = get.var.ncdf(test,"fittedvalues")
resv = get.var.ncdf(test,"residuals")

fittedvals[[f]] = aperm(apply(fitv,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))
resvals[[f]] = aperm(apply(resv,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

close.ncdf(test)

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

test = open.ncdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_GCMweights.nc",sep=""))
GCMweights = get.var.ncdf(test,"GCMweights")
close.ncdf(test)

test = open.ncdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_DSweights.nc",sep=""))
DSweights = get.var.ncdf(test,"DSweights")
close.ncdf(test)

#############
# Get unique GCM, DS, and scenario names from filenames

filesplit = do.call(rbind,strsplit(files,"_",fixed=TRUE))
GCMs = unique(filesplit[,2])
DSs = unique(filesplit[,3])
scens = unique(filesplit[,4])

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

if(varname=="tmax" | varname=="tmin"){
dataunits1 = "C^2"
dataunits2 = "C"
}

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25"){
dataunits1 = "days^2"
dataunits2 = "days"
}

if(varname=="pr"){
dataunits1 = "%^2"
dataunits2 = "%"
}

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
dimT <- dim.def.ncdf( "time","year",time)

# Make varables of various dimensionality, for illustration purposes
mv <- -999 # missing value to use

fileout = paste("uncertainty/SE/",anfilepath,"/",varname,"_components_fix2.nc",sep="")

var1d <- var.def.ncdf("V", dataunits1, list(dimX,dimY), mv )
var2d <- var.def.ncdf("M", dataunits1, list(dimX,dimY,dimT), mv )
var3d <- var.def.ncdf("M_smooth", dataunits1, list(dimX,dimY,dimT), mv )
var4d <- var.def.ncdf("D", dataunits1, list(dimX,dimY,dimT), mv )
var5d <- var.def.ncdf("D_smooth", dataunits1, list(dimX,dimY,dimT), mv )
var6d <- var.def.ncdf("S", dataunits1, list(dimX,dimY,dimT), mv )
var7d <- var.def.ncdf("S_smooth", dataunits1, list(dimX,dimY,dimT), mv )
var8d <- var.def.ncdf("G", dataunits2, list(dimX,dimY,dimT), mv )
var9d <- var.def.ncdf("G_smooth", dataunits2, list(dimX,dimY,dimT), mv )
var10d <- var.def.ncdf("T", dataunits1, list(dimX,dimY,dimT), mv )
var11d <- var.def.ncdf("T_smooth", dataunits1, list(dimX,dimY,dimT), mv )


# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d,var11d) )
put.var.ncdf( nc1, var1d, V) # no start or count: write all values\
put.var.ncdf( nc1, var2d, M) # no start or count: write all values\
put.var.ncdf( nc1, var3d, M_smooth) # no start or count: write all values\
put.var.ncdf( nc1, var4d, D) # no start or count: write all values\
put.var.ncdf( nc1, var5d, D_smooth) # no start or count: write all values\
put.var.ncdf( nc1, var6d, S) # no start or count: write all values\
put.var.ncdf( nc1, var7d, S_smooth) # no start or count: write all values\
put.var.ncdf( nc1, var8d, G) # no start or count: write all values\
put.var.ncdf( nc1, var9d, G_smooth) # no start or count: write all values\
put.var.ncdf( nc1, var10d, T) # no start or count: write all values\
put.var.ncdf( nc1, var11d, T_smooth) # no start or count: write all values\

close.ncdf(nc1)

