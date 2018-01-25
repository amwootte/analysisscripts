################
#
# Uncertainty analysis gridded

############
# library and function load

############
# load libraries

library(ncdf4)
library(fields)
library(maps)

############
# load calculator functions

source("/home/woot0002/scripts/uncertainty/componentfunctions_2018.R")

###
# set variable name

varname = "tasmax"
anfilepath = "analysis"

###
# set extra options

printPDF = TRUE # if you want some plots printed to PDF set this to TRUE

##############
# Read in components from netcdf file

compfile = paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_components_fix2.nc",sep="")

# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file

test = nc_open(compfile)

lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
time = 1981:2099

V = ncvar_get(test,"V")
M = ncvar_get(test,"M_smooth")
D = ncvar_get(test,"D_smooth")
S = ncvar_get(test,"S_smooth")
G = ncvar_get(test,"G_smooth")
T = ncvar_get(test,"T_smooth")
nc_close(test)

################
# Subset to future period only

timeidx = which(time>=2006)

M=M[,,timeidx]
D=D[,,timeidx]
S=S[,,timeidx]
G=G[,,timeidx]
T=T[,,timeidx]

################
# Fractional Uncertainty Calcs

FM = FD = FS = FV = FT = M

for(i in 1:dim(M)[3]){
FM[,,i] = (1.65*sqrt(M[,,i]))/G[,,i]
FD[,,i] = (1.65*sqrt(D[,,i]))/G[,,i]
FS[,,i] = (1.65*sqrt(S[,,i]))/G[,,i]
FV[,,i] = (1.65*sqrt(V))/G[,,i]
FT[,,i] = (1.65*sqrt(T[,,i]))/G[,,i]
}

#################
# Signal to Noise Ratio Calcs

SM = 1/FM
SD = 1/FD
SS = 1/FS
SV = 1/FV
ST = 1/FT

################
# Contribution to the Total Uncertainty per component calc

CM = CD = CS = CV = M

for(i in 1:dim(M)[3]){
CM[,,i] = (M[,,i]/T[,,i])*100
CD[,,i] = (D[,,i]/T[,,i])*100
CS[,,i] = (S[,,i]/T[,,i])*100
CV[,,i] = (V/T[,,i])*100
}

###############
# Write NetCDF with analysis output

dataunits1=dataunits2=""

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )
dimT <- ncdim_def( "time","year",time[timeidx])

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

fileout = paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_analysis_fix2.nc",sep="")

var1d <- ncvar_def("FM", dataunits1, list(dimX,dimY,dimT), mv )
var2d <- ncvar_def("FD", dataunits1, list(dimX,dimY,dimT), mv )
var3d <- ncvar_def("FS", dataunits1, list(dimX,dimY,dimT), mv )
var4d <- ncvar_def("FV", dataunits1, list(dimX,dimY,dimT), mv )
var5d <- ncvar_def("FT", dataunits1, list(dimX,dimY,dimT), mv )

var6d <- ncvar_def("SM", dataunits1, list(dimX,dimY,dimT), mv )
var7d <- ncvar_def("SD", dataunits1, list(dimX,dimY,dimT), mv )
var8d <- ncvar_def("SS", dataunits2, list(dimX,dimY,dimT), mv )
var9d <- ncvar_def("SV", dataunits2, list(dimX,dimY,dimT), mv )
var10d <- ncvar_def("ST", dataunits1, list(dimX,dimY,dimT), mv )

var11d <- ncvar_def("CM", dataunits1, list(dimX,dimY,dimT), mv )
var12d <- ncvar_def("CD", dataunits1, list(dimX,dimY,dimT), mv )
var13d <- ncvar_def("CS", dataunits1, list(dimX,dimY,dimT), mv )
var14d <- ncvar_def("CV", dataunits1, list(dimX,dimY,dimT), mv )


# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d,var11d,var12d,var13d,var14d) )

FM = ifelse(FM==Inf,NA,FM)
FV = ifelse(FV==Inf,NA,FV)
FD = ifelse(FD==Inf,NA,FD)
FS = ifelse(FS==Inf,NA,FS)
FT = ifelse(FT==Inf,NA,FT)

ncvar_put( nc1, var1d, FM) # no start or count: write all values\
ncvar_put( nc1, var2d, FD) # no start or count: write all values\
ncvar_put( nc1, var3d, FS) # no start or count: write all values\
ncvar_put( nc1, var4d, FV) # no start or count: write all values\
ncvar_put( nc1, var5d, FT) # no start or count: write all values\

ncvar_put( nc1, var6d, SM) # no start or count: write all values\
ncvar_put( nc1, var7d, SD) # no start or count: write all values\
ncvar_put( nc1, var8d, SS) # no start or count: write all values\
ncvar_put( nc1, var9d, SV) # no start or count: write all values\
ncvar_put( nc1, var10d, ST) # no start or count: write all values\

ncvar_put( nc1, var11d, CM) # no start or count: write all values\
ncvar_put( nc1, var12d, CD) # no start or count: write all values\
ncvar_put( nc1, var13d, CS) # no start or count: write all values\
ncvar_put( nc1, var14d, CV) # no start or count: write all values\

nc_close(nc1)

##############33
# PDF contribution plots - if printPDF==TRUE

if(printPDF==TRUE){

earlyidx = which(time[timeidx]==2020)
mididx = which(time[timeidx]==2055)
lateidx = which(time[timeidx]==2090)

testsfcCV1 = list(x=lon-360,y=lat,z=CV[,,earlyidx])
testsfcCV2 = list(x=lon-360,y=lat,z=CV[,,mididx])
testsfcCV3 = list(x=lon-360,y=lat,z=CV[,,lateidx])

testsfcCM1 = list(x=lon-360,y=lat,z=CM[,,earlyidx])
testsfcCM2 = list(x=lon-360,y=lat,z=CM[,,mididx])
testsfcCM3 = list(x=lon-360,y=lat,z=CM[,,lateidx])

testsfcCS1 = list(x=lon-360,y=lat,z=CS[,,earlyidx])
testsfcCS2 = list(x=lon-360,y=lat,z=CS[,,mididx])
testsfcCS3 = list(x=lon-360,y=lat,z=CS[,,lateidx])

testsfcCD1 = list(x=lon-360,y=lat,z=CD[,,earlyidx])
testsfcCD2 = list(x=lon-360,y=lat,z=CD[,,mididx])
testsfcCD3 = list(x=lon-360,y=lat,z=CD[,,lateidx])

zlims = c(0,100)
brks = seq(0,100,by=5)
colbar = colorRampPalette(c("darkblue","lightblue","green","yellow","orange","red","darkred"))(length(brks)-1)

pdf(paste("/home/woot0002/uncertainty/",anfilepath,"/",varname,"_contributionplots_fix2.pdf",sep=""),height=12,width=12,onefile=TRUE)

par(mfrow=c(3,4))

surface(testsfcCV1,type="I",xlab="Longitude",ylab="Latitude",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCM1,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCS1,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCD1,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)

surface(testsfcCV2,type="I",xlab="Longitude",ylab="Latitude",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCM2,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCS2,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCD2,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)

surface(testsfcCV3,type="I",xlab="Longitude",ylab="Latitude",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCM3,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCS3,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)
surface(testsfcCD3,type="I",xlab="Longitude",ylab="",zlim=zlims,col=colbar,breaks=brks)
map("state",add=TRUE)

dev.off()

}

