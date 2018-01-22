################
#
# Uncertainty analysis gridded

############
# library and function load

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
anfilepath = "analysis"

###
# set extra options

printPDF = TRUE # if you want some plots printed to PDF set this to TRUE


##############
# Read in components from netcdf file

compfile = paste("uncertainty/SE/",anfilepath,"/",varname,"_components_fix.nc",sep="")

# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file

test = open.ncdf(compfile)

lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
time = 1950:2099

V = get.var.ncdf(test,"V")
M = get.var.ncdf(test,"M_smooth")
D = get.var.ncdf(test,"D_smooth")
S = get.var.ncdf(test,"S_smooth")
G = get.var.ncdf(test,"G_smooth")
T = get.var.ncdf(test,"T_smooth")
close.ncdf(test)

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

dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
dimT <- dim.def.ncdf( "time","year",time[timeidx])

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

fileout = paste("uncertainty/SE/",anfilepath,"/",varname,"_analysis_fix.nc",sep="")

var1d <- var.def.ncdf("FM", dataunits1, list(dimX,dimY,dimT), mv )
var2d <- var.def.ncdf("FD", dataunits1, list(dimX,dimY,dimT), mv )
var3d <- var.def.ncdf("FS", dataunits1, list(dimX,dimY,dimT), mv )
var4d <- var.def.ncdf("FV", dataunits1, list(dimX,dimY,dimT), mv )
var5d <- var.def.ncdf("FT", dataunits1, list(dimX,dimY,dimT), mv )

var6d <- var.def.ncdf("SM", dataunits1, list(dimX,dimY,dimT), mv )
var7d <- var.def.ncdf("SD", dataunits1, list(dimX,dimY,dimT), mv )
var8d <- var.def.ncdf("SS", dataunits2, list(dimX,dimY,dimT), mv )
var9d <- var.def.ncdf("SV", dataunits2, list(dimX,dimY,dimT), mv )
var10d <- var.def.ncdf("ST", dataunits1, list(dimX,dimY,dimT), mv )

var11d <- var.def.ncdf("CM", dataunits1, list(dimX,dimY,dimT), mv )
var12d <- var.def.ncdf("CD", dataunits1, list(dimX,dimY,dimT), mv )
var13d <- var.def.ncdf("CS", dataunits1, list(dimX,dimY,dimT), mv )
var14d <- var.def.ncdf("CV", dataunits1, list(dimX,dimY,dimT), mv )


# Create the test file
# Write some data to the file
# close ncdf

nc1 <- create.ncdf(fileout , list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d,var11d,var12d,var13d,var14d) )

FM = ifelse(FM==Inf,NA,FM)
FV = ifelse(FV==Inf,NA,FV)
FD = ifelse(FD==Inf,NA,FD)
FS = ifelse(FS==Inf,NA,FS)
FT = ifelse(FT==Inf,NA,FT)

put.var.ncdf( nc1, var1d, FM) # no start or count: write all values\
put.var.ncdf( nc1, var2d, FD) # no start or count: write all values\
put.var.ncdf( nc1, var3d, FS) # no start or count: write all values\
put.var.ncdf( nc1, var4d, FV) # no start or count: write all values\
put.var.ncdf( nc1, var5d, FT) # no start or count: write all values\

put.var.ncdf( nc1, var6d, SM) # no start or count: write all values\
put.var.ncdf( nc1, var7d, SD) # no start or count: write all values\
put.var.ncdf( nc1, var8d, SS) # no start or count: write all values\
put.var.ncdf( nc1, var9d, SV) # no start or count: write all values\
put.var.ncdf( nc1, var10d, ST) # no start or count: write all values\

put.var.ncdf( nc1, var11d, CM) # no start or count: write all values\
put.var.ncdf( nc1, var12d, CD) # no start or count: write all values\
put.var.ncdf( nc1, var13d, CS) # no start or count: write all values\
put.var.ncdf( nc1, var14d, CV) # no start or count: write all values\

close.ncdf(nc1)

##############33
# PDF contribution plots - if printPDF==TRUE

if(printPDF==TRUE){

earlyidx = which(time[timeidx]==2020)
mididx = which(time[timeidx]==2055)
lateidx = which(time[timeidx]==2090)

testsfcCV1 = list(x=lon,y=lat,z=CV[,,earlyidx])
testsfcCV2 = list(x=lon,y=lat,z=CV[,,mididx])
testsfcCV3 = list(x=lon,y=lat,z=CV[,,lateidx])

testsfcCM1 = list(x=lon,y=lat,z=CM[,,earlyidx])
testsfcCM2 = list(x=lon,y=lat,z=CM[,,mididx])
testsfcCM3 = list(x=lon,y=lat,z=CM[,,lateidx])

testsfcCS1 = list(x=lon,y=lat,z=CS[,,earlyidx])
testsfcCS2 = list(x=lon,y=lat,z=CS[,,mididx])
testsfcCS3 = list(x=lon,y=lat,z=CS[,,lateidx])

testsfcCD1 = list(x=lon,y=lat,z=CD[,,earlyidx])
testsfcCD2 = list(x=lon,y=lat,z=CD[,,mididx])
testsfcCD3 = list(x=lon,y=lat,z=CD[,,lateidx])

zlims = c(0,100)
brks = seq(0,100,by=5)
colbar = colorRampPalette(c("darkblue","lightblue","green","yellow","orange","red","darkred"))(length(brks)-1)

pdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_contributionplots_fix.pdf",sep=""),height=10,width=17,onefile=TRUE)

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

