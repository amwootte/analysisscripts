###############
# 3^5 Projected Change Analyses - ens plotter

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##############
# User supplied inputs

varname = "tmax95" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname # don't change this
#noleap=TRUE # if your data has leap days change this to FALSE, otherwise leave it alone.  Should be false for 3^5.

if(varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd") varin="pr"

difftype="absolute" # type of difference to take, can be either absolute or percent

applymask = NA # if you don't want to apply a state mask, leave this alone

colorchoicediff = "yellowtored" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 30 # maximum number of color bins allowed for plotting the projected changes

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

###########
# 1. Data Gather and conversion

filename = paste("/data2/3to5/I35/ens_means/",varname,"_ensmean_absolute_2041-2070.nc",sep="")

test = nc_open(filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
projmean_rcp26 = ncvar_get(test,"projmean_rcp26")
projmean_rcp85 = ncvar_get(test,"projmean_rcp85")
nc_close(test)

################
# 2. Plotting

vals = as.vector(c(projmean_rcp26,projmean_rcp85))
vals = vals[-which(is.na(vals)==TRUE)]

lon = lon + 360

diffcolorbar = colorramp(vals,colorchoice=colorchoicediff,Blimit=20,type="raw")

pdf(paste("Group1_",varname,"_projvals.pdf",sep=""),onefile=TRUE,width=10,height=5)

par(mfrow=c(1,2))

  testsfc1 = list(x=lon,y=lat,z=projmean_rcp26)
  surface(testsfc1,type="I",main="Projected Climate\nScen: rcp26",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc1 = list(x=lon,y=lat,z=projmean_rcp85)
  surface(testsfc1,type="I",main="Projected Climate\nScen: rcp85",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)

dev.off()

