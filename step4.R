########
# Step 4 analysis 3^5
#
source("analysisfunctions.R")
source("colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

# must pass in the following: step2_filename, varname, difftype, futureperiod,tempperiod,varunits,changeunits,projnotes
# colorchoicedata, colorchoicediff, BINLIMIT, obsbartype, diffbartype

test = nc_open(step2_filename)
rcp26dat = ncvar_get(test,"projmeandiff_rcp26")
rcp45dat = ncvar_get(test,"projmeandiff_rcp45")
rcp85dat = ncvar_get(test,"projmeandiff_rcp85")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

####

pdf(paste("/data2/3to5/I35/plots/ens_mem/Group1_",varname,"_",difftype,"_change.pdf",sep=""),onefile=TRUE,width=10,height=5)
diffcolorbar = diff_colorramp(c(rcp26dat,rcp45dat,rcp85dat),colorchoice=colorchoicediff,Blimit=BINLIMIT)

par(mfrow=c(1,2))

  testsfc1 = list(x=lon,y=lat,z=rcp26dat)
  surface(testsfc1,type="I",main="Projected Difference from Historical Climate\nScen: RCP2.6",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc1 = list(x=lon,y=lat,z=rcp85dat)
  surface(testsfc1,type="I",main="Projected Difference from Historical Climate\nScen: RCP8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)

dev.off()


