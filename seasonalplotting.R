##############################
#
# plot seasonal differences

library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("analysisfunctions.R")

##################

varname = "pr" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50

colorchoicediff = "browntogreen" # colorramps for difference plots, choices include "bluetored","redtoblue","browntogreen","greentobrown"
BINLIMIT = 25 # maximum number of color bins allowed for plotting the projected changes

##################

load("FileBreakdowns.Rdata")

load(file=paste("projchange_",varname,"_DJF_absolute.Rdata",sep=""))
DJFhist = histlist
DJFproj = projlist
DJFdiffs = diffs
DJFdiffsg1 = diffsg1

load(file=paste("projchange_",varname,"_MAM_absolute.Rdata",sep=""))
MAMhist = histlist
MAMproj = projlist
MAMdiffs = diffs
MAMdiffsg1 = diffsg1

load(file=paste("projchange_",varname,"_JJA_absolute.Rdata",sep=""))
JJAhist = histlist
JJAproj = projlist
JJAdiffs = diffs
JJAdiffsg1 = diffsg1

load(file=paste("projchange_",varname,"_SON_absolute.Rdata",sep=""))
SONhist = histlist
SONproj = projlist
SONdiffs = diffs
SONdiffsg1 = diffsg1

#################
scens = unique(projfilebreakdown$scen)

scensin = scens[c(1,3)]
DJFdiffsg1_sort = DJFdiffsg1[[3]][,,c(1,3)]
MAMdiffsg1_sort = MAMdiffsg1[[3]][,,c(1,3)]
JJAdiffsg1_sort = JJAdiffsg1[[3]][,,c(1,3)]
SONdiffsg1_sort = SONdiffsg1[[3]][,,c(1,3)]

diffcolorbar = diff_colorramp(c(DJFdiffsg1_sort,MAMdiffsg1_sort,JJAdiffsg1_sort,SONdiffsg1_sort),colorchoice=colorchoicediff,Blimit=BINLIMIT)

pdf(paste("Seasonal_Group1_",varname,"_absolute_projdiff.pdf",sep=""),onefile=TRUE,width=10,height=10)

par(mfrow=c(2,2))

for(i in 1:length(scensin)){
  testsfc1 = list(x=DJFdiffsg1[[1]],y=DJFdiffsg1[[2]],z=DJFdiffsg1_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference DJF \nScen: ",scensin[i],sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc1 = list(x=MAMdiffsg1[[1]],y=MAMdiffsg1[[2]],z=MAMdiffsg1_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference MAM \nScen: ",scensin[i],sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc1 = list(x=JJAdiffsg1[[1]],y=JJAdiffsg1[[2]],z=JJAdiffsg1_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference JJA \nScen: ",scensin[i],sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc1 = list(x=SONdiffsg1[[1]],y=SONdiffsg1[[2]],z=SONdiffsg1_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference SON \nScen: ",scensin[i],sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
}

dev.off()
