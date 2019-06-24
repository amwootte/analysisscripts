source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

setwd("/home/woot0002/")

var = "cwd"

files = system(paste("ls ",var,"*OFvOH_KSK_0masks.nc",sep=""),intern=TRUE)

pdf(paste(var,"_KSK_changemaps.pdf",sep=""),onefile=TRUE,width=10,height=10)

for(i in 1:length(files)){

test = nc_open(files[i])
lon = ncvar_get(test,"lon")-360
lat = ncvar_get(test,"lat")
dpku = ncvar_get(test,"dpku")
dpks = ncvar_get(test,"dpks")
nc_close(test)

filesplit = strsplit(files[i],"_",fixed=TRUE)[[1]]

#if(filesplit[2]=="DvL"){
#  titlepart = "Daymet vs. Livneh"
#}

#if(filesplit[2]=="PvL"){
#  titlepart = "PRISM vs. Livneh"
#}

#if(filesplit[2]=="PvD"){
#  titlepart = "PRISM vs. Daymet"
#}

titlepart = paste(filesplit[2],filesplit[3],filesplit[4],sep=" ")

diffcolorbar = colorramp(c(dpku,dpks),colorchoice="redtowhite",type="raw",Blimit=20,use_fixed_scale = TRUE,fixed_scale=c(0,1))

testsfc1 = list(x=lon,y=lat,z=dpks)
surface(testsfc1,type="I",main=paste(var," ",titlepart," change signal K-S test probability",sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon,y=lat,z=dpku)
surface(testsfc1,type="I",main=paste(var," ",titlepart," change signal Kuiper test probability",sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

}

dev.off()



