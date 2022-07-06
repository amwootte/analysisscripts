######
# Mussels map, climate data with mussels locations


#ptm= Sys.time()
#.libPaths('/home/amwootte/R/x86_64-redhat-linux-gnu-library/3.3/')
library(iterators)
#library(parallel)
#library(foreach)
#library(doParallel)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

#cl <- makePSOCKcluster(20,outfile="")#makePSOCKcluster(20)
#registerDoParallel(cl)

##########
# Set input information

message("Setting inputs")

shapefile1 = "/home/woot0002/shapefiles/WBD_HU8_SanSabaLlano"
shapefile2 = "/home/woot0002/shapefiles/MusselStreams"

#shapefile = "/home/woot0002/shapefiles/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01"

varname = "pr" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2070,2099)

colorchoicechange = "browntogreen"
BINSchange = 20
USEFS = FALSE
colorchoicehist = "whitetogreen"
BINShist = 20

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

############

varin = varname
if(varname=="tmax85" | varname=="tmax90" | varname=="tmax95" | varname=="tmax100" | varname=="gsl") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="prcptot" | varname=="prmean" | varname=="prmed" | varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd" | varname=="R99" | varname=="R95" | varname=="R90"  | varname=="R99freq" | varname=="R95freq" | varname=="R90freq") varin="pr"

########
# Find all file names

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

projnotes = paste(projnotes,collapse=",")
histnotes = paste(histnotes,collapse=",")

histlist = paste(histfilelist,collapse=",")
projlist = paste(projfilelist,collapse=",")

###########

test = nc_open(paste("/data2/3to5/I35/all_mems/",varname,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_ann.nc",sep=""))
memhist = ncvar_get(test,"histmean")
memproj = ncvar_get(test,"projmeandiff")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

###
# ensemble mean precipitation change


test = nc_open("/data2/3to5/I35/ens_means/pr_ensmean_absolute_2070-2099_ann.nc")
prprojrcp85 = ncvar_get(test,"projmeandiff_rcp85")
prprojrcp26 = ncvar_get(test,"projmeandiff_rcp26")
nc_close(test)

###
# ensemble mean tasmax change

test = nc_open("/data2/3to5/I35/ens_means/tasmax_ensmean_absolute_2070-2099_ann.nc")
tasmaxprojrcp85 = ncvar_get(test,"projmeandiff_rcp85")
tasmaxprojrcp26 = ncvar_get(test,"projmeandiff_rcp26")
nc_close(test)


###
# get shapefile
message("opening URG shapefile")
test = readShapePoly(shapefile1)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

message("opening river shapefile")
test2 = readShapeLines(shapefile2)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")


######
# Mussels locations

locname = c("FM126","CR340","Nobles","Charlies","5mile","10mile","BoisD'Arc")
latloc = c(31.23869444,31.19016833,30.99433833,30.91122222,30.913,30.88311111,30.90107944)
lonloc = c(-98.60208333,-98.90301639,-99.29857583,-99.51830556,-99.70875,-99.63033333,-99.91475917)

loctable = data.frame(locname,latloc,lonloc)

#10 FM126 Mussel Site 31.23869444 -98.60208333
#9 CR340 Mussel Site 31.19016833 -98.90301639
#8 Nobles Mussel Site 30.99433833 -99.29857583
#7 Charlies Mussel Site 30.91122222 -99.51830556
#5 5 mile Mussel Site 30.913 -99.70875
#6 10 mile Mussel Site 30.88311111 -99.63033333
#4 Bois D'Arc Mussel Site 30.90107944 -99.91475917

gaugetable = gauges=read.table("/home/woot0002/usgsgaugelocations.csv",sep=",",header=TRUE)

# Streamflow gages
#1 San Saba Rv at Menard, TX Gaging Station 30.91913605 -99.78561401
#2 San Saba Rv nr Brady, TX Gaging Station 31.00392342 -99.26885986
#3 San Saba Rv at San Saba, TX Gaging Station 31.2131691 -98.71947479


locname = c("FM864","ClearCreek","RamboRoad")
latloc = c(30.83476387,30.90370955,30.92)
lonloc = c(-100.0943805,-99.92433729,-99.722222)

TNCtable = data.frame(locname,latloc,lonloc)

# TNC locations
#11 FM 864 crossing TNC Rating Site 30.83476387 -100.0943805
#12 Clear Creek at US HWY 190 TNC Rating Site 30.90370955 -99.92433729
#13 Rambo Road TNC Rating Site 30.92 -99.722222

######
# Weather Stations

locname = c("USC00411138","USC00415650","USC00411875","USC00418449")
latloc = c(31.73831,30.74765,31.82766,30.583)
lonloc = c(-98.94545,-99.23067,-99.43192,-100.6502)

WStable = data.frame(locname,latloc,lonloc)


######
##########
# Plot maps

pdf(paste("/home/woot0002/mussels_histprecip_and_sites2.pdf",sep=""),width=10,height=5,onefile=TRUE)
par(mfrow=c(1,2))
s=1
  
    tmpV = memhist[,,3]
    tmpdat = memhist[,,3]
    title = "PRISM Historical"
  
  modrasV = raster(t(tmpV)[length(lat):1,])
  #rm(tmpV)
  if(all(lon>0)){
    extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
  } else {
    extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
  }
  areas = eval(parse(text=paste("test@data$","HUC_8",sep="")))

    #diffcolorbar = colorramp(tmpdat,colorchoice=colorchoicehist,Blimit=BINShist,type="raw",use_fixed_scale = USEFS)
  
  #testsfc = list(x=lon-360,y=lat,z=tmpdat)
  #surface(testsfc,type="I",main="Sites overlaid with Gridded Precipitation Data",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  #plot(test,add=TRUE)
  #map("state",add=TRUE)
  
  
  tmpdat2 = tmpdat[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]
  diffcolorbar = colorramp(tmpdat2,colorchoice=colorchoicehist,Blimit=BINShist,type="raw",use_fixed_scale = TRUE,fixed_scale = c(400,900))
  testsfc = list(x=lon-360,y=lat,z=tmpdat)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Gridded Precipitation Data",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  
  tmpdat2 = tmpdat[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]/25.4
  diffcolorbar = colorramp(tmpdat2,colorchoice=colorchoicehist,Blimit=15,type="raw",use_fixed_scale = TRUE,fixed_scale = c(10,40))
  testsfc = list(x=lon-360,y=lat,z=tmpdat/25.4)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Gridded Precipitation Data",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  dev.off()
  

  pdf(paste("/home/woot0002/mussels_sites_prchange.pdf",sep=""),width=10,height=5,onefile=TRUE)
  par(mfrow=c(1,2))
  
  tmpdat2.1 = prprojrcp26[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]
  tmpdat2.2 = prprojrcp85[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]
  diffcolorbar = colorramp(c(tmpdat2.1,tmpdat2.2),colorchoice="browntogreen",Blimit=20,type="difference",use_fixed_scale = TRUE,c(-100,100))
  
  testsfc = list(x=lon-360,y=lat,z=prprojrcp26)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nPrecipitation Change RCP2.6 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  testsfc = list(x=lon-360,y=lat,z=prprojrcp85)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nPrecipitation Change RCP8.5 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  
  tmpdat2.1 = prprojrcp26[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]/25.4
  tmpdat2.2 = prprojrcp85[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]/25.4
  diffcolorbar = colorramp(c(tmpdat2.1,tmpdat2.2),colorchoice="browntogreen",Blimit=20,type="difference",use_fixed_scale = USEFS)
  
  testsfc = list(x=lon-360,y=lat,z=prprojrcp26/25.4)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nPrecipitation Change RCP2.6 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  testsfc = list(x=lon-360,y=lat,z=prprojrcp85/25.4)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nPrecipitation Change RCP8.5 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  
  dev.off()
  
  
  
  pdf(paste("/home/woot0002/mussels_sites_tasmaxchange.pdf",sep=""),width=10,height=5,onefile=TRUE)
  par(mfrow=c(1,2))
  
  tmpdat2.1 = tasmaxprojrcp26[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]
  tmpdat2.2 = tasmaxprojrcp85[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]
  diffcolorbar = colorramp(c(tmpdat2.1,tmpdat2.2),colorchoice="bluetored",Blimit=20,type="difference",use_fixed_scale = USEFS)
  
  testsfc = list(x=lon-360,y=lat,z=tasmaxprojrcp26)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nHigh Temperature Change RCP2.6 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  testsfc = list(x=lon-360,y=lat,z=tasmaxprojrcp85)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nHigh Temperature Change RCP8.5 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  
  tmpdat2.1 = tasmaxprojrcp26[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]*9/5
  tmpdat2.2 = tasmaxprojrcp85[which((lon-360)>=-100.7 & (lon-360)<=-98), which(lat>=29.5 & lat<=32) ]*9/5
  diffcolorbar = colorramp(c(tmpdat2.1,tmpdat2.2),colorchoice="bluetored",Blimit=20,type="difference",use_fixed_scale = USEFS)
  
  testsfc = list(x=lon-360,y=lat,z=tasmaxprojrcp26*9/5)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nHigh Temperature Change RCP2.6 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  testsfc = list(x=lon-360,y=lat,z=tasmaxprojrcp85*9/5)
  surface(testsfc,type="I",xlim=c(-100.7,-98),ylim=c(29.5,32),main="Sites overlaid with Mean Projected \nHigh Temperature Change RCP8.5 2070-2099",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test2,add=TRUE,col="blue")
  plot(test,add=TRUE,lwd=2)
  points(latloc~lonloc,data=loctable,pch=19,col="yellow")
  points(latitude~longitude,data=gaugetable,pch=19,col="orange")
  points(latloc~lonloc,data=TNCtable,pch=19,col="red")
  points(latloc~lonloc,data=WStable,pch=19,col="purple")
  legend("bottomright",legend=c("Mussel Site","Gaging Station","TNC Rating Site","Weather Stations"),bg="white",cex=0.8,col=c("yellow","orange","red","purple"),pch=19)
  
  dev.off()
  
  plot(NA,xlim=c(-100.7,-98),ylim=c(29.5,32))
  plot(test2,col="blue",add=TRUE)
  plot(test,add=TRUE,lwd=2)
  
  
  
