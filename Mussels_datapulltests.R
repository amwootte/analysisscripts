
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

shapefile = "/home/woot0002/shapefiles/WBD_HU8_SanSabaLlano"

varname = "tmax100" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2070,2099)

colorchoicechange = "purpletoorange"
BINSchange = 20
USEFS = FALSE
colorchoicehist = "whitetored"
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

test = nc_open(paste("/data2/3to5/I35/ens_means/",varname,"_ensmean_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_ann.nc",sep=""))
meanproj85 = ncvar_get(test,"projmeandiff_rcp85")
meanproj45 = ncvar_get(test,"projmeandiff_rcp45")
meanproj26 = ncvar_get(test,"projmeandiff_rcp26")
nc_close(test)

rcp26idx = which(projfilebreakdown$scen=="rcp26")
rcp45idx = which(projfilebreakdown$scen=="rcp45")
rcp85idx = which(projfilebreakdown$scen=="rcp85")

message("Gathered file list")


message("opening URG shapefile")
test = readShapePoly(shapefile)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

test@data$hist=NA
test@data$rcp26=NA
test@data$rcp45=NA
test@data$rcp85=NA

######
##########
# Calculate area weighted mean for each HRU.

pdf(paste("/home/woot0002/mussels_",varname,".pdf",sep=""),width=10,height=5,onefile=TRUE)
par(mfrow=c(1,2))
for(s in 1:4){
  
  if(s==1){
    tmpV = memhist[,,3]
    tmpdat = memhist[,,3]
    title = "PRISM Historical"
  }
  
  if(s==2){
    tmpV = meanproj26
    tmpdat = meanproj26
    title = "2070-2099 RCP2.6 Projected Change"
  }
  
  if(s==3){
    tmpV = meanproj45
    tmpdat = meanproj45
    title = "2070-2099 RCP4.5 Projected Change"
  }
  
  if(s==4){
    tmpV = meanproj85
    tmpdat = meanproj85
    title = "2070-2099 RCP8.5 Projected Change"
  }
  
  modrasV = raster(t(tmpV)[length(lat):1,])
  #rm(tmpV)
  if(all(lon>0)){
    extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
  } else {
    extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
  }
  areas = eval(parse(text=paste("test@data$","HUC_8",sep="")))
  
  #plot(modrasV)
  #plot(test,add=TRUE)
  #map("state",add=TRUE)
  if(s>1){
    diffcolorbar = colorramp(c(meanproj85,meanproj45,meanproj26),colorchoice=colorchoicechange,Blimit=BINSchange,use_fixed_scale = USEFS)
  } else {
    diffcolorbar = colorramp(tmpdat,colorchoice=colorchoicehist,Blimit=BINShist,type="raw",use_fixed_scale = USEFS)
  }
  
  testsfc = list(x=lon-360,y=lat,z=tmpdat)
  surface(testsfc,type="I",main=paste(title,varname,sep=" "),xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
  map("state",add=TRUE)
  
  surface(testsfc,type="I",xlim=c(-100.7,-98.4),ylim=c(29.8,31.5),main=paste(title,varname,sep=" "),xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
  plot(test,add=TRUE)
  
  tmp=c()
  
  for(a in 1:length(eval(parse(text=paste("test@data$","HUC_8",sep=""))))){ # total is 1021
    test.sub <- test[as.character(eval(parse(text=paste("test@data$","HUC_8",sep="")))) %in% areas[a], ] # in here are two important details. 1) column name in the data array, 2) item in the data array to subset by
    valueout=extract(modrasV,test.sub,weights=TRUE,fun=mean,na.rm=TRUE)  
    if(s==1){
      test@data$hist[a] = valueout
    } 
    if(s==2){
      test@data$rcp26[a] = valueout
    } 
    if(s==3){
      test@data$rcp45[a] = valueout
    } 
    if(s==4){
      test@data$rcp85[a] = valueout
    } 
  }
}

dev.off()

#plot(test,col=test@data$hist)
