source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(mapdata)
library(maptools)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(ggplot2)

setwd("/home/woot0002/DS_ind/")

var = varin = "pr"
type="ann"

if(var=="tmin"){
  varin="tasmin"
}
if(var=="tmax"){
  varin="tasmax"
}
if(var=="tmax95"){
  varin="tasmax"
}
if(var=="tmin32"){
  varin="tasmin"
}
if(var=="r1mm"){
  varin="pr"
}
if(var=="pr50"){
  varin="pr"
}
if(var=="rx1day"){
  varin="pr"
}
if(var=="rx5day"){
  varin="pr"
}

GCMfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)
LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*histclimo*.nc",sep=""),intern=TRUE)

GCMprojfiles = system(paste("ls /home/woot0002/GCMs/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)
LOCAprojfiles = system(paste("ls /home/woot0002/LOCA/regrid/",varin,"_*projclimo*.nc",sep=""),intern=TRUE)

LIVNEHfile = system(paste("ls /home/woot0002/monthlyclimo/",varin,"_day*livneh*.nc",sep=""),intern=TRUE)

load("/home/woot0002/DS_ind/GCMlist.Rdata")

GCM_hfiles = GCM_pfiles = LOCA_hfiles = LOCA_pfiles = c()

for(i in 1:length(GCMlist)){
  
  GCM_hfiles[i] = GCMfiles[grep(paste(GCMlist[i],"_",sep=""),GCMfiles)]
  GCM_pfiles[i] = GCMprojfiles[grep(paste(GCMlist[i],"_",sep=""),GCMprojfiles)]
  
  LOCA_hfiles[i] = LOCAfiles[grep(paste(GCMlist[i],"_",sep=""),LOCAfiles)]
  LOCA_pfiles[i] = LOCAprojfiles[grep(paste(GCMlist[i],"_",sep=""),LOCAprojfiles)]
  
}

###
# create full filelist + metadata table - historical

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_hfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMhdat = filelist2[,c(2,3,4,6)]
names(GCMhdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_hfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCAhdat = filelist2[,c(2,3,4,6)]
names(LOCAhdat) = names(GCMhdat)

#All metadata
GCM = rep(NA,1)
exp = rep(NA,1)
DS = rep(NA,1)
training = "LIVNEH"
obsdat = data.frame(GCM,exp,DS,training)

GCMhdat = rbind(GCMhdat,obsdat)
LOCAhdat= rbind(LOCAhdat,obsdat)

# all files
GCMgroup = c(GCM_hfiles,LIVNEHfile)
LOCAgroup = c(LOCA_hfiles,LIVNEHfile)

###
# create full filelist + metadata table - projected

#GCMs
filelist1 = do.call("rbind",strsplit(GCM_pfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "NA"
GCMpdat = filelist2[,c(2,3,4,6)]
names(GCMpdat) = c("GCM","exp","DS","training")

#LOCA
filelist1 = do.call("rbind",strsplit(LOCA_pfiles,"/",fixed=TRUE))
filelist2 = do.call("rbind",strsplit(filelist1[,6],"_",fixed=TRUE))
filelist2 = as.data.frame(filelist2)
filelist2$training = "Livneh"
LOCApdat = filelist2[,c(2,3,4,6)]
names(LOCApdat) = names(GCMpdat)

# all files
GCMpgroup = GCM_pfiles
LOCApgroup = LOCA_pfiles

### LOCA historical + Livneh

  nctest = nc_open(LOCAgroup[1])
  tmp = ncvar_get(nctest,"prclimo",start=c(1,1,1),count=c(-1,-1,1))
  lat = ncvar_get(nctest,"lat")
  lon=ncvar_get(nctest,"lon")
  nc_close(nctest)


#######

######
# state mask

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

statelist = c("oklahoma","texas","louisiana","new mexico")
for(s in 1:length(statelist)){
  
  regionname=statelist[s]
  states <- map_data("state")
  
  stateregion = subset(states,region==regionname)
  
  pointcheck = point.in.polygon(modelgrid$lon,modelgrid$lat,stateregion$long,stateregion$lat)
  modelgrid$mask = pointcheck
  regionmask = tmp
  for(R in 1:190){
    for(C in 1:140){
      regionmask[R,C] = modelgrid$mask[which(modelgrid$R==R & modelgrid$C==C)]
    }
  }
  out_filename = paste(statelist[s],"_mask.nc",sep="")
  
  dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def("mask","",longname="statemask", list(dimX,dimY), mv )
  
  #######
  # Create netcdf file
  
  nc <- nc_create(paste("/home/woot0002/DS_ind/",out_filename,sep="") ,  var1d )
  
  # Write some data to the file
  ncvar_put(nc, var1d, regionmask) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  
}


#library(ggplot2)
#ggplot(meansframe, aes(x=DS, y=histmeans)) + geom_point(aes(colour = factor(group)),size=5) +geom_hline(yintercept=mean(OBS[regionidx],na.rm=TRUE))
#ggplot(meansframe, aes(x=DS, y=changemeans)) + geom_point(aes(colour = factor(group)),size=5) +geom_hline(yintercept=0)
