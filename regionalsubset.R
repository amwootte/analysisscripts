#################
# Subsetting 3^5 to URG Head Waters

library(ncdf4)
library(maps)
library(fields)
library(sp)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

# set arguments
varname = "tasmin" # these are all required arguments for step 1

lonbox = c(-107.5,-106.375)+360 # information needed for step 6 if regiontype = "box"
latbox = c(37.375,38.0)
regionname = "RGheadwater"

#####

if(varname == "tasmax" | varname == "tasmin"){
  appfunc = "mean"
}

if(varname == "pr"){
  appfunc = "sum"
}

#####

histfilelist = system(paste("ls /data2/3to5/I35/",varname,"/*/",varname,"_day_*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varname,"/*/",varname,"_day_*rcp*.nc",sep=""),intern=T)

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

#######

test = nc_open(histfilelist[1])
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

#####

Xextent = lonbox
Yextent = latbox

print(Xextent)
print(Yextent)

if(all(Xextent>0) & all(modelgrid$lon<0)) Xextent=Xextent-360

# get closest points in the grid which match the box
point1 = distfunc(Xextent[1],Yextent[1],modelgrid) #Xmin, Ymin
point2 = distfunc(Xextent[2],Yextent[1],modelgrid) #Xmax, Ymin
point3 = distfunc(Xextent[2],Yextent[2],modelgrid) #Xmax, Ymax
point4 = distfunc(Xextent[1],Yextent[2],modelgrid) #Xmin, Ymax

print(point1)
print(point3)

locstart = c(point1$R[1],point1$C[1])
locend = c(point3$R[1],point3$C[1])

##########

for(i in 1:length(histfilelist)){
test = nc_open(histfilelist[i])
tempdata = ncvar_get(test,varname,start=c(locstart[1],locstart[2],1),count=c((locend[1]-locstart[1])+1,(locend[2]-locstart[2])+1,-1))
lons = ncvar_get(test,"lon",start=locstart[1],count=(locend[1]-locstart[1])+1)
lats = ncvar_get(test,"lat",start=locstart[2],count=(locend[2]-locstart[2])+1)
times = ncvar_get(test,"time",start=1,count= -1)

varunits = test$var[[varname]]$units
lonunits = test$dim$lon$units
latunits = test$dim$lat$units
timeunits = test$dim$time$units
nc_close(test)

###

dimX <- ncdim_def( "lon",lonunits, lons)
dimY <- ncdim_def( "lat",latunits, lats)
dimT <- ncdim_def( "time",timeunits, times)

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

if(varname=="tasmax"){
  longvarname = "Daily High Temperature"
}
if(varname=="tasmin"){
  longvarname = "Daily Low Temperature"
}
if(varname=="pr"){
  longvarname = "Daily Total Precipitation"
}

var1d <- ncvar_def(varname,varunits,longname=longvarname, list(dimX,dimY,dimT), mv ,compression=9)

#######
# Create netcdf file

filesplit = do.call("c",strsplit(histfilelist[i],"/"))
filename = paste(substr(filesplit[length(filesplit)],1,nchar(filesplit[length(filesplit)])-3),"_",regionname,".nc",sep="")

nc <- nc_create(paste("/home/woot0002/RGheadwaters/",filename,sep="") ,  list(var1d) )

# Write some data to the file
ncvar_put(nc, var1d, tempdata) # no start or count: write all values\

# close ncdf
nc_close(nc)

message("Finished subset for file: ",histfilelist[i])

}

###################


for(i in 1:length(projfilelist)){
  test = nc_open(projfilelist[i])
  tempdata = ncvar_get(test,varname,start=c(locstart[1],locstart[2],1),count=c((locend[1]-locstart[1])+1,(locend[2]-locstart[2])+1,-1))
  lons = ncvar_get(test,"lon",start=locstart[1],count=(locend[1]-locstart[1])+1)
  lats = ncvar_get(test,"lat",start=locstart[2],count=(locend[2]-locstart[2])+1)
  times = ncvar_get(test,"time",start=1,count= -1)
  
  varunits = test$var[[varname]]$units
  lonunits = test$dim$lon$units
  latunits = test$dim$lat$units
  timeunits = test$dim$time$units
  nc_close(test)
  
  ###
  
  dimX <- ncdim_def( "lon",lonunits, lons)
  dimY <- ncdim_def( "lat",latunits, lats)
  dimT <- ncdim_def( "time",timeunits, times)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  if(varname=="tasmax"){
    longvarname = "Daily High Temperature"
  }
  if(varname=="tasmin"){
    longvarname = "Daily Low Temperature"
  }
  if(varname=="pr"){
    longvarname = "Daily Total Precipitation"
  }
  
  var1d <- ncvar_def(varname,varunits,longname=longvarname, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  filesplit = do.call("c",strsplit(projfilelist[i],"/"))
  filename = paste(substr(filesplit[length(filesplit)],1,nchar(filesplit[length(filesplit)])-3),"_",regionname,".nc",sep="")
  
  nc <- nc_create(paste("/home/woot0002/RGheadwaters/",filename,sep="") ,  list(var1d) )
  
  # Write some data to the file
  ncvar_put(nc, var1d, tempdata) # no start or count: write all values\
  
  # close ncdf
  nc_close(nc)
  
  message("Finished subset for file: ",projfilelist[i])
  
}

