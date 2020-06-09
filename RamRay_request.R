
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(mailR)
library(raster)
library(rasterVis)
library(maptools)
library(ggplot2)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

# email settings
emadd = "amwootte@ou.edu"
pswd = "D0wnSc2l!ng"

# set basic arguments
varname = "tasmax" # these are all required arguments for step 1
varunits = "degrees_C"
seasonin = "ann"

# set regional calculation arguments 
lon = c(-101,-94)+360 # information needed for step 6 if regiontype = "box"
lat = c(33,35)
regiontype = "states_and_domain"
regionname = "texas"

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

############

###########
# 1. Data Gather and conversion
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

histidx = which(histfilebreakdown$obs=="Daymet")
projidx = which(projfilebreakdown$obs=="Daymet")

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

projnotes = projnotes[projidx]
projfilelist = projfilelist[projidx]
projfilebreakdown = projfilebreakdown[projidx,]

histnotes = histnotes[histidx]
histfilelist = histfilelist[histidx]
histfilebreakdown = histfilebreakdown[histidx,]

#############
# make Texas mask

test = nc_open(histfilelist[1])
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

###
# create model grid

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360

print(summary(modelgrid$lon))

states <- map_data("state")

check = length(grep("-",regionname))
if(check==1){
  namesplit = strsplit(regionname,split="-")[[1]]
  stateregion = subset(states,region==paste(namesplit[1],namesplit[2],sep=" "))
} else {
  stateregion = subset(states,region==regionname)
}

pointcheck = point.in.polygon(modelgrid$lon,modelgrid$lat,stateregion$long,stateregion$lat)
modelgrid$mask = pointcheck
regionmask = matrix(NA,nrow=length(lon),ncol=length(lat))
for(R in 1:190){
  for(C in 1:140){
    regionmask[R,C] = modelgrid$mask[which(modelgrid$R==R & modelgrid$C==C)]
  }
}

testsfc=list(x=lon-360,y=lat,z=regionmask)
surface(testsfc,type="I")
map("state",add=TRUE)

#############

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  if(histfilebreakdown$DS[i]!="DeltaSD"){
    if(histfilebreakdown$GCM[i]=="MPI-ESM-LR"){
      noleap=FALSE
    } else {
      noleap=TRUE
    }
  } else {
    noleap=TRUE
  }  
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  #dates = datesin
  #print(dates[1:10])
  tempdata = array(NA,dim=c(length(lon),length(lat),length(datesin)))
  
  for(d in 1:length(datesin)){
    test = nc_open(histfilelist[i])
    tmp = ncvar_get(test,varname,start=c(1,1,d),count=c(-1,-1,1))
    tmp = ifelse(regionmask==1,tmp,NA)
    if(varname=="pr"){
      tempdata[,,d]=tmp*86400
    } else {
      tempdata[,,d]=tmp-273.15
    }
    nc_close(test)
    message("Finished date grab for date ",d," / ",length(datesin))
  }
    
    step1_filename = paste(varname,"_",regionname,"_",histfilebreakdown$DS[i],"_",histfilebreakdown$GCM[i],"_",histfilebreakdown$obs[i],"_",histfilebreakdown$scen[i],"_1981-2005.nc",sep="")
    
    dimX <- ncdim_def( "lon", "degrees_east", lon)
    dimY <- ncdim_def( "lat", "degrees_north", lat)
    dimT <- ncdim_def("time","days since 1981-01-01",as.numeric(datesin-as.Date("1981-01-01")))
    
    # Make varables of various dimensionality, for illustration purposes
    mv <- 1E20 # missing value to use
    
    var1d <- ncvar_def(varname,varunits, list(dimX,dimY,dimT), mv )
    
    #######
    # Create netcdf file
    nc <- nc_create(paste("/home/woot0002/RamRay/netcdfs/",step1_filename,sep="") ,  var1d )
    # Write some data to the file
    ncvar_put(nc, var1d, tempdata) # no start or count: write all values\
    # close ncdf
    nc_close(nc)
  rm(tempdata)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

###############
# 1b- Future calcs

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  message("Starting work on file ",projfilelist[i])
  if(projfilebreakdown$DS[i]!="DeltaSD"){
    if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
      noleap=FALSE
    } else {
      noleap=TRUE
    }
  } else {
    noleap=TRUE
  }  
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  #dates = datesin
  #print(dates[1:10])
  tempdata = array(NA,dim=c(length(lon),length(lat),length(datesin)))
  
  for(d in 1:length(datesin)){
    test = nc_open(projfilelist[i])
    tmp = ncvar_get(test,varname,start=c(1,1,d),count=c(-1,-1,1))
    tmp = ifelse(regionmask==1,tmp,NA)
    if(varname=="pr"){
      tempdata[,,d]=tmp*86400
    } else {
      tempdata[,,d]=tmp-273.15
    }
    nc_close(test)
    message("Finished date grab for date ",d," / ",length(datesin))
  }
  
  step1_filename = paste(varname,"_",regionname,"_",projfilebreakdown$DS[i],"_",projfilebreakdown$GCM[i],"_",projfilebreakdown$obs[i],"_",projfilebreakdown$scen[i],"_2006-2099.nc",sep="")
  
  dimX <- ncdim_def( "lon", "degrees_east", lon)
  dimY <- ncdim_def( "lat", "degrees_north", lat)
  dimT <- ncdim_def("time","days since 1981-01-01",as.numeric(datesin-as.Date("1981-01-01")))
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varname,varunits, list(dimX,dimY,dimT), mv )
  
  #######
  # Create netcdf file
  nc <- nc_create(paste("/home/woot0002/RamRay/netcdfs/",step1_filename,sep="") ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, tempdata) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(tempdata)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

send.mail(from = emadd,
          to = emadd,
          subject = "message from R on climatedata",
          body = paste("RamRay_request.R has finished running for ",varname,sep=""), 
          authenticate = TRUE,
          smtp = list(host.name = "smtp.office365.com", port = 587,
                      user.name = emadd, passwd = pswd, tls = TRUE))
