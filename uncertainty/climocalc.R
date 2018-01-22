setwd("wicci/res15km")

dataset = "wicci"

files = system("ls *.nc",intern=TRUE)
files
#####################

library(ncdf)
library(sp)
library(fields)

# Select correct file based on the above checks

for(f in 1:length(files)){
  ptm1 = proc.time()
  
  message("Climo calcs began at ",Sys.time()," for file: ",files[f])
  
  splitstring = strsplit(files[f],"_",fixed=TRUE)
  
if(dataset=="PRISM"){
  modelname = substr(splitstring[[1]][3],1,nchar(splitstring[[1]][3])-3)
  } 

if(dataset=="dcp") {
  splitstring = strsplit(files[f],"-",fixed=TRUE)
  
  if(length(grep("gfdl",splitstring[[1]]))==0){
  modelname = splitstring[[1]][3]
  } else {
  modelname = splitstring[[1]][4]
}
}

if(dataset=="cmip5_bcca"){
  modelname = splitstring[[1]][3]
}

if(dataset=="wicci") {
  splitstring = strsplit(files[f],"-",fixed=TRUE)
  modelname = splitstring[[1]][3]
}


  
  #########
  # Data prep
  
  # read file
  test = open.ncdf(files[f])
  
    lon = get.var.ncdf(test,"lon")
    lat = get.var.ncdf(test,"lat")
    time = get.var.ncdf(test,"time")
    
    newmatrix = array(data=NA,dim=c(length(lon),length(lat),12))
      
  testvar = get.var.ncdf(test,modelname)
 
  dataunits = test$var[[1]]$units
  timeunits = test$var[[1]]$dim[[3]]$units

  if(dataset=="PRISM"){
  unitdate = substr(timeunits,nchar(timeunits)-15,nchar(timeunits)-6)
  } else {
  unitdate = substr(timeunits,12,nchar(timeunits))
  }

if(dataset == "cmip5_bcca") unitdate = substr(timeunits,11,nchar(timeunits))
if(dataset == "wicci") unitdate = substr(timeunits,12,nchar(timeunits)-8)

  startdate = as.Date(unitdate)+time[1]
  dates = seq(startdate,length.out=length(time),by="day")  
  missingval = 1E+20

  close.ncdf(test)
  
  dateidx = which(as.numeric(substr(dates,1,4))>=1981 & as.numeric(substr(dates,1,4))<=2010)
  
  dates = dates[dateidx]
  testvar = testvar[,,dateidx]

if(dataset=="cmip5_bcca" & (modelname=="tasmax" | modelname=="tasmin")) testvar=  ifelse(testvar <= -50,NA,testvar)
if(dataset=="cmip5_bcca" & modelname=="pr") testvar = ifelse(testvar<0,NA,testvar)
#if(dataset=="wicci" & modelname=="prcp") testvar = ifelse(testvar<0,NA,testvar)

  #########
  # Calculate monthly climatological means
  
  for(m in 1:12){
    dateidx2 = which(as.numeric(substr(dates,6,7))==m)

    years= unique(substr(dates[dateidx2],1,4))
    temp = array(data=NA,dim=c(length(lon),length(lat),length(years)))

    for(y in 1:length(years)){
	if(dataset=="PRISM" & modelname=="ppt") temp[,,y] = apply(testvar[,,dateidx2[which(substr(dates[dateidx2],1,4)==years[y])]],c(1,2),sum,na.rm=TRUE)
	if(dataset=="PRISM" & (modelname=="tmax" | modelname=="tmin")) temp[,,y] = apply(testvar[,,dateidx2[which(substr(dates[dateidx2],1,4)==years[y])]],c(1,2),mean,na.rm=TRUE)

	if(dataset=="dcp" & modelname=="pr") temp[,,y] = apply(testvar[,,dateidx2[which(substr(dates[dateidx2],1,4)==years[y])]],c(1,2),sum,na.rm=TRUE)
	if(dataset=="dcp" & (modelname=="tmax" | modelname=="tmin")) temp[,,y] = apply(testvar[,,dateidx2[which(substr(dates[dateidx2],1,4)==years[y])]],c(1,2),mean,na.rm=TRUE)

	if(dataset=="cmip5_bcca" & modelname=="pr") temp[,,y] = apply(testvar[,,dateidx2[which(substr(dates[dateidx2],1,4)==years[y])]],c(1,2),sum,na.rm=TRUE)
	if(dataset=="cmip5_bcca" & (modelname=="tasmax" | modelname=="tasmin")) temp[,,y] = apply(testvar[,,dateidx2[which(substr(dates[dateidx2],1,4)==years[y])]],c(1,2),mean,na.rm=TRUE)
   }

    newmatrix[,,m] = apply(temp,c(1,2),mean,na.rm=TRUE)
  }
  
if(dataset=="dcp" & modelname=="pr") newmatrix=ifelse(newmatrix<0,NA,newmatrix)
if(dataset=="cmip5_bcca" & modelname=="pr") newmatrix=ifelse(newmatrix<0,NA,newmatrix)
if(dataset=="dcp" & (modelname=="tmax" | modelname=="tmin")) newmatrix=ifelse(newmatrix== -999,NA,newmatrix)/100

  #########
  # Write netcdf files
  
    dimX <- dim.def.ncdf( "lon", "degrees_east", lon)
    dimY <- dim.def.ncdf( "lat", "degrees_north", lat )
    dimT <- dim.def.ncdf( "time","month",1:12)
    # Make varables of various dimensionality, for illustration purposes
    mv <- missingval # missing value to use
    
    var1d <- var.def.ncdf(modelname, dataunits, list(dimX,dimY,dimT), mv )
    
    fileout = paste(strsplit(files[f],".",fixed=TRUE)[[1]][1],"_climo.nc",sep="")
    
    nc <- create.ncdf(fileout , list(var1d) )
    
    # Write some data to the file
    put.var.ncdf( nc, var1d, newmatrix ) # no start or count: write all values\
   
    # close ncdf
    close.ncdf(nc)
  
  message("Finished climo calcs for file ",f," / ",length(files))
  ptm1end = proc.time()-ptm1
  message("Time = ",ptm1end[3]," secs")
  message("Calcs complete at ",Sys.time())
  
}
  
  
  
  
  
