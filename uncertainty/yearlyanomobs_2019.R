####################
#
# Yearly Anomaly calculator
#
# This script will do the following with each dataset
#
# 1) Calculate Yearly averages (or sums) - also doing scale factor corrections if needed
# 2) Calculate 1981-2000 climatological average
# 3) Calculate Yearly anomalies (absolute or percent difference)
# 4) Combine historical and future runs into one file
# 5) Make variable, scenario, GCM, and DS names match to common naming in master lists
# 6) Write new yearly netcdfs containing the yearly values, climatology, and averages
# 
#####################
source("/data2/3to5/I35/scripts/analysisfunctions.R")

###
# Arguments set

datasetname = "METDATA" # common name of dataset to start with
varname = "tasmin" # common name of variable

###
# load libraries and master list files
library(ncdf4)
library(maps) # these two just to check plotting issues if need be
library(fields)

###
# determine files to use in calculation

# this is section is helping with the right file path on convection
if(datasetname=="METDATA") inname = "METDATA" 
if(datasetname=="PRISM") inname = "PRISM" 
if(datasetname=="cmip5.bcca") inname="cmip5_bcca"
if(datasetname=="regcm") inname="regcmdata"

if(varname=="tasmax" | varname=="tasmin" | varname=="pr"){
invar = invar2 = varname
if(datasetname=="METDATA" & varname=="tasmax"){ invar="tmmx"; invar2 = "tasmax";}
if(datasetname=="METDATA" & varname=="tasmin"){ invar="tmmn"; invar2 = "tasmin"; message("if worked right");}
if(datasetname=="METDATA" & varname=="pr"){ invar = "pr"; invar2 = "pr";}
} else {
if(varname=="tmax95") invar="tasmax"
if(datasetname=="METDATA" & varname=="tmax95"){
  invar="tmmx"; invar2 = "tasmax";
} 
if(varname=="tmin32") invar="tasmin"
if(datasetname=="METDATA" & varname=="tmin32"){
  invar="tmmn"; invar2 = "tasmin";
} 
if(varname=="pr25"){
  invar="pr"; invar2 = "pr"
} 
if(varname=="r1mm"){
  invar="pr"; invar2 = "pr"
} 
if(varname=="pr50"){
  invar="pr"; invar2 = "pr"
} 
if(varname=="rx1day"){
  invar="pr"; invar2 = "pr"
} 
if(varname=="rx5day"){
  invar="pr"; invar2 = "pr"
} 
}

# find files to work with
files = system(paste("ls /data4/data/OBS/",inname,"/",invar2,"/",invar,"*.nc",sep=""),intern=TRUE)
files = files[3:27]
#if(datasetname=="PRISM") files=files[-grep("climo",files)]


###########
# Start calculating yearly averages for each file

outputdata = list()
yearlist = list()

for(f in 1:length(files)){
    
    test = nc_open(files[f])
    
    if(f==1){
      lon = ncvar_get(test,"lon")
      lat = ncvar_get(test,"lat")
      #lat = rev(lat)
      yeararray = array(NA,dim=c(length(lon),length(lat),length(files)))
    }
    time = ncvar_get(test,"day")
    timeunits = test$var[[1]]$dim[[3]]$units
    
    vardata = ncvar_get(test,test$var[[1]]$name)
    vardatause = array(NA,dim=c(length(lon),length(lat),length(time)))
    nc_close(test)
    for(d in 1:length(time)){
      tmp2 = t(vardata[,,d])
      tmp3 = tmp2[,length(lat):1]
      vardatause[,,d] = tmp3
      rm(tmp2)
      rm(tmp3)
      message("Finished data gather for day ",d," / ",length(time))
    }
    vardata= vardatause
    rm(vardatause)
    if(f==1) lat = rev(lat)
    
    startdate = as.Date(substr(timeunits,12,21))
    times = startdate+time
    years = unique(as.numeric(substr(times,1,4)))
    
    yearlist[[f]] = years
    
      if(varname=="tasmax" | varname=="tasmin"){
        
          tmp = apply(vardata,c(1,2),mean,na.rm=TRUE)
          yeararray[,,f] = tmp
        
      } else {
        if(varname=="tmax95") {
          temp = ifelse((vardata-273.15)>=35,1,0)
        }
        if(varname=="tmin32"){
          temp = ifelse((vardata-273.15)<=0,1,0)
        } 
        if(varname=="pr25"){
          message("Calculating pr25")
          temp = ifelse((vardata)>=25.4,1,0)
        } 
        if(varname=="r1mm") {
          temp = ifelse((vardata)>=1,1,0)
        }
        if(varname=="pr" | varname=="rx1day" | varname=="rx5day"){
          temp = vardata
        } 
        
        if(varname!="rx1day" & varname!="rx5day"){
          message("summing number of days")
          yeararray[,,f] = apply(temp,c(1,2),sum,na.rm=TRUE)
          if(varname=="pr") yeararray[,,f] = ifelse(yeararray[,,f]==0,NA,yeararray[,,f])
          yeararray[,,f] = ifelse(is.na(vardata[,,1])==FALSE,yeararray[,,f],NA)
        } else {
          message("taking the maximum")
          if(varname=="rx1day"){
          yeararray[,,f] = apply(temp,c(1,2),max,na.rm=TRUE)
          yeararray[,,f] = ifelse(is.na(vardata[,,1])==FALSE,yeararray[,,f],NA)
          }
          if(varname=="rx5day"){
            yeararray[,,f] = apply(temp,c(1,2),calcrollsum,size=5)
            yeararray[,,f] = ifelse(is.na(vardata[,,1])==FALSE,yeararray[,,f],NA)
          }
        }
      }
    message("Finished yearly averages for file ",f," / ",length(files))
    }
    
    outputdata[[1]] = yeararray
    rm(yeararray)
    
if(datasetname=="METDATA"){
years = 1981:2005
outputarray = array(NA,dim=c(length(lon),length(lat),length(years)))
for(y in 1:length(years)) outputarray[,,y]=outputdata[[1]][,,y]
outputdata = list()
outputdata[[1]]=outputarray
}

#if(varname=="pr") outputdata[[1]] = ifelse(outputdata[[1]]==0,NA,outputdata[[1]]/10)

rm(vardata)

#############
# Calculate historical climatology

outputclimodat = list()

for(i in 1:length(outputdata)){
outputclimodat[[i]] = apply(outputdata[[i]],c(1,2),mean,na.rm=TRUE)
message("Finished climo calcs for file ",i," / ",length(outputdata))
}

#############
# Calculate anomalies

outputanomdat = list()

for(i in 1:length(outputdata)){
outputanomdat[[i]] = outputdata[[i]]
for(t in 1:length(years)){

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25"  | varname=="rx1day" | varname=="rx5day"  | varname=="r1mm" | varname=="tasmax" | varname=="tasmin"){
outputanomdat[[i]][,,t] = outputdata[[i]][,,t]-outputclimodat[[i]]
}

if(varname=="pr"){
outputanomdat[[i]][,,t] = ((outputdata[[i]][,,t]-outputclimodat[[i]])/(outputclimodat[[i]]))*100
}
}
message("Finished anomaly calcs for file ",i," / ",length(outputdata))
}

#############
# Regrid yearly anomalies, climatology, and yearly values to match desired grid.

regrid=TRUE

if(regrid==TRUE){

test=nc_open("/home/woot0002/uncertainty/anomdat/old/tasmax_MPI-ESM-LR_EDQMP0_rcp85_anom.nc")
newlon = ncvar_get(test,"lon")
newlat = ncvar_get(test,"lat")
nc_close(test)

if(all(newlon>0)==TRUE) lon=lon+360

newmatrix = newmatrix2 = array(NA,dim=c(length(newlon),length(newlat),length(years)))

for(y in 1:length(years)){
  
  test2tavg = interp.surface.grid(list(x=lon,y=lat,z=outputdata[[1]][,,y]),grid.list=list(x=newlon,y=newlat))
  newmatrix[,,y] = matrix(test2tavg[[3]],nrow=length(newlon),ncol=length(newlat))
  test2tavg = interp.surface.grid(list(x=lon,y=lat,z=outputanomdat[[1]][,,y]),grid.list=list(x=newlon,y=newlat))
  newmatrix2[,,y] = matrix(test2tavg[[3]],nrow=length(newlon),ncol=length(newlat))
  
}

test2tavg = interp.surface.grid(list(x=lon,y=lat,z=outputclimodat[[1]]),grid.list=list(x=newlon,y=newlat))
newmatrix3 = matrix(test2tavg[[3]],nrow=length(newlon),ncol=length(newlat))

outputdata[[1]] = newmatrix
outputanomdat[[1]] = newmatrix2
outputclimodat[[1]] = newmatrix3
lon = newlon
lat = newlat

}

#############
# Write data with new file and variable names to uncertainty folder

if(varname=="tasmax" | varname=="tasmin"){
  dataunits1 = dataunits2 = "K"
}

if(varname=="tmax95" | varname=="tmin32" | varname=="pr25" | varname=="r1mm"){
  dataunits1 = "days"
  dataunits2 = "days"
}

if(varname=="pr"){
  dataunits1 = "mm"
  dataunits2 = "%"
}

if(varname=="rx1day" | varname=="rx5day"){
  dataunits1 = "mm"
  dataunits2 = "mm"
}

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )
dimT <- ncdim_def( "time","year",years)

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

for(i in 1:length(outputdata)){

fileout = paste("/home/woot0002/uncertainty/anomdat/",varname,"_",datasetname,"_anom.nc",sep="")

var1d <- ncvar_def(varname, dataunits1, list(dimX,dimY,dimT), mv )
var2d <- ncvar_def(paste(varname,"_climo",sep=""), dataunits1, list(dimX,dimY), mv )
var3d <- ncvar_def(paste(varname,"_anom",sep=""), dataunits2, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , list(var1d,var2d,var3d) )
ncvar_put( nc1, var1d, outputdata[[i]]) # no start or count: write all values\
ncvar_put( nc1, var2d, outputclimodat[[i]]) # no start or count: write all values\
ncvar_put( nc1, var3d, outputanomdat[[i]]) # no start or count: write all values\

nc_close(nc1)

message("Finished writing file ",i," / ",length(outputdata))
}




















 
