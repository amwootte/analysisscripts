####################
#
# METDATA historical baseline to 3^5 domain
# 
#####################
###
# Arguments set

datasetname = "METDATA" # common name of dataset to start with

###
# load libraries and master list files
library(ncdf4)
library(maps) # these two just to check plotting issues if need be
library(fields)
source("/home/woot0002/scripts/analysisfunctions.R")

###
# setting variable names

varname = "tmax100"

if(varname == "tasmax"){ varin="tmmx"; dataunits1 = "degrees_K";}

if(varname == "tasmin"){ varin="tmmn"; dataunits1 = "degrees_K";}
if(varname == "pr"){ varin="pr"; dataunits1 = "mm";}

if(varname == "tmax100" | varname == "tmax95"){ varin="tmmx"; dataunits1 = "days";}
if(varname == "tmin32" | varname == "tmin28"){ varin="tmmn"; dataunits1 = "days";}
if(varname == "pr25" | varname == "pr50" | varname == "mdrn" | varname == "rx1day" | varname == "rx5day" | varname == "cwd" | varname == "cdd"){ varin="pr"}

if(varname == "pr25" | varname == "pr50" | varname == "mdrn"  | varname == "cwd" | varname == "cdd"){ dataunits1="days"}
if(varname == "rx1day" | varname == "rx5day"){ dataunits1="mm"}

###########
# Start calculating yearly averages for each file

years = 1981:2005

for(y in 1:length(years)){
    test = nc_open(paste("/data4/data/DS_proj/METDATA/",varin,"_",years[y],".nc",sep=""))
  
      if(y==1){
        lon = ncvar_get(test,"lon")
        lat = ncvar_get(test,"lat")
        lat = rev(lat)
        yeararray = array(NA,dim=c(length(lon),length(lat),length(years)))
      }
    
    vardata = ncvar_get(test,test$var[[1]]$name)
    
    time = ncvar_get(test,"day")
    timeunits = test$var[[1]]$dim[[3]]$units
    nc_close(test)
    
    startdate = as.Date(substr(timeunits,12,21))
    times = startdate+time
    
    if(varname == "tasmax" | varname == "tasmin"){
    temp2 = apply(vardata,c(1,2),mean,na.rm=TRUE)
    }
    
    if(varname == "pr"){
      temp2 = apply(vardata,c(1,2),sum,na.rm=TRUE)
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
    }
    
    if(varname == "tmax95" | varname == "tmax100" | varname == "tmin32" | varname == "tmin28" | varname == "pr25" | varname == "pr50" | varname == "mdrn"){
      if(varname == "tmax100")  temp = ifelse(vardata>=310.928,1,0)
      if(varname == "tmax95")  temp = ifelse(vardata>=308.15,1,0)
      if(varname == "tmin32")  temp = ifelse(vardata<=273.15,1,0)
      if(varname == "tmin28")  temp = ifelse(vardata<=270.928,1,0)
      if(varname == "pr25")  temp = ifelse(vardata>=50.8,1,0)
      if(varname == "pr50")  temp = ifelse(vardata>=25.4,1,0)
      if(varname == "mdrn")  temp = ifelse(vardata>=0.254,1,0)
      temp2 = apply(temp,c(1,2),sum,na.rm=TRUE)
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      rm(temp)
    }
    
    if(varname=="rx1day"){
      temp2 = apply(vardata,c(1,2),max,na.rm=TRUE)
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
    }
    if(varname=="rx5day"){
      temp2 = apply(vardata,c(1,2),calcrollsum,size=5)
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
    }
    if(varname=="frd"){
      #test1 = apply(tempdata,c(1,2),lastfreeze,startdate=startdate,inputtimes=times)
      ptm = proc.time()
      temp2 = matrix(NA,nrow=length(lat),ncol=length(lon))
      for(r in 1:length(lat)){
        for(c in 1:length(lon)){
          message("Working on r ",r," and c ",c)
         temp2[r,c] = lastfreeze(vardata[r,c,],startdate=startdate,inputtimes=times[yearidx])
        }
      }
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      ptm.end = proc.time()
    }
    if(varname=="cdd"){
      temp2 = apply(vardata,c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=0.254,outtype="max")
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
    }
    if(varname=="cwd"){
      temp2 = apply(vardata,c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=0.254,outtype="max")
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
    }
    
    temp3 = t(temp2)
    yeararray[,,y] = temp3[,length(lat):1]
    
    #testsfc = list(x=lon,y=lat,z=temp3[,length(lat):1])
    #surface(testsfc,type="I")
    
    rm(vardata)
    rm(temp2)
    rm(temp3)
    message("Finished data grab for year: ",years[y])
    }
    
climovalues = apply(yeararray,c(1,2),mean,na.rm=TRUE)
rm(yeararray)

#############
# Regrid climatology

test=nc_open("/home/woot0002/uncertainty/anomdat/tasmax_MPI-ESM-LR_EDQMP0_rcp85_anom.nc")
newlon = ncvar_get(test,"lon")
newlat = ncvar_get(test,"lat")
temp = ncvar_get(test,"tasmax",start=c(1,1,1),count=c(-1,-1,1))
nc_close(test)

dmask = ifelse(is.na(temp)==FALSE,1,0)

if(all(newlon>0)==TRUE) lon=lon+360

  test2tavg = interp.surface.grid(list(x=lon,y=lat,z=climovalues),grid.list=list(x=newlon,y=newlat))
  newmatrix = matrix(test2tavg[[3]],nrow=length(newlon),ncol=length(newlat))
  newmatrix = ifelse(dmask==1,newmatrix,NA)
  
  testsfc = list(x=newlon,y=newlat,z=newmatrix)
  surface(testsfc,type="I")
  
lon = newlon
lat = newlat

#############
# Write data with new file and variable names to uncertainty folder

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat )

# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

fileout = paste("/data2/3to5/I35/METDATA/",varname,"_histclimo.nc",sep="")

var1d <- ncvar_def(varname, dataunits1, list(dimX,dimY), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , var1d)
ncvar_put( nc1, var1d, newmatrix) # no start or count: write all values\

nc_close(nc1)





















 
