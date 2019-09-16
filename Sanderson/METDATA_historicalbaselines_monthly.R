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
library(mailR)
source("/home/woot0002/scripts/analysisfunctions.R")

###
# setting variable names

varname = "rx5day"

if(varname == "tasmax"){ varin2="tasmax"; varin="tmmx"; dataunits1 = "degrees_C";}

if(varname == "tasmin"){ varin2="tasmin"; varin="tmmn"; dataunits1 = "degrees_C";}
if(varname == "pr"){ varin2="pr"; varin="pr"; dataunits1 = "mm";}
if(varname == "frd"){varin="tmmn"; varin2="tasmin"; dataunits1="days"}
if(varname == "gsl"){varin="tmmx"; varin2="tasmax"; dataunits1="days"}


if(varname == "tmax100" | varname == "tmax95"){ varin="tmmx"; varin2="tasmax"; dataunits1 = "days";}
if(varname == "tmin32" | varname == "tmin28"){ varin="tmmn"; varin2="tasmin"; dataunits1 = "days";}
if(varname == "pr25" | varname == "pr50" | varname == "r1mm" | varname == "rx1day" | varname == "rx5day" | varname == "cwd" | varname == "cdd"){ varin2="pr"; varin="pr";}

if(varname == "pr25" | varname == "pr50" | varname == "r1mm"  | varname == "cwd" | varname == "cdd"){ dataunits1="days"}
if(varname == "rx1day" | varname == "rx5day"){ dataunits1="mm"}

###########
# Start calculating yearly averages for each file

years = 1981:2005

for(y in 1:length(years)){
    test = nc_open(paste("/data4/data/OBS/METDATA/",varin2,"/",varin,"_",years[y],".nc",sep=""))
  
      if(y==1){
        lon = ncvar_get(test,"lon")
        lat = ncvar_get(test,"lat")
        lat = rev(lat)
        yeararray = array(NA,dim=c(length(lon),length(lat),length(years)*12))
      }
    
    vardata = ncvar_get(test,test$var[[1]]$name)
    
    time = ncvar_get(test,"day")
    timeunits = test$var[[1]]$dim[[3]]$units
    nc_close(test)
    
    startdate = as.Date(substr(timeunits,12,21))
    times = startdate+time
    for(m in 1:12){
      monidx = which(as.numeric(substr(times,6,7))==m)
      if(varname == "tasmax" | varname == "tasmin"){
        temp2 = apply(vardata[,,monidx],c(1,2),mean,na.rm=TRUE)-273.15
      }
      
      if(varname=="gsl"){ 
        test = nc_open(paste("/data4/data/DS_proj/METDATA/tmmn_",years[y],".nc",sep=""))
        vardata2 = ncvar_get(test,test$var[[1]]$name)
        nc_close(test)
        temp2 = matrix(NA,nrow=dim(vardata)[1],ncol=dim(vardata)[2])
        for(r in 1:length(lat)){
          for(c in 1:length(lon)){
            temp2[r,c] = GSLcalc(vardata[r,c,],vardata2[r,c,],inputtimes=time,startdate=startdate)
            #message("Finished Calculation R ",r," and C ",c)
          }
        }  
      }
    
      if(varname == "pr"){
        temp2 = apply(vardata[,,monidx],c(1,2),sum,na.rm=TRUE)
        temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      }
      
      if(varname == "tmax95" | varname == "tmax100" | varname == "tmin32" | varname == "tmin28" | varname == "pr25" | varname == "pr50" | varname == "r1mm"){
        if(varname == "tmax100")  temp = ifelse(vardata[,,monidx]>=310.928,1,0)
        if(varname == "tmax95")  temp = ifelse(vardata[,,monidx]>=308.15,1,0)
        if(varname == "tmin32")  temp = ifelse(vardata[,,monidx]<=273.15,1,0)
        if(varname == "tmin28")  temp = ifelse(vardata[,,monidx]<=270.928,1,0)
        if(varname == "pr25")  temp = ifelse(vardata[,,monidx]>=25.4,1,0)
        if(varname == "pr50")  temp = ifelse(vardata[,,monidx]>=50.8,1,0)
        if(varname == "r1mm")  temp = ifelse(vardata[,,monidx]>=1,1,0)
        temp2 = apply(temp,c(1,2),sum,na.rm=TRUE)
        temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
        rm(temp)
      }
      
      if(varname=="rx1day"){
        temp2 = apply(vardata[,,monidx],c(1,2),max,na.rm=TRUE)
        temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      }
      if(varname=="rx5day"){
        temp2 = apply(vardata[,,monidx],c(1,2),calcrollsum,size=5)
        temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      }
      if(varname=="frd"){
        #test1 = apply(tempdata,c(1,2),lastfreeze,startdate=startdate,inputtimes=times)
        ptm = proc.time()
        temp2 = matrix(NA,nrow=length(lat),ncol=length(lon))
        for(r in 1:length(lat)){
          for(c in 1:length(lon)){
            message("Working on r ",r," and c ",c)
            temp2[r,c] = lastfreeze(vardata[r,c,monidx],startdate=startdate,inputtimes=times)
          }
        }
        temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
        ptm.end = proc.time()
      }
      if(varname=="cdd"){
        temp2 = apply(vardata[,,monidx],c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=0.254,outtype="max")
        temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      }
      if(varname=="cwd"){
        temp2 = apply(vardata[,,monidx],c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=0.254,outtype="max")
        temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      }
      
      temp3 = t(temp2)
      idxout = m+(12*(y-1))
      yeararray[,,idxout] = temp3[,length(lat):1]
      
    }
    
    #testsfc = list(x=lon,y=lat,z=temp3[,length(lat):1])
    #surface(testsfc,type="I")
    
    rm(vardata)
    rm(temp2)
    rm(temp3)
    message("Finished data grab for year: ",years[y])
    }
    
monyear = seq(as.Date("1981-01-15"),as.Date("2005-12-15"),by="month")
climovalues = array(NA,dim=c(length(lon),length(lat),12))
for(m in 1:12){
  monidx = which(as.numeric(substr(monyear,6,7))==m)
  climovalues[,,m] = apply(yeararray[,,monidx],c(1,2),mean,na.rm=TRUE)
}

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

newmatrix = array(NA,dim=c(length(newlon),ncol=length(newlat),12))

for(m in 1:12){
  test2tavg = interp.surface.grid(list(x=lon,y=lat,z=climovalues[,,m]),grid.list=list(x=newlon,y=newlat))
  newvar = matrix(test2tavg[[3]],nrow=length(newlon),ncol=length(newlat))
  newvar = ifelse(dmask==1,newvar,NA)
  newmatrix[,,m]=newvar
  
 
}

#testsfc = list(x=newlon,y=newlat,z=newmatrix[,,m])
#surface(testsfc,type="I")

#lon = newlon
#lat = newlat

#############
# Write data with new file and variable names to uncertainty folder

dimX <- ncdim_def( "lon", "degrees_east", newlon)
dimY <- ncdim_def( "lat", "degrees_north", newlat )
dimT <- ncdim_def( "time", "month", 1:12 )
# Make varables of various dimensionality, for illustration purposes
mv <- 1e20 # missing value to use

fileout = paste("/data2/3to5/I35/METDATA/",varname,"_histclimo_mon.nc",sep="")

var1d <- ncvar_def(varname, dataunits1, list(dimX,dimY,dimT), mv )

# Create the test file
# Write some data to the file
# close ncdf

nc1 <- nc_create(fileout , var1d)
ncvar_put( nc1, var1d, newmatrix) 

nc_close(nc1)

send.mail(from = "amwootte@ou.edu",
          to = "amwootte@ou.edu",
          subject = "Notification from R on climatedata",
          body = paste("METDATA climatology calculations complete for var: ",varname," \n ",Sys.time(),sep=""), 
          authenticate = TRUE,
          smtp = list(host.name = "smtp.office365.com", port = 587,
                      user.name = "amwootte@ou.edu", passwd = "UNc3r!2inty", tls = TRUE))



















 
