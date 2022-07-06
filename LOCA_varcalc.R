source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(multiApply)

# set basic arguments
varname = "rx5day" # these are all required arguments for step 1
varunits = "mm"
seasonin = "ann"

varin = varname
if(varname=="tmax85" | varname=="tmax90" | varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd" | varname=="gsl2") varin="tasmin"
if(varname=="prcptot" | varname=="prmean" | varname=="prmed" | varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd" | varname=="R99" | varname=="R95" | varname=="R90"  | varname=="R99freq" | varname=="R95freq" | varname=="R90freq") varin="pr"

TC = FALSE 
TH = NA
cond=NA

if(varname == "tasmax"){
  TC=FALSE
  appfunc = "mean"
}

if(varname == "tasmin"){
  TC=FALSE
  appfunc = "mean"
}

if(varname == "pr"){
  TC=FALSE
  appfunc = "sum"
}

if(substr(varname,1,3) == "SPI"){
  TC=FALSE
  appfunc = varname
  varin = "pr"
}

if(varname == "heatwaves"){
  TC=FALSE
  appfunc = "heatwaves"
  varin = "tasmax"
}

if(varname == "gsl"){
  TC=FALSE
  appfunc = "growing_season_length"
  varin = "tasmax"
}

if(varname == "gsl2"){
  TC=FALSE
  appfunc = "gsl2"
  varin = "tasmin"
}

if(varname == "R90"){
  TC=FALSE
  appfunc = "R90"
  varin = "pr"
}

if(varname == "R95"){
  TC=FALSE
  appfunc = "R95"
  varin = "pr"
}

if(varname == "R99"){
  TC=FALSE
  appfunc = "R99"
  varin = "pr"
}

if(varname == "R90freq"){
  TC=FALSE
  appfunc = "R90freq"
  varin = "pr"
}

if(varname == "R95freq"){
  TC=FALSE
  appfunc = "R95freq"
  varin = "pr"
}

if(varname == "R99freq"){
  TC=FALSE
  appfunc = "R99freq"
  varin = "pr"
}

if(varname == "tmax85"){
  TC=TRUE
  TH = 302.59
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmax90"){
  TC=TRUE
  TH = 305.372
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmax95"){
  TC=TRUE
  TH = 308.15
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmax100"){
  TC=TRUE
  TH = 310.928
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmin32"){
  TC=TRUE
  TH = 273.15
  cond = "lte"
  appfunc = "sum"
}

if(varname == "tmin28"){
  TC=TRUE
  TH = 270.928
  cond = "lte"
  appfunc = "sum"
}

if(varname == "frd"){
  TC=FALSE
  appfunc = "lastfreeze"
}

if(varname == "prcptot"){
  TC=FALSE
  appfunc = "prcptot"
}

if(varname == "prmean"){
  TC=FALSE
  appfunc = "prmean"
}

if(varname == "prmed"){
  TC=FALSE
  appfunc = "prmed"
}

if(varname == "mdrn"){
  TC=TRUE
  TH = 0.254
  cond = "gte"
  appfunc = "sum"
}

if(varname == "r1mm"){
  TC=TRUE
  TH = 1
  cond = "gte"
  appfunc = "sum"
}

if(varname == "pr25"){
  TC=TRUE
  TH = 25.4
  cond = "gte"
  appfunc = "sum"
}

if(varname == "pr50"){
  TC=TRUE
  TH = 50.8
  cond = "gte"
  appfunc = "sum"
}

if(varname == "rx1day"){
  TC=FALSE
  appfunc = "max"
}

if(varname == "rx5day"){
  TC=FALSE
  appfunc = "rx5day"
}

if(varname == "cdd"){
  TC=FALSE
  appfunc = "maxdryspell"
}

if(varname == "cwd"){
  TC=FALSE
  appfunc = "maxwetspell"
}

######

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",varin,"/historical/*.nc",sep=""),intern=TRUE)

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesin[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  varinuse = paste(filesplit2[1],filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
  nctest = nc_open(filesin[i])
  
  if(i==1){
    dataunits = nctest$var[[1]]$units
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    time=ncvar_get(nctest,"time")
    timeunits = nctest$var[[1]]$dim[[3]]$units
    latunits = nctest$var[[1]]$dim[[2]]$units
    lonunits = nctest$var[[1]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    timeout= time
    
    years = unique(substr(dates,1,4))
    timeoutyears = timeout[which(substr(dates,6,10)=="01-01")]
    
  }
  
  tmp = array(NA,dim=c(length(lon),length(lat),length(years)))
  start_time <- Sys.time()
  for(y in 1:length(years)){
    idx = which(substr(dates,1,4)==years[y])
    testdat = ncvar_get(nctest,varinuse,start=c(1,1,idx[1]),count=c(-1,-1,length(idx))) #length(idx)
    if(mean(testdat,na.rm=TRUE)<1 & varin=="pr"){
      testdat = testdat*86400
    }
    if((varin=="tasmax" | varin=="tasmin") & mean(testdat,na.rm=TRUE)<200){
      testdat = testdat+273.15
    }
    
    domainmask = ifelse(is.na(testdat[,,1])==FALSE,1,0)
    
    if(TC==TRUE){
      message("Doing threshold calculations")
      message("The threshold being used is ",TH)
      message("The condition is ",cond)
      
      if(cond=="gte"){
        testdat=ifelse(testdat >= as.numeric(TH),1,0)
      } 
      if(cond=="gt"){
        testdat=ifelse(testdat > as.numeric(TH),1,0)
      } 
      if(cond=="lte"){
        testdat=ifelse(testdat <= as.numeric(TH),1,0)
      } 
      if(cond=="lt"){
        testdat=ifelse(testdat < as.numeric(TH),1,0)
      } 
    }
    
    if(appfunc=="max"){
      tmp1 = apply(testdat,c(1,2),max,na.rm=TRUE)
      tmp1 = ifelse(tmp1== -Inf | tmp1==Inf,NA,tmp1)
      tmp[,,y] = ifelse(domainmask==1,tmp1,NA)
    }
    
    if(appfunc=="sum") {
      tmp1 = apply(testdat,c(1,2),sum,na.rm=TRUE)
      tmp[,,y] = ifelse(domainmask==1,tmp1,NA)
      message("Min value = ",min(tmp1,na.rm=TRUE)," Max value = ",max(tmp1,na.rm=TRUE))
    }
    
    if(appfunc=="mean") {
      tmp1 = apply(testdat,c(1,2),mean,na.rm=TRUE)
      tmp[,,y] = ifelse(domainmask==1,tmp1,NA)
      message("Min value = ",min(tmp1,na.rm=TRUE)," Max value = ",max(tmp1,na.rm=TRUE))
    }
    
    if(appfunc=="rx5day"){
      tmp1 = apply(testdat,c(1,2),calcrollsum,size=5)
      tmp[,,y] = ifelse(domainmask==1,tmp1,NA)
    }
    
    if(appfunc=="maxdryspell"){
      tmp1 = apply(testdat,c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=1,outtype="max") # changed to 1 mm 11/4
      tmp[,,y] = ifelse(domainmask==1,tmp1,NA)
    }
    if(appfunc=="maxwetspell"){
      tmp1 = apply(testdat,c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=1,outtype="max") # changed threshold to 1 mm 11/4
      tmp[,,y] = ifelse(domainmask==1,tmp1,NA)
    }
    rm(testdat)
    gc()
    message("Finished with year ",y," / ",length(years))
  }
  nc_close(nctest)
  end_time <- Sys.time()
  end_time - start_time
  testdat=tmp
  
  filenamepart = paste(varname,filesplit2[2],filesplit2[3],filesplit2[4],sep="_")
  filenameout = paste("/data4/data/DS_proj/LOCA/",varname,"/historical/",filenamepart,sep="")
  varout = paste(varname,filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
  dimX <- ncdim_def( "lon", lonunits, lon)
  dimY <- ncdim_def( "lat", latunits, lat)
  dimT <- ncdim_def("time",timeunits,timeoutyears)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varout,varunits,longname=varout, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, testdat) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(testdat)
  gc()
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(filesin))
    
  }
  
 
 

