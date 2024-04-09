#########################
#
# Analysis functions - a group of helper functions built up by AMW (University of Oklahoma)

source("/ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore02/3to5/I35/scripts/climdex.R")
source("/ourdisk/hpc/southcentralcasc/auto_archive_notyet/tape_2copies/tmp/cdstore02/3to5/I35/scripts/SPIcalcs.R")
####

distfunc = function(lon,lat,modelgrid){
  #model grid must have the lats and lons in the same format as the location lat and lon. modelgrid should also have r and c per grid point
  dist=3963 * acos(sin(lat/57.2958)*sin(modelgrid$lat/57.2958)+cos(lat/57.2958)*cos(modelgrid$lat/57.2958)*cos((modelgrid$lon/57.2958)-(lon/57.2958)))
  minidx = which(dist==min(dist,na.rm=TRUE))
  modelgrid[minidx,]
}

####

netcdftodailypointseries = function(filename,varname,dimnames=c("lon","lat","time"),datesdataperiod,loc_lat,loc_lon){
  
  #filename = histfilelist[i]
  #varname=varin
  #dimnames=c("lon","lat","time")
  #datesdataperiod=datesin
  #loc_lat=loc_lat
  #loc_lon=loc_lon
  yearlydat = c()
  #dates = datesdataperiod
  for(y in 1:length(datesdataperiod)){
    
    test = nc_open(filename)
    
    if(y==1){
      ###
      # create model grid
      lon = ncvar_get(test,dimnames[1])
      lat = ncvar_get(test,dimnames[2])
      LON = rep(lon,each=length(lat))
      LAT = rep(lat,length(lon))
      R = rep(1:length(lon),each=length(lat))
      C = rep(1:length(lat),length(lon))
      modelgrid = data.frame(R,C,LON,LAT)
      names(modelgrid) = c("R","C","lon","lat")
      if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
      
      ###
      # get cells to use
      loc_lon = as.numeric(loc_lon)
      loc_lat = as.numeric(loc_lat)
      if(loc_lon>0) loc_lon=loc_lon-360
      
      pointarea = distfunc(loc_lon,loc_lat,modelgrid)
      locstart = c(pointarea$R-1,pointarea$C-1)
    }
  
    tempdata = ncvar_get(test,varname,start=c(locstart[1],locstart[2],y),count=c(3,3,1))
    if(varname=="pr"){
      tempdata=tempdata*86400
      #message("Did pr unit conversion from mks to mm")  
      #message("Min value = ",min(tempdata,na.rm=TRUE)," Max value = ",max(tempdata,na.rm=TRUE))
    } 
    
    yearlydat[y] = mean(tempdata,na.rm=TRUE)
    nc_close(test)
    rm(tempdata)
    message("Finished Calcs for day ",datesdataperiod[y])
  }
  list(lon,lat,yearlydat)
}

  
####

netcdfcombocalcs = function(filename,varname="tasmax",dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),combofunction="growing_season_length"){
  #filename="/data2/3to5/I35/pr/DeltaSD/pr_day_I35prp1-DeltaSD-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
  #varname="pr"
  #dimnames=c("lon","lat","time")
  #yearlydataperiod=c(1981,2005)
  #datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  #datesdataperiod = datesdataperiod[-which(substr(datesdataperiod,6,10)=="02-29")]
  #combofunction="R99"
  
  file1 = filename
  #message(file1)
  #message(varname)
                            
  if(varname=="tasmax"){
    message("tasmax provided, figuring out correct tasmin file")
    #file1 = "/data2/3to5/I35/tasmax/PARM/tasmax_day_I35txdetrp1-PARM-B10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
    filesplit = do.call("c",strsplit(file1,"/",fixed=TRUE))
    filesplit2 = strsplit(filesplit[length(filesplit)],"_")[[1]]
    splitpart = strsplit(filesplit2[3],"-")[[1]]
    if(filesplit[6]=="PARM"){
      splitpart[3] = paste(substr(splitpart[3],1,1),"d",substr(splitpart[3],2,nchar(splitpart[3])),sep="")
    }
    findthis = paste("tasmin_day_*-",splitpart[2],"-",splitpart[3],"_",paste(filesplit2[4:length(filesplit2)],collapse="_"),sep="")
    file2 = system(paste("ls /",filesplit[2],"/",filesplit[3],"/",filesplit[4],"/tasmin/",filesplit[6],"/",findthis,sep=""),intern=TRUE)
    message("file 1 is tasmax, file 2 is tasmin")
  }
                            
  if(varname=="tasmin"){
    filesplit = do.call("c",strsplit(file1,"/",fixed=TRUE))
    filesplit2 = strsplit(filesplit[length(filesplit)],"_")[[1]]
    splitpart = strsplit(filesplit2[3],"-")[[1]]
    findthis = paste("tasmax_day_*-",splitpart[2],"-",splitpart[3],"_",paste(filesplit2[4:length(filesplit2)],collapse="_"),sep="")
    file2 = system(paste("ls /",filesplit[2],"/",filesplit[3],"/",filesplit[4],"/tasmax/",filesplit[6],"/",findthis,sep=""),intern=TRUE)
    message("file 1 is tasmin, file 2 is tasmax, switching the order")
    tmpfile1 = file2
    tmpfile2 = file1
    file1 = tmpfile1
    file2 = tmpfile2
    message("Switch complete")
  }
  
  if(varname=="pr"){
    message("pr only")
    #file1 = "/data2/3to5/pr/DeltaSD/pr_day_I35prp1-DeltaSD-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
    filesplit = do.call("c",strsplit(file1,"/",fixed=TRUE))
    filesplit2 = strsplit(filesplit[length(filesplit)],"_")[[1]]
    splitpart = strsplit(filesplit2[3],"-")[[1]]
    file2 = file1
    message("only one file used with pr")
  }
  
  if(combofunction=="heatwaves"){
    #splitname = do.call("c",strsplit(nameend,"_",fixed=TRUE))
    #nameuse = paste(substr(splitname[3],9,15),"0",substr(splitname[3],17,nchar(splitname[3])),"_historical_",splitname[5],"_",splitname[6],"_19810101-20051231",sep="")
    #q95list = system("ls /data2/3to5/I35/q95/*",intern=TRUE)
    if(filesplit[6]!="PARM"){
      findthis2 = paste("*-",splitpart[2],"-",substr(splitpart[3],1,2),"*",substr(splitpart[3],4,4),"*.nc",sep="")
      q95list = system(paste("ls /data2/3to5/I35/q95/",findthis2,sep=""),intern=TRUE)
      file3 = q95list[grep("tasmax",q95list)] # tasmax historical q95
      file4 = q95list[grep("tasmin",q95list)] # tasmin historical q95
    } else {
      splitpart = strsplit(filesplit2[3],"-")[[1]]
      findthis2 = paste("*-",splitpart[2],"-",substr(splitpart[3],1,2),"*",substr(splitpart[3],4,4),"*.nc",sep="")
      q95list = system(paste("ls /data2/3to5/I35/q95/",findthis2,sep=""),intern=TRUE)
      file3 = q95list[grep("tasmax",q95list)] # tasmax historical q95
      
      splitpart[3] = paste(substr(splitpart[3],1,1),"d",substr(splitpart[3],2,nchar(splitpart[3])),sep="")
      findthis2 = paste("*-",splitpart[2],"-",substr(splitpart[3],1,3),"*",substr(splitpart[3],5,5),"*.nc",sep="")
      q95list = system(paste("ls /data2/3to5/I35/q95/",findthis2,sep=""),intern=TRUE)
      file4 = q95list[grep("tasmin",q95list)] # tasmin historical q95
    }
    
  }
    
  if(combofunction=="R99" | combofunction=="R95" | combofunction=="R90" | combofunction=="R99freq" | combofunction=="R95freq" | combofunction=="R90freq"){
    #splitname = do.call("c",strsplit(nameend,"_",fixed=TRUE))
    #nameuse = paste(substr(splitname[3],9,15),"0",substr(splitname[3],17,nchar(splitname[3])),"_historical_",splitname[5],"_",splitname[6],"_19810101-20051231",sep="")
    #q95list = system("ls /data2/3to5/I35/q95/*",intern=TRUE)
    quant = as.numeric(substr(combofunction,2,3))
    findthis2 = paste("*-",splitpart[2],"-",substr(splitpart[3],1,2),"*",substr(splitpart[3],4,4),"*q",quant,".nc",sep="")
    qlist = system(paste("ls /data2/3to5/I35/prq/",findthis2,sep=""),intern=TRUE)
    file3 = qlist
  }
                          
  years = yearlydataperiod[1]:yearlydataperiod[2]
                            
  dates = datesdataperiod
                            
    if(combofunction=="heatwaves" | combofunction=="growing_season_length" | combofunction=="freeze_thaw_cycles"){              
      for(y in 1:length(years)){
        message("Combo function is: ",combofunction)
        yearidx = which(as.numeric(substr(dates,1,4))==years[y])
        
        test = nc_open(file1) # should always be tasmax file
        tempdata1 = ncvar_get(test,"tasmax",start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
                            
        if(y==1){
          yearlydat = array(NA,dim=c(dim(tempdata1)[1],dim(tempdata1)[2],length(years)))
          lat = ncvar_get(test,dimnames[2])
          lon = ncvar_get(test,dimnames[1])
          times = ncvar_get(test,dimnames[3])
          domainmask = ifelse(is.na(tempdata1[,,1])==FALSE,1,0)
          startdate = substr(test$dim[[dimnames[3]]]$units,12,21)
        }
        nc_close(test)
                            
        test = nc_open(file2)
        tempdata2 = ncvar_get(test,"tasmin",start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
        nc_close(test)
                            
        if(combofunction=="growing_season_length") {
          for(r in 1:length(lon)){
            for(c in 1:length(lat)){
              yearlydat[r,c,y] = GSLcalc(tempdata1[r,c,],tempdata2[r,c,],inputtimes=times[yearidx],startdate=startdate)
              #message("Finished Calculation R ",r," and C ",c)
            }
          }
          yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
        }
        
        if(combofunction=="freeze_thaw_cycles") {
          tmpdat = ifelse(tempdata1>273.15 & tempdata2<= 272.15,1,0)
          yearlydat[,,y] = apply(tmpdat,c(1,2),sum,na.rm=TRUE)
          yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
          message("min value = ",min(yearlydat[,,y],na.rm=TRUE)," max value = ",max(yearlydat[,,y],na.rm=TRUE))
        }
        
        if(combofunction=="heatwaves"){
          test = nc_open(file3)
          q95tmax = ncvar_get(test,"tmaxq95",start=c(1,1),count=c(-1,-1))
          nc_close(test)
          test = nc_open(file4)
          q95tmin = ncvar_get(test,"tminq95",start=c(1,1),count=c(-1,-1))
        nc_close(test)
                            
          for(r in 1:length(lon)){
            for(c in 1:length(lat)){
              yearlydat[r,c,y] = heatwave.calc(tempdata1[r,c,],tempdata2[r,c,],q95tmax[r,c],q95tmin[r,c])
            #message("Calcs complete for r ",r," and c ",c)
            }
          }
          yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
        }
        rm(tempdata1)
        rm(tempdata2)
      }
    }
  
    if(combofunction=="R99" | combofunction=="R95" | combofunction=="R90" | combofunction=="R99freq" | combofunction=="R95freq" | combofunction=="R90freq"){
      for(y in 1:length(years)){
        yearidx = which(as.numeric(substr(dates,1,4))==years[y])
        
        test = nc_open(file1) # should always be a pr file
        tempdata1 = ncvar_get(test,"pr",start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))*86400
      
        if(y==1){
          yearlydat = array(NA,dim=c(dim(tempdata1)[1],dim(tempdata1)[2],length(years)))
          lat = ncvar_get(test,dimnames[2])
          lon = ncvar_get(test,dimnames[1])
          times = ncvar_get(test,dimnames[3])
          domainmask = ifelse(is.na(tempdata1[,,1])==FALSE,1,0)
          startdate = substr(test$dim[[dimnames[3]]]$units,12,21)
        }
        nc_close(test)
      
        test = nc_open(file3)
        RQvalue = ncvar_get(test,paste("prq",quant,sep=""))
        nc_close(test)
      
        for(r in 1:length(lon)){
          for(c in 1:length(lat)){
            if(all(is.na(tempdata1[r,c,])==TRUE)==FALSE){
              if(combofunction=="R99" | combofunction=="R95" | combofunction=="R90"){
                yearlydat[r,c,y] = quantsum(tempdata1[r,c,],qval=RQvalue[r,c])  
              }
              if(combofunction=="R99freq" | combofunction=="R95freq" | combofunction=="R90freq"){
                yearlydat[r,c,y] = quantfreqsum(tempdata1[r,c,],qval=RQvalue[r,c])
              }
            }
            #message("Finished Calculation R ",r," and C ",c)
          }
        }
        yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
        message("Finished Calcs for year ",years[y])
        rm(tempdata1)
        }
    }
list(lon,lat,yearlydat)
}
                            

####

####

netcdftopointcomboclimo = function(filenames,varnames=c(),dimnames=c("lon","lat","time"),locstart =c(1,1),loccount=c(-1,-1),yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),combofunction="gsl"){
  
  # filenames = 2-4 unique files each with a different variable needed, signified by varnames, separated by commas
  # varnames should be the unique netcdf variable names needed for creating compound variables
  
  #filenames=filenames;varnames=varnamesin;dimnames=c("lon","lat","time");yearlydataperiod=c(1981,2005);datesdataperiod=datesin;combofunction=appfunc;locstart=locstart;loccount=loccount
  
  files = unlist(strsplit(filenames,",",fixed=TRUE))
  
  years = yearlydataperiod[1]:yearlydataperiod[2]
  
  dates = datesdataperiod
  
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(dates,1,4))==years[y])
    
    test = nc_open(files[1])
    tempdata1 = ncvar_get(test,varnames[1],start=c(locstart[1],locstart[2],yearidx[1]),count=c(loccount[1],loccount[2],length(yearidx)))
    if(varnames[1]=="pr") tempdata1=tempdata1*86400
    
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata1)[1],dim(tempdata1)[2],length(years)))
      times = ncvar_get(test,dimnames[3])
      startdate = substr(test$dim[[4]]$units,12,21)
    }
    nc_close(test)
    
    test = nc_open(files[2])
    tempdata2 = ncvar_get(test,varnames[2],start=c(locstart[1],locstart[2],yearidx[1]),count=c(loccount[1],loccount[2],length(yearidx)))
    if(varnames[2]=="pr") tempdata2=tempdata2*86400
    nc_close(test)
    
    if(combofunction=="growing_season_length") {
      for(r in 1:loccount[1]){
        for(c in 1:loccount[2]){
          yearlydat[r,c,y] = GSLcalc(tempdata1[r,c,],tempdata2[r,c,],inputtimes=times[yearidx],startdate=startdate)
          message("Finished Calculation R ",r," and C ",c)
        }
      }
    }
    if(combofunction=="heatwaves"){
      test = nc_open(files[3])
      q95tmax = ncvar_get(test,varnames[3],start=c(locstart[1],locstart[2]),count=c(loccount[1],loccount[2]))
      nc_close(test)
      test = nc_open(files[4])
      q95tmin = ncvar_get(test,varnames[4],start=c(locstart[1],locstart[2]),count=c(loccount[1],loccount[2]))
      nc_close(test)
      
      for(r in 1:loccount[1]){
        for(c in 1:loccount[2]){
          yearlydat[r,c,y] = heatwave.calc(tempdata1[r,c,],tempdata2[r,c,],q95tmax[r,c],q95tmin[r,c])
          message("Calcs complete for r ",r," and c ",c)
        }
      }
    }
    rm(tempdata1)
    rm(tempdata2)
    message("Finished Calcs for year ",years[y])
  }
  climocalc = apply(yearlydat,c(1,2),mean,na.rm=TRUE)
  climolocvalue = mean(climocalc,na.rm=TRUE)
  climolocvalue
}

####

netcdftopointtoclimocalcs = function(filename,varname,dimnames=c("lon","lat","time"),locstart =c(1,1),loccount=c(-1,-1),threscalc=FALSE,thres=NA,condition="gte",yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),appliedfunction="mean"){
  
  #filename=histfilelist[7];varname=varin;dimnames=c("lon","lat","time");locstart=locstart;loccount=loccount;threscalc=TC;thres=TH;condition="gte";yearlydataperiod=c(1981,2005);datesdataperiod=datesin;appliedfunction=appfunc
  
  years = yearlydataperiod[1]:yearlydataperiod[2]
  
  dates = datesdataperiod
  
  test = nc_open(filename)
  inputdata = ncvar_get(test,varname,start=c(locstart[1],locstart[2],1),count=c(loccount[1],loccount[2],-1))
  if(varname=="pr") inputdata=inputdata*86400
    times = ncvar_get(test,dimnames[3])
    domainmask = ifelse(is.na(inputdata[,,1])==FALSE,1,0)
    startdate = substr(test$dim[[4]]$units,12,21)
  nc_close(test)
  
  yearlydat = array(NA,dim=c(loccount[1],loccount[2],length(years)))
  
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(dates,1,4))==years[y])
    tempdata=inputdata[,,yearidx]
    
    if(threscalc==TRUE){
      dmask = array(ifelse(is.na(tempdata[,,1])==FALSE,1,0),dim=dim(tempdata))
      if(condition=="gte") tempdata=ifelse(tempdata>=thres,1,0)
      if(condition=="gt") tempdata=ifelse(tempdata>thres,1,0)
      if(condition=="lte") tempdata=ifelse(tempdata<=thres,1,0)
      if(condition=="lt") tempdata=ifelse(tempdata<thres,1,0)
      tempdata = ifelse(dmask==1,tempdata,NA)
    }
    
    if(appliedfunction=="mean") yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
    if(appliedfunction=="prcptot") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="prmean") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="prmed") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),median,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="sum") {
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="max"){
      yearlydat[,,y] = apply(tempdata,c(1,2),max,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="rx5day"){
      yearlydat[,,y] = apply(tempdata,c(1,2),calcrollsum,size=5)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="lastfreeze"){
      #test1 = apply(tempdata,c(1,2),lastfreeze,startdate=startdate,inputtimes=times)
      ptm = proc.time()
      for(r in 1:length(loccount[1])){
        for(c in 1:length(loccount[2])){
          message("Working on r ",r," and c ",c)
          yearlydat[r,c,y] = lastfreeze(tempdata[r,c,],startdate=startdate,inputtimes=times)
        }
      }
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
      ptm.end = proc.time()
    }
    if(appliedfunction=="maxdryspell"){
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=1,outtype="max") # changed from trace to 1mm 11/4
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="maxwetspell"){
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=1,outtype="max") # changed from trace to 1mm 11/4
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    
    rm(tempdata)
    message("Finished Calcs for year ",years[y])
  }
  climocalc = apply(yearlydat,c(1,2),mean,na.rm=TRUE)
  climolocvalue = mean(climocalc,na.rm=TRUE)
  climolocvalue
}


####

statemask = function(inputdata,inputlon,inputlat,state="oklahoma") {

  outline = data.frame(map("state",regions=state,plot=FALSE)[c("x","y")])
  rows = 1:length(inputlon)
  cols = 1:length(inputlat)
  LON = rep(inputlon,each=length(inputlat))
  LAT = rep(inputlat,length(inputlon))
  R = rep(rows,each=length(inputlat))
  C = rep(cols,length(inputlon))
  modgrid = data.frame(R,C,LON,LAT)
  library(sp)
  modin = point.in.polygon(modgrid$LON,modgrid$LAT,outline$x,outline$y)
  modgrid = cbind(modgrid,modin)
  mask = array(NA,dim=dim(inputdata))
  
  lonrange = range(modgrid$R[which(modgrid$modin>=1)]) + c(-2,2)
  latrange = range(modgrid$C[which(modgrid$modin>=1)]) + c(-2,2)
  
  if(any(lonrange<0)==TRUE) lonrange[which(lonrange<0)]=0
  if(any(latrange<0)==TRUE) latrange[which(latrange<0)]=0
  
  if(any(lonrange>length(inputlon))==TRUE) lonrange[which(lonrange>length(inputlon))]=length(inputlon)
  if(any(latrange>length(inputlat))==TRUE) latrange[which(latrange>length(inputlat))]=length(inputlat)
  
  outputlon = inputlon[lonrange[1]:lonrange[2]]
  outputlat = inputlat[latrange[1]:latrange[2]]
  
  for(r in 1:length(lon)){
    for(c in 1:length(lat)){
      mask[r,c,] = modgrid$modin[which(modgrid$R==r & modgrid$C==c)]
    }
  }
  outputdata = ifelse(mask>=1,inputdata,NA)
  outputdata = outputdata[lonrange[1]:lonrange[2],latrange[1]:latrange[2],]
  output = list(lon=outputlon,lat=outputlat,outputdata=outputdata)
  output
}
####

netcdfdatagrab = function(filename,varname,dimnames=c("lon","lat","time"),threscalc=FALSE,thres=NA,condition="gte",subsettime=FALSE,subsetpoints=c(1,-1)){
  require(ncdf4)
  test = nc_open(filename)
  if(subsettime==TRUE){
    tempdata = ncvar_get(test,varname,start=c(1,1,subsetpoints[1]),count=c(-1,-1,subsetpoints[2]))
  } else {
    tempdata = ncvar_get(test,varname)
  }
  lat = ncvar_get(test,dimnames[2])
  lon = ncvar_get(test,dimnames[1])
  times = ncvar_get(test,dimnames[3])
  domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
  nc_close(test)
  
  if(threscalc==TRUE){
    dmask = array(ifelse(is.na(tempdata[,,1])==FALSE,1,0),dim=dim(tempdata))
    if(condition=="gte") tempdata=ifelse(tempdata>=thres,1,0)
    if(condition=="gt") tempdata=ifelse(tempdata>thres,1,0)
    if(condition=="lte") tempdata=ifelse(tempdata<=thres,1,0)
    if(condition=="lt") tempdata=ifelse(tempdata<thres,1,0)
    tempdata = ifelse(dmask==1,tempdata,NA)
  }
  
  output = list(varname,lon,lat,times,domainmask,tempdata)
  output
}

###

netcdftoseasonalcalcs = function(filename,varname,dimnames=c("lon","lat","time"),season="JJA",threscalc=FALSE,thres=NA,condition="gte",yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),appliedfunction="mean"){
  
  years = yearlydataperiod[1]:yearlydataperiod[2]
  
  dates = datesdataperiod
  
  for(y in 1:length(years)){
    
    if(season=="DJF" | season=="MAM" | season=="JJA" | season=="SON" | season=="mon" | season=="snowmelt"){
      if(season!="DJF") {
        if(season=="mon") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=7 & as.numeric(substr(dates,6,7))<=10)
        if(season=="snowmelt") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=2 & as.numeric(substr(dates,6,7))<=6)
        if(season=="JJA") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=6 & as.numeric(substr(dates,6,7))<=8)
        if(season=="SON") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=9 & as.numeric(substr(dates,6,7))<=11)
        if(season=="MAM") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=3 & as.numeric(substr(dates,6,7))<=5)
      } else {
        yearidx = which(substr(dates,1,4)==years[y] & (as.numeric(substr(dates,6,7))==12 | as.numeric(substr(dates,6,7))<=2))
      }
    } else {
      if(season=="Jan") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==1)
      if(season=="Feb") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==2)
      if(season=="Mar") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==3)
      if(season=="Apr") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==4)
      if(season=="May") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==5)
      if(season=="Jun") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==6)
      if(season=="Jul") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==7)
      if(season=="Aug") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==8)
      if(season=="Sep") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==9)
      if(season=="Oct") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==10)
      if(season=="Nov") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==11)
      if(season=="Dec") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))==12)
    }
    
    test = nc_open(filename)
    tempdata = ncvar_get(test,varname,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    if(varname=="pr") tempdata=tempdata*86400
    
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata)[1],dim(tempdata)[2],length(years)))
      lat = ncvar_get(test,dimnames[2])
      lon = ncvar_get(test,dimnames[1])
      times = ncvar_get(test,dimnames[3])
      domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
    }
    nc_close(test)
    
    if(threscalc==TRUE){
      dmask = array(ifelse(is.na(tempdata[,,1])==FALSE,1,0),dim=dim(tempdata))
      if(condition=="gte") tempdata=ifelse(tempdata>=thres,1,0)
      if(condition=="gt") tempdata=ifelse(tempdata>thres,1,0)
      if(condition=="lte") tempdata=ifelse(tempdata<=thres,1,0)
      if(condition=="lt") tempdata=ifelse(tempdata<thres,1,0)
      tempdata = ifelse(dmask==1,tempdata,NA)
    }
    
    if(appliedfunction=="mean") yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
    if(appliedfunction=="prcptot") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="prmean") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="prmed") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),median,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="sum") {
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="max"){
      yearlydat[,,y] = apply(tempdata,c(1,2),max,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    rm(tempdata)
  }
  
  list(lon,lat,yearlydat)
}
###

netcdftoyearlycombocalcs = function(filenames,varnames=c(),dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),combofunction="gsl"){
  
  # filenames = 2-4 unique files each with a different variable needed, signified by varnames, separated by commas
  # varnames should be the unique netcdf variable names needed for creating compound variables
  
  #filenames=filenames;varnames=varnamesin;dimnames=c("lon","lat","time");yearlydataperiod=c(1981,2005);datesdataperiod=datesin;combofunction=appfunc
  
  files = unlist(strsplit(filenames,",",fixed=TRUE))
  
  years = yearlydataperiod[1]:yearlydataperiod[2]
  
  dates = datesdataperiod
  
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(dates,1,4))==years[y])
    
    test = nc_open(files[1])
    tempdata1 = ncvar_get(test,varnames[1],start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    if(varnames[1]=="pr") tempdata1=tempdata1*86400
    
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata1)[1],dim(tempdata1)[2],length(years)))
      lat = ncvar_get(test,dimnames[2])
      lon = ncvar_get(test,dimnames[1])
      times = ncvar_get(test,dimnames[3])
      domainmask = ifelse(is.na(tempdata1[,,1])==FALSE,1,0)
      startdate = substr(test$dim[[dimnames[3]]]$units,12,21)
    }
    nc_close(test)
    
    test = nc_open(files[2])
    tempdata2 = ncvar_get(test,varnames[2],start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    if(varnames[2]=="pr") tempdata2=tempdata2*86400
    nc_close(test)
    
    if(combofunction=="growing_season_length") {
      for(r in 1:length(lon)){
        for(c in 1:length(lat)){
          yearlydat[r,c,y] = GSLcalc(tempdata1[r,c,],tempdata2[r,c,],inputtimes=times[yearidx],startdate=startdate)
          message("Finished Calculation R ",r," and C ",c)
          }
      }
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(combofunction=="heatwaves"){
      test = nc_open(files[3])
      q95tmax = ncvar_get(test,varnames[3],start=c(1,1),count=c(-1,-1))
      nc_close(test)
      test = nc_open(files[4])
      q95tmin = ncvar_get(test,varnames[4],start=c(1,1),count=c(-1,-1))
      nc_close(test)

        for(r in 1:length(lon)){
          for(c in 1:length(lat)){
            yearlydat[r,c,y] = heatwave.calc(tempdata1[r,c,],tempdata2[r,c,],q95tmax[r,c],q95tmin[r,c])
          message("Calcs complete for r ",r," and c ",c)
            }
        }
        yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    rm(tempdata1)
    rm(tempdata2)
    message("Finished Calcs for year ",years[y])
  }
  
  list(lon,lat,yearlydat)
}


###

netcdffreqcalcs = function(filename,varname,dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),appliedfunction="freq2011heat"){
  
  years = yearlydataperiod[1]:yearlydataperiod[2]
  
  dates = datesdataperiod
  
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(dates,1,4))==years[y])
    
    test = nc_open(filename)
    tempdata = ncvar_get(test,varname,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    if(varname=="pr"){
      tempdata=tempdata*86400
    } 
    
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata)[1],dim(tempdata)[2],length(years)))
      lat = ncvar_get(test,dimnames[2])
      lon = ncvar_get(test,dimnames[1])
      times = ncvar_get(test,dimnames[3])
      domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
      startdate = substr(test$dim[[dimnames[3]]]$units,12,21)
    }
    nc_close(test)
    
    dmask = array(ifelse(is.na(tempdata[,,1])==FALSE,1,0),dim=dim(tempdata))
    if(appfunc=="freq2011heat"){
      tempdata=ifelse(tempdata >= 310.928,1,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
      
      if(y==1){
        test = nc_open("/data2/3to5/I35/METDATA_tmax100.nc")
        thresvals = ncvar_get(test,"tmax100_2011")
        nc_close(test)
      }
      yearlydat[,,y]=ifelse(yearlydat[,,y]>thresvals,1,0)
    } 
    
    rm(tempdata)
    message("Finished Calcs for year ",years[y])
  }
  
  list(lon,lat,yearlydat)
}

###

netcdfdroughtcalcs = function(filename,varname,yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),appliedfunction="SPI3"){
  
  #filename = "/data2/3to5/I35/pr/DeltaSD/pr_day_I35prp1-DeltaSD-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
  #yearlydataperiod=c(1981,2005)
  #datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
  #appliedfunction="SPI3"
  periodin = as.numeric(substr(appliedfunction,4,4))
  years = yearlydataperiod[1]:yearlydataperiod[2]
  
  SPIdat = spi.calc(filename,period=periodin,dates =datesdataperiod)
  lon = SPIdat[[1]]
  lat = SPIdat[[2]]
  months = SPIdat[[3]]
  SPIout = SPIdat[[4]]
  
  #testsfc = list(x=lon,y=lat,z=SPIout[,,48])
  #surface(testsfc,type="I")
  
  yearlydat = array(NA,dim=c(length(lon),length(lat),length(years)))
  
  for(y in 1:length(years)){
    idx = which(as.numeric(substr(months,1,4))==years[y])
    yearlydat[,,y] = apply(SPIout[,,idx],c(1,2),mean,na.rm=TRUE)
    yearlydat[,,y] = ifelse(is.na(SPIout[,,48])==FALSE,yearlydat[,,y],NA)
    message("Finished SPI calcs for year ",years[y])
  }
  
  #testsfc = list(x=lon,y=lat,z=yearlydat[,,4])
  #surface(testsfc,type="I")
  
  list(lon,lat,yearlydat)
}



###

netcdftoyearlycalcs = function(filename,varname,dimnames=c("lon","lat","time"),threscalc=FALSE,thres=NA,condition="gte",yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),appliedfunction="mean"){
  
  #filename=histfilelist[i];varname=varin;dimnames=c("lon","lat","time");threscalc=TC;thres=TH;condition=cond;yearlydataperiod=c(1981,2005);datesdataperiod=datesin;appliedfunction=appfunc
  
  years = yearlydataperiod[1]:yearlydataperiod[2]
  #print(years)
  dates = datesdataperiod
  #print(dates[1:10])
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(dates,1,4))==years[y])
    #print(yearidx)
    #print(length(yearidx))
    test = nc_open(filename)
      tempdata = ncvar_get(test,varname,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
      if(varname=="pr"){
        tempdata=tempdata*86400
        #message("Did pr unit conversion from mks to mm")  
        #message("Min value = ",min(tempdata,na.rm=TRUE)," Max value = ",max(tempdata,na.rm=TRUE))
      } 
      
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata)[1],dim(tempdata)[2],length(years)))
      lat = ncvar_get(test,dimnames[2])
      lon = ncvar_get(test,dimnames[1])
      times = ncvar_get(test,dimnames[3])
      domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
      startdate = substr(test$dim[[dimnames[3]]]$units,12,21)
    }
    nc_close(test)
    
    if(threscalc==TRUE){
      message("Doing threshold calculations")
      message("The threshold being used is ",thres)
      message("The condition is ",condition)
      
      dmask = array(ifelse(is.na(tempdata[,,1])==FALSE,1,0),dim=dim(tempdata))
      if(condition=="gte"){
        tempdata=ifelse(tempdata >= as.numeric(thres),1,0)
      } 
      if(condition=="gt"){
        tempdata=ifelse(tempdata > as.numeric(thres),1,0)
      } 
      if(condition=="lte"){
        tempdata=ifelse(tempdata <= as.numeric(thres),1,0)
      } 
      if(condition=="lt"){
        tempdata=ifelse(tempdata < as.numeric(thres),1,0)
      } 
      tempdata = ifelse(dmask==1,tempdata,NA)
    }
    
    if(appliedfunction=="mean") yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
    if(appliedfunction=="prmean") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),mean,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="prmed") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),median,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="prcptot") {
      tempdata=ifelse(tempdata>=0.254,tempdata,0)
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="sum") {
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
      message("Min value = ",min(yearlydat[,,y],na.rm=TRUE)," Max value = ",max(yearlydat[,,y],na.rm=TRUE))
    }
    if(appliedfunction=="max"){
      yearlydat[,,y] = apply(tempdata,c(1,2),max,na.rm=TRUE)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="rx5day"){
      yearlydat[,,y] = apply(tempdata,c(1,2),calcrollsum,size=5)
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="gsl2"){
      ptm = proc.time()
      for(r in 1:length(lon)){
        for(c in 1:length(lat)){
          message("Working on r ",r," and c ",c)
          yearlydat[r,c,y] = GSLcalc2(tempdata[r,c,],inputtimes=times[yearidx],startdate=startdate)
        }
      }
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
      ptm.end = proc.time()
    }
    if(appliedfunction=="lastfreeze"){
      #test1 = apply(tempdata,c(1,2),lastfreeze,startdate=startdate,inputtimes=times)
      ptm = proc.time()
      for(r in 1:length(lon)){
        for(c in 1:length(lat)){
          message("Working on r ",r," and c ",c)
          yearlydat[r,c,y] = lastfreeze(tempdata[r,c,],startdate=startdate,inputtimes=times[yearidx])
        }
      }
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
      ptm.end = proc.time()
    }
    if(appliedfunction=="maxdryspell"){
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=1,outtype="max") # changed to 1 mm 11/4
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="maxwetspell"){
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=1,outtype="max") # changed threshold to 1 mm 11/4
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    
    rm(tempdata)
    message("Finished Calcs for year ",years[y])
  }
  
  list(lon,lat,yearlydat)
}

###

yearlycalcs = function(dailydata,yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day"),appliedfunction="mean"){
  years = yearlydataperiod[1]:yearlydataperiod[2]
  yearlydat = array(NA,dim=c(dim(dailydata)[1],dim(dailydata)[2],length(years)))
  for(y in 1:length(years)){
    yearidx = which(substr(datesdataperiod,1,4)==years[y])
    tempdat = dailydata[,,yearidx]
    if(appliedfunction=="mean") totalvals = apply(tempdat,c(1,2),mean,na.rm=TRUE)
    if(appliedfunction=="sum") totalvals = apply(tempdat,c(1,2),sum,na.rm=TRUE)
    yearlydat[,,y]=totalvals
    rm(tempdat)
  }
  yearlydat
}

####

diffcalc = function(projclimo,histclimo,type="absolute"){
  if(type=="absolute") output = projclimo-histclimo
  if(type=="percent") output = ((projclimo-histclimo)/histclimo)*100
  if(type=="relative") output = (projclimo/histclimo)
  output
}


####

climocalc = function(yearlydata,yearlydataperiod=c(1981,2005),climoperiod=c(1981,2005)){
  # yearlydata = array with means or sums by year
  # yearlydataperiod = period of the yearly output
  # climoperiod = period for the climotology calculation
  years = yearlydataperiod[1]:yearlydataperiod[2]
  yearidx = which(years<=climoperiod[2] & years>=climoperiod[1])
  climo = apply(yearlydata[,,yearidx],c(1,2),mean,na.rm=TRUE)
  message("Climo min value = ",min(climo,na.rm=TRUE)," max value = ",max(climo,na.rm=TRUE))
  climo
}


#####

obs_colorramp = function(obs,colorchoice,Blimit){
  
  diffrange = range(obs,na.rm=TRUE)
  diffrange[1]=floor(diffrange[1])
  diffrange[2]=ceiling(diffrange[2])
  
  breakcheck = 1
  breaklist = c(0.01,0.02,0.025,0.05,0.1,0.2,0.25,0.3,0.5,1,2,3,4,5,10,20,25,30,50)
  
  actualbins = diff(diffrange)/breaklist
  actidx = which(actualbins<Blimit)
  diffact = actualbins[actidx]-floor(actualbins[actidx])
  
  if(any(diffact==0)==TRUE){
    diffidx = which(diffact==0)
    breakcheck=actidx[diffidx[1]]
  } else {
    diffrange[1] = floor(diffrange[1]/10)*10
    diffrange[2] = ceiling(diffrange[2]/10)*10
    actualbins = diff(diffrange)/breaklist
    actidx = which(actualbins<Blimit)
    diffact = actualbins[actidx]-floor(actualbins[actidx])
    diffidx = which(diffact==0)
    breakcheck=actidx[diffidx[1]]
  }
  
  zlimdiff = diffrange
  breaksdiff = c(seq(diffrange[1],diffrange[2],by=breaklist[breakcheck]))
  
  if(colorchoice == "redtowhite") colorbardiff = colorRampPalette(c("red","grey96"))(length(breaksdiff)-1)
  if(colorchoice == "whitetored") colorbardiff = colorRampPalette(c("grey96","red"))(length(breaksdiff)-1)
  if(colorchoice == "whitetogreen") colorbardiff = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksdiff)-1)

  output = list(zlimdiff,breaksdiff,colorbardiff)
  output
  
}

####

colorramp = function(inputdata,colorchoice,Blimit,type = "difference",use_fixed_scale = FALSE, fixed_scale=c(-100,100)){
  
  #inputdata=c(subseterr,tmperr)
  #colorchoice="browntoblue"
  #Blimit = 20
  #type="difference"
  #use_fixed_scale = FALSE
  #fixed_scale = c(0,1)
  
  if(use_fixed_scale==FALSE){
    message("Not using a fixed scale")
    datarange = range(inputdata,na.rm=TRUE)
    if(all(abs(datarange)<1)==TRUE){
      i=1
      while(all(abs(datarange)<1)==TRUE){
        datarange=datarange*(10*i)
        
        if(all(abs(datarange)<1)==TRUE) i=i+1
        if(all(abs(datarange)<1)==FALSE) break;
      }
      datarange[1]=floor(datarange[1])/(10*i)
      datarange[2]=ceiling(datarange[2])/(10*i)
    } else {
      datarange[1]=floor(datarange[1])
      if(datarange[1] %% 2 != 0) datarange[1]=datarange[1]-1
      datarange[2]=ceiling(datarange[2])
      if(datarange[2] %% 2 != 0) datarange[2]=datarange[2]+1
      
      if(datarange[1]<0 & type=="difference"){
        if(abs(datarange[1])>10 & abs(datarange[1])<100){
          datarange[1]=-(ceiling(abs(datarange[1])/10)*10)
        }
        if(abs(datarange[1])>100 & abs(datarange[1])<1000){
          datarange[1]=-(ceiling(abs(datarange[1])/10)*10)
        }
        if(datarange[2]>1000 & datarange[2]<10000){
          datarange[2]=-(ceiling(abs(datarange[1])/100)*100)
        }
      }
      
      if(datarange[1]>0 & type=="difference"){
        datarange[1]==0
      }
      if(datarange[2]>0){
        if(datarange[2]>10 & datarange[2]<100){
          datarange[2]=ceiling(datarange[2]/10)*10
        }
        if(datarange[2]>100 & datarange[2]<1000){
         datarange[2]=ceiling(datarange[2]/10)*10
        }
        if(datarange[2]>1000 & datarange[2]<10000){
          datarange[2]=ceiling(datarange[2]/100)*100
        }
      }
      
    }
  } else {
    message("Using a fixed scale")
    #tmp = strsplit(fixed_scale,",")
    #datarange = c(as.numeric(tmp[[1]][1]),as.numeric(tmp[[1]][2]))
    datarange = fixed_scale
  }
  
  if(datarange[1]>=0 & type=="difference"){centerpoint = 0; startpoint=0; datarange[1]=0; message("type=difference")}
  if(datarange[1]<0 & type=="difference"){centerpoint = 0; startpoint=datarange[1]; message("type=difference"); if(datarange[2]<0){datarange[2]=0}}
  
  if(type=="ratio"){centerpoint=1; startpoint=datarange[1]; message("type=ratio")}
  if(type=="raw"){centerpoint=datarange[1]; startpoint=datarange[1]; message("type=raw")}
  
  breakcheck = 1
  breaklist = c(0.0001,0.0002,0.00025,0.0003,0.0005,0.001,0.002,0.0025,0.003,0.005,0.01,0.02,0.025,0.03,0.05,0.1,0.2,0.25,0.3,0.5,1,2,3,4,5,10,20,25,30,40,50,75,80,100,200,250,300,400,500,750,800,1000)
  
  actualbins = diff(datarange)/breaklist
  #print(actualbins)
  actidx = which(actualbins<=Blimit)
  #print(actidx)
  #dataact = actualbins[actidx]-floor(actualbins[actidx])
  dataact = actualbins[actidx]-round(actualbins[actidx],0)
  #print(dataact)
  
  if(any(abs(dataact)<=1E-6)==TRUE){
    message("exact match for bins")
    dataidx = which(dataact<=1E-6)
    breakcheck=actidx[dataidx[1]]
    #breakcheck=actidx[dataidx[length(dataidx)]]
  } else {
    message("no exact match going through while loop")
    checkpoint = any(dataact<=1E-6)
    counter=1
    while(checkpoint==FALSE){
      datarange[1] = floor(datarange[1]/(10^counter))*10^counter
      datarange[2] = ceiling(datarange[2]/(10^counter))*10^counter
      actualbins = diff(datarange)/breaklist
      actidx = which(actualbins<Blimit)
      dataact = actualbins[actidx]-floor(actualbins[actidx])
      dataidx = which(dataact<=1E-6)
      
      if(length(dataidx)>=1){
        breakcheck=actidx[dataidx[1]]
        checkpoint = any(dataact<=1E-6)
        break
      } else {
        counter=counter+1
        checkpoint = any(dataact<=1E-6)
      }
      
    }
  }
  
  if(datarange[1]==0 & colorchoice=="redtoblue") colorchoice="whitetoblue" # check this line
  if(datarange[2]==0 & colorchoice=="redtoblue") colorchoice="redtowhite" # check this line
  if(datarange[2]==0 & colorchoice=="greentobrown") colorchoice="greentowhite"
  if(datarange[2]==0 & colorchoice=="browntogreen") colorchoice="browntowhite"
  if(datarange[2]==0 & colorchoice=="orangetopurple") colorchoice="orangetowhite"
  if(datarange[2]==0 & colorchoice=="yellowtopurple") colorchoice="yellowtowhite"
  if(startpoint==0 & centerpoint==0 & colorchoice=="bluetored") colorchoice="whitetored"
  #if(startpoint==0 & centerpoint==0 & colorchoice=="browntogreen") colorchoice="whitetogreen"
  if(startpoint==0 & colorchoice=="purpletoorange") colorchoice="whitetoorange"
  if(startpoint==0 & colorchoice=="orangetopurple") colorchoice="whitetopurple"
  if(startpoint==0 & colorchoice=="browntoblue") colorchoice="whitetoblue2"
  
  zlimdiff = datarange
  breaksdiff = c(seq(datarange[1],datarange[2],by=breaklist[breakcheck]))
  
  if(any(breaksdiff==centerpoint)==FALSE & zlimdiff[1]<centerpoint & zlimdiff[2]>centerpoint){
  idx = which(abs(breaksdiff)==min(abs(breaksdiff))) # changed from min to max
  if(length(idx)==1) breaksdiff[idx]=centerpoint
  if(length(idx)>1){
  breaksdiff = c(breaksdiff[1:(idx[1]-1)],centerpoint,breaksdiff[(idx[2]+1):length(breaksdiff)])
  }
  } 
  
  if(max(breaksdiff,na.rm=TRUE)<max(datarange,na.rm=TRUE)){
    breaksdiff = c(breaksdiff,(breaksdiff[length(breaksdiff)]+diff(breaksdiff)[1]))
    datarange[2]=max(breaksdiff,na.rm=TRUE)
  }
  
  if(min(breaksdiff,na.rm=TRUE)>min(datarange,na.rm=TRUE)){
    breaksdiff = c((breaksdiff[1]-diff(breaksdiff)[1]),breaksdiff)
    datarange[1]=min(breaksdiff,na.rm=TRUE)
  }
  
  message("zlimdiff: ",zlimdiff)
  message("centerpoint: ",centerpoint)
  message("startpoint: ",startpoint)
  message("colorchoice: ",colorchoice)
  
  if(startpoint==centerpoint){
  message("startpoint matches centerpoint")
  if(colorchoice=="rainbow") colorbardiff = colorRampPalette(c("#2d004b","#4393c3","#35978f","#ffffb2","#e6550d","#b2182b","#67001f"))(length(breaksdiff)-1)
  
  if(colorchoice == "whitetored") colorbardiff = colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-1)
  if(colorchoice == "redtowhite") colorbardiff = colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(length(breaksdiff)-1)
  if(colorchoice == "yellowtored") colorbardiff = colorRampPalette(c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"))(length(breaksdiff)-1)
  
  if(colorchoice == "whitetogreen") colorbardiff = colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-1)
  if(colorchoice == "greentowhite") colorbardiff = colorRampPalette(c("#003c30","#01665e","#35978f","#80cdc1","#c7eae5","#f5f5f5"))(length(breaksdiff)-1)
  
  if(colorchoice == "whitetobrown") colorbardiff = colorRampPalette(c("#f5f5f5","#f6e8c3","#dfc27d","#bf812d","#8c510a","#543005"))(length(breaksdiff)-1)
  
  if(colorchoice == "bluetored") colorbardiff = colorRampPalette(c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-1)
  
  if(colorchoice == "browntogreen") colorbardiff = colorRampPalette(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-1)
  if(colorchoice == "greentobrown") colorbardiff = colorRampPalette(c("#003c30","#01665e","#35978f","#80cdc1","#c7eae5","#f5f5f5","#f6e8c3","#dfc27d","#bf812d","#8c510a","#543005"))(length(breaksdiff)-1)
  
  if(colorchoice == "whitetoorange") colorbardiff = colorRampPalette(c("#f5f5f5","#feedde","#fdbe85","#fd8d3c","#e6550d","#a63603"))(length(breaksdiff)-1)
  if(colorchoice == "whitetopurple") colorbardiff = colorRampPalette(c("#f5f5f5","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(breaksdiff)-1)  
  if(colorchoice == "whitetogreenblue") colorbardiff = colorRampPalette(c("#f5f5f5","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))(length(breaksdiff)-1)
  
  if(colorchoice == "whitetoblue") colorbardiff = colorRampPalette(c("#f5f5f5","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(length(breaksdiff)-1)
  if(colorchoice == "whitetoblue2") colorbardiff = colorRampPalette(c("#f5f5f5","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))(length(breaksdiff)-1)
  if(colorchoice == "bluetowhite") colorbardiff = colorRampPalette(c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f5f5f5"))(length(breaksdiff)-1)
  
  if(colorchoice == "grayscale") colorbardiff = colorRampPalette(c("#ffffff","#d9d9d9","#bdbdbd","969696","#636363","#252525"))(length(breaksdiff)-1)
  
  } else {
    if(datarange[2]>centerpoint){
      message("datarange[2] > centerpoint")
      message("colorchoice = ",colorchoice)
      zeroidx = which(breaksdiff==centerpoint)
      if(colorchoice == "bluetored"){
        colorbardiff = c(colorRampPalette(c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "redtoblue"){
        colorbardiff = c(colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "browntogreen"){
        colorbardiff = c(colorRampPalette(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "greentobrown"){
        colorbardiff = c(colorRampPalette(c("#003c30","#01665e","#35978f","#80cdc1","#c7eae5","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#f6e8c3","#dfc27d","#bf812d","#8c510a","#543005"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "purpletoorange"){
        colorbardiff = c(colorRampPalette(c("#2d004b","#542788","#8073ac","#b2abd2","#d8daeb","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#feedde","#fdbe85","#fd8d3c","#e6550d","#a63603"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "orangetopurple"){
        colorbardiff = c(colorRampPalette(c("#a63603","#e6550d","#fd8d3c","#fdbe85","#feedde","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(breaksdiff)-zeroidx))
      } 
      
      if(colorchoice == "browntoblue"){
        message("making brown to blue color ramp")
        #colorbardiff = c(colorRampPalette(c("lightgoldenrod4","lightgoldenrod3","lightgoldenrod2","lightgoldenrod1","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))(length(breaksdiff)-zeroidx))
        colorbardiff = c(colorRampPalette(c("burlywood4","burlywood3","burlywood2","burlywood1","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))(length(breaksdiff)-zeroidx))
        #colorbardiff = c(colorRampPalette(c("tan4","tan3","tan2","tan1","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))(length(breaksdiff)-zeroidx))
        #colorbardiff = c(colorRampPalette(c("wheat4","wheat3","wheat2","wheat1","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))(length(breaksdiff)-zeroidx))
      } 
      
      if(colorchoice == "yellowtopurple"){
        colorbardiff = c(colorRampPalette(c("lightgoldenrod4","lightgoldenrod3","lightgoldenrod2","lightgoldenrod1","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(breaksdiff)-zeroidx))
      } 
      
      
    } else {
      message("odd if")
      if(colorchoice == "redtowhite"){
        colorbardiff = colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(length(breaksdiff)-1)
      } 
      if(colorchoice=="yellowtowhite"){
        colorbardiff = colorRampPalette(c("lightgoldenrod4","lightgoldenrod3","lightgoldenrod2","lightgoldenrod1","#f5f5f5"))(length(breaksdiff)-1)
      }
      
      if(colorchoice=="orangetowhite"){
        colorbardiff = colorRampPalette(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f5f5f5"))(length(breaksdiff)-1)
      }
      if(colorchoice == "greentowhite"){
        colorbardiff = colorRampPalette(c("#003c30","#01665e","#35978f","#80cdc1","#c7eae5","#f5f5f5"))(length(breaksdiff)-1)
      }
      if(colorchoice == "browntowhite"){
        colorbardiff = colorRampPalette(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5"))(length(breaksdiff)-1)
      }
    }
  }

#output = list(zlimdiff,breaksdiff,colorbardiff)
  output = list(datarange,breaksdiff,colorbardiff)
output

}



####
diff_colorramp = function(diffs,colorchoice,Blimit){
  
  #diffs=diffs[[3]];colorchoice=colorchoicediff;Blimit=BINLIMIT
  #colorchoice = "browntogreen"
  #Blimit = 20
  
  diffrange = range(diffs,na.rm=TRUE)
  diffrange[1]=floor(diffrange[1])
  if(diffrange[1] %% 2 != 0) diffrange[1]=diffrange[1]-1
  diffrange[2]=ceiling(diffrange[2])
  if(diffrange[2] %% 2 != 0) diffrange[2]=diffrange[2]+1
  
  if(diffrange[1]>=0){centerpoint = 0; startpoint=0; diffrange[1]=0}
  if(diffrange[1]<0){centerpoint = 0; startpoint=diffrange[1]; if(diffrange[2]<0){diffrange[2]=0}}
  breakcheck = 1
  breaklist = c(0.01,0.02,0.025,0.05,0.1,0.2,0.25,0.3,0.5,1,2,3,4,5,10,20,25,30,50)
  
  actualbins = diff(diffrange)/breaklist
  actidx = which(actualbins<Blimit)
  diffact = actualbins[actidx]-floor(actualbins[actidx])
  
  if(any(diffact==0)==TRUE){
    diffidx = which(diffact==0)
    breakcheck=actidx[diffidx[1]]
  } else {
    diffrange[1] = floor(diffrange[1]/10)*10
    diffrange[2] = ceiling(diffrange[2]/10)*10
    actualbins = diff(diffrange)/breaklist
    actidx = which(actualbins<Blimit)
    diffact = actualbins[actidx]-floor(actualbins[actidx])
    diffidx = which(diffact==0)
    breakcheck=actidx[diffidx[1]]
  }
  
  if(diffrange[2]==0 & colorchoice=="redtoblue") colorchoice="redtowhite"
  if(diffrange[2]==0 & colorchoice=="orangetopurple") colorchoice="orangetowhite"
  if(startpoint==0 & colorchoice=="bluetored") colorchoice="whitetored"
  if(startpoint==0 & colorchoice=="browntogreen") colorchoice="whitetogreen"
  if(startpoint==0 & colorchoice=="purpletoorange") colorchoice="whitetoorange"
  if(startpoint==0 & colorchoice=="orangetopurple") colorchoice="whitetopurple"
  
  zlimdiff = diffrange
  breaksdiff = c(seq(diffrange[1],diffrange[2],by=breaklist[breakcheck]))
  
  if(any(breaksdiff==0)==FALSE & zlimdiff[1]<0 & zlimdiff[2]>0){
    idx = which(abs(breaksdiff)==min(abs(breaksdiff)))
    if(length(idx)==1) breaksdiff[idx]=0
    if(length(idx)>1){
      breaksdiff = c(breaksdiff[1:(idx[1]-1)],0,breaksdiff[(idx[2]+1):length(breaksdiff)])
    }
  } 
  
  if(startpoint==centerpoint){
    if(colorchoice == "whitetored") colorbardiff = colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-1)
    if(colorchoice == "whitetogreen") colorbardiff = colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-1)
    if(colorchoice == "whitetoorange") colorbardiff = colorRampPalette(c("#f5f5f5","#fee0b6","#fdb863","#e08214","#b35806","#7f3b08"))(length(breaksdiff)-1)
    if(colorchoice == "whitetopurple") colorbardiff = colorRampPalette(c("#f5f5f5","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(breaksdiff)-1)
  } else {
    if(diffrange[2]>0){
      zeroidx = which(breaksdiff==0)
      if(colorchoice == "bluetored"){
        colorbardiff = c(colorRampPalette(c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-zeroidx))
        if((zeroidx-1)==1)colorbardiff = c("#f5f5f5",colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-zeroidx))
        } 
      if(colorchoice == "redtoblue"){
        colorbardiff = c(colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "browntogreen"){
        #colorbardiff = c(colorRampPalette(c("chocolate4","grey96"))(zeroidx-1),colorRampPalette(c("grey96","darkgreen"))(length(breaksdiff)-zeroidx))
        colorbardiff = c(colorRampPalette(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "greentobrown"){
        #colorbardiff = c(colorRampPalette(c("darkgreen","grey96"))(zeroidx-1),colorRampPalette(c("grey96","chocolate4"))(length(breaksdiff)-zeroidx))
        colorbardiff = c(colorRampPalette(c("#003c30","#01665e","#35978f","#80cdc1","#c7eae5","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#f6e8c3","#dfc27d","#bf812d","#8c510a","#543005"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "purpletoorange"){
        colorbardiff = c(colorRampPalette(c("#2d004b","#542788","#8073ac","#b2abd2","#d8daeb","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#fee0b6","#fdb863","#e08214","#b35806","#7f3b08"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "orangetopurple"){
        colorbardiff = c(colorRampPalette(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(breaksdiff)-zeroidx))
      } 
      
    } else {
      if(colorchoice == "redtowhite"){
        colorbardiff = colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(length(breaksdiff)-1)
      } 
      if(colorchoice == "orangetowhite"){
        colorbardiff = colorRampPalette(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f5f5f5"))(length(breaksdiff)-1)
      } 
      
      
    }
  }
  
  output = list(zlimdiff,breaksdiff,colorbardiff)
  output
  
}
