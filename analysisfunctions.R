#########################
#
# Analysis functions

source("/home/woot0002/scripts/climdex.R")

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
                            
  file1 = filename
  #message(file1)
  #message(varname)
                            
  if(varname=="tasmax"){
    message("tasmax provided, figuring out correct tasmin file")
    filesplit = do.call("c",strsplit(file1,"/",fixed=TRUE))
    nameend = substr(filesplit[length(filesplit)],nchar(varname)+1,nchar(filesplit[length(filesplit)]))
    #print(nameend)
    file2 = paste("/",filesplit[2],"/",filesplit[3],"/",filesplit[4],"/tasmin/",filesplit[6],"/tasmin",nameend,sep="")
    file2split = do.call("c",strsplit(file2,"txp",fixed=TRUE))
    file2 = paste(file2split[1],"tnp",file2split[2],sep="")
    message("file 1 is tasmax, file 2 is tasmin")
  }
                            
  if(varname=="tasmin"){
    filesplit = do.call("c",strsplit(file1,"/",fixed=TRUE))
    nameend = substr(filesplit[length(filesplit)],nchar(varname)+1,nchar(filesplit[length(filesplit)]))
    file2 = paste("/",filesplit[2],"/",filesplit[3],"/",filesplit[4],"/tasmax/",filesplit[6],"/tasmax",nameend,sep="")
    file2split = do.call("c",strsplit(file2,"tnp",fixed=TRUE))
    file2 = paste(file2split[1],"txp",file2split[2],sep="")
    message("file 1 is tasmin, file 2 is tasmax, switching the order")
    tmpfile1 = file2
    tmpfile2 = file1
    file1 = tmpfile1
    file2 = tmpfile2
    message("Switch complete")
  }
                            
  if(combofunction=="heatwaves"){
    splitname = do.call("c",strsplit(nameend,"_",fixed=TRUE))
    nameuse = paste(substr(splitname[3],9,15),"0",substr(splitname[3],17,nchar(splitname[3])),"_historical_",splitname[5],"_",splitname[6],"_19810101-20051231",sep="")
    q95list = system("ls /data2/3to5/I35/q95/*",intern=TRUE)
    idx = grep(nameuse,q95list)
    file3 = q95list[idx[1]] # tasmax historical q95
    file4 = q95list[idx[2]] # tasmin historical q95
  }
                            
  years = yearlydataperiod[1]:yearlydataperiod[2]
                            
  dates = datesdataperiod
                            
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(dates,1,4))==years[y])
                            
    test = nc_open(file1) # should always be tasmax file
    tempdata1 = ncvar_get(test,"tasmax",start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
                            
    if(y==1){
      yearlydat = array(NA,dim=c(dim(tempdata1)[1],dim(tempdata1)[2],length(years)))
      lat = ncvar_get(test,dimnames[2])
      lon = ncvar_get(test,dimnames[1])
      times = ncvar_get(test,dimnames[3])
      domainmask = ifelse(is.na(tempdata1[,,1])==FALSE,1,0)
      startdate = substr(test$dim[[4]]$units,12,21)
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
  message("Finished Calcs for year ",years[y])
  }
list(lon,lat,yearlydat)
}
                            

####

netcdfdroughtcalcs = function(filenames,varnames=c(),dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=seq(as.Date("1981-01-15"),as.Date("2005-12-15"),by="month"),appfunc="spi",scale=1){
  #filenames=filenames;varnames=varnamesin;dimnames=c("lon","lat","time");yearlydataperiod=c(1981,2005);datesdataperiod=datesin;appfunc=appfunc;scale=1
  
  files = unlist(strsplit(filenames,",",fixed=TRUE))
  
  test = nc_open(files[1])
  lat = ncvar_get(test,dimnames[2])
  lon = ncvar_get(test,dimnames[1])
  times = ncvar_get(test,dimnames[3])
  nc_close(test)
  
  years = yearlydataperiod[1]:yearlydataperiod[2]
  dates = datesdataperiod
  
  dsetout = array(NA,dim=c(length(lon),length(lat),length(dates)))
  
  for(r in 1:length(lon)){
    for(c in 1:length(lat)){
      
      if(appfunc=="spi"){
        test = nc_open(files[1])
        tempdata1 = ncvar_get(test,varnames[1],start=c(r,c,1),count=c(1,1,-1))
        tempdata1=tempdata1*86400
        nc_close(test)
        if(all(is.na(tempdata1)==TRUE)==FALSE){
          #tmp = spi(tempdata1,scale=scale,na.rm=TRUE)
          dsetout[r,c,]=spi(tempdata1,scale=scale,na.rm=TRUE)[[2]]
        }
      }
      
      if(appfunc=="spei"){
        test = nc_open(files[1]) #precip
        tempdata1 = ncvar_get(test,varnames[1],start=c(r,c,1),count=c(1,1,-1))
        if(varnames[1]=="pr") tempdata1=tempdata1*86400
        if(varnames[1]=="tasmax" | varnames[1]=="tasmin") tempdata1=tempdata1-273.15
        nc_close(test)
        
        test = nc_open(files[2]) #tasmax
        tempdata2 = ncvar_get(test,varnames[2],start=c(r,c,1),count=c(1,1,-1))
        if(varnames[2]=="pr") tempdata2=tempdata2*86400
        if(varnames[2]=="tasmax" | varnames[2]=="tasmin") tempdata1=tempdata1-273.15
        nc_close(test)
        
        test = nc_open(files[3]) #tasmin
        tempdata3 = ncvar_get(test,varnames[3],start=c(r,c,1),count=c(1,1,-1))
        if(varnames[3]=="pr") tempdata3=tempdata3*86400
        if(varnames[3]=="tasmax" | varnames[3]=="tasmin") tempdata1=tempdata1-273.15
        nc_close(test)
        
        if(all(is.na(tempdata1)==TRUE)==FALSE & all(is.na(tempdata2)==TRUE)==FALSE & all(is.na(tempdata3)==TRUE)==FALSE){
          ET = hargreaves(Tmin=matrix(tempdata3),Tmax=matrix(tempdata2),Pre=matrix(tempdata1),lat=lat[c],na.rm=TRUE)
          ########### Error in the above line
          
          DIFF = matrix(tempdata1)-matrix(ET)
          dsetout[r,c,]=spei(DIFF,scale=scale,na.rm=TRUE)[[2]]
        }
      }
     message("Finished calcs for row ",r," and column ",c) 
    }
  }
  list(lon,lat,times,dsetout)
}

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
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=0.254,outtype="max")
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="maxwetspell"){
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=0.254,outtype="max")
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
    if(season!="DJF") {
      if(season=="JJA") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=6 & as.numeric(substr(dates,6,7))<=8)
      if(season=="SON") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=9 & as.numeric(substr(dates,6,7))<=11)
      if(season=="MAM") yearidx = which(as.numeric(substr(dates,1,4))==years[y] & as.numeric(substr(dates,6,7))>=3 & as.numeric(substr(dates,6,7))<=5)
    } else {
      yearidx = which(substr(dates,1,4)==years[y] & (as.numeric(substr(dates,6,7))==12 | as.numeric(substr(dates,6,7))<=2))
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
    if(appliedfunction=="sum") {
      yearlydat[,,y] = apply(tempdata,c(1,2),sum,na.rm=TRUE)
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
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=0.254,outtype="max")
      yearlydat[,,y] = ifelse(domainmask==1,yearlydat[,,y],NA)
    }
    if(appliedfunction=="maxwetspell"){
      yearlydat[,,y] = apply(tempdata,c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=0.254,outtype="max")
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
  
  #inputdata=diffs[[3]]
  #colorchoice="bluetored"
  #Blimit = 30
  #type="difference"
  #use_fixed_scale = FALSE
  
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
  breaklist = c(0.0001,0.0002,0.00025,0.0005,0.001,0.002,0.0025,0.005,0.01,0.02,0.025,0.05,0.1,0.2,0.25,0.3,0.5,1,2,3,4,5,10,20,25,30,50,100,200,250,300,500,1000)
  
  actualbins = diff(datarange)/breaklist
  actidx = which(actualbins<Blimit)
  dataact = actualbins[actidx]-floor(actualbins[actidx])
  
  if(any(dataact<=1E-14)==TRUE){
    message("exact match for bins")
    dataidx = which(dataact<=1E-14)
    breakcheck=actidx[dataidx[1]]
  } else {
    message("no exact match going through while loop")
    checkpoint = any(dataact<=1E-14)
    counter=1
    while(checkpoint==FALSE){
      datarange[1] = floor(datarange[1]/(10^counter))*10^counter
      datarange[2] = ceiling(datarange[2]/(10^counter))*10^counter
      actualbins = diff(datarange)/breaklist
      actidx = which(actualbins<Blimit)
      dataact = actualbins[actidx]-floor(actualbins[actidx])
      dataidx = which(dataact==0)
      
      if(length(dataidx)>=1){
        breakcheck=actidx[dataidx[1]]
        checkpoint = any(dataact==0)
        break
      } else {
        counter=counter+1
        checkpoint = any(dataact==0)
      }
      
    }
  }
  
  if(datarange[1]==0 & colorchoice=="redtoblue") colorchoice="whitetoblue" # check this line
  if(datarange[2]==0 & colorchoice=="redtoblue") colorchoice="redtowhite" # check this line
  if(datarange[2]==0 & colorchoice=="orangetopurple") colorchoice="orangetowhite"
  if(startpoint==0 & centerpoint==0 & colorchoice=="bluetored") colorchoice="whitetored"
  if(startpoint==0 & centerpoint==0 & colorchoice=="browntogreen") colorchoice="whitetogreen"
  if(startpoint==0 & colorchoice=="purpletoorange") colorchoice="whitetoorange"
  if(startpoint==0 & colorchoice=="orangetopurple") colorchoice="whitetopurple"
  
  zlimdiff = datarange
  breaksdiff = c(seq(datarange[1],datarange[2],by=breaklist[breakcheck]))
  
  if(any(breaksdiff==centerpoint)==FALSE & zlimdiff[1]<centerpoint & zlimdiff[2]>centerpoint){
  idx = which(abs(breaksdiff)==min(abs(breaksdiff)))
  if(length(idx)==1) breaksdiff[idx]=centerpoint
  if(length(idx)>1){
  breaksdiff = c(breaksdiff[1:(idx[1]-1)],centerpoint,breaksdiff[(idx[2]+1):length(breaksdiff)])
  }
  } 
  
  message("zlimdiff: ",zlimdiff)
  message("centerpoint: ",centerpoint)
  message("startpoint: ",startpoint)
  message("colorchoice: ",colorchoice)
  
  if(startpoint==centerpoint){
  message("startpoint matches centerpoint")
  if(colorchoice == "whitetored") colorbardiff = colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-1)
  if(colorchoice == "yellowtored") colorbardiff = colorRampPalette(c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"))(length(breaksdiff)-1)
  
  if(colorchoice == "whitetogreen") colorbardiff = colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-1)
  if(colorchoice == "whitetoorange") colorbardiff = colorRampPalette(c("#f5f5f5","#fee0b6","#fdb863","#e08214","#b35806","#7f3b08"))(length(breaksdiff)-1)
  if(colorchoice == "whitetopurple") colorbardiff = colorRampPalette(c("#f5f5f5","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(breaksdiff)-1)  
  
  if(colorchoice == "whitetoblue") colorbardiff = colorRampPalette(c("#f5f5f5","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(length(breaksdiff)-1)
  
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
        colorbardiff = c(colorRampPalette(c("#2d004b","#542788","#8073ac","#b2abd2","#d8daeb","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#fee0b6","#fdb863","#e08214","#b35806","#7f3b08"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "orangetopurple"){
        colorbardiff = c(colorRampPalette(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(breaksdiff)-zeroidx))
      } 
      
    } else {
      message("odd if")
      if(colorchoice == "redtowhite"){
        colorbardiff = colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(length(breaksdiff)-1)
      } 
      
      if(colorchoice=="orangetowhite"){
        colorbardiff = colorRampPalette(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f5f5f5"))(length(breaksdiff)-1)
      }
      
    }
  }

output = list(zlimdiff,breaksdiff,colorbardiff)
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