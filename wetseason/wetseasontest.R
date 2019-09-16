library(ncdf4)
source("/data2/3to5/I35/scripts/analysisfunctions.R")

varname = "pr"
#loc_name = "Deming"
#loc_lon = 360-107.7586
#loc_lat = 32.2687

#loc_name = "SanAntonito"
#loc_lon = 360-106.8659
#loc_lat = 33.9178

loc_name = "Albuquerque"
loc_lon = 360-106.6504
loc_lat = 35.0844

blocksize = 3

load("/data2/3to5/I35/scripts/prhistfiles.Rdata")
load("/data2/3to5/I35/scripts/prprojfiles.Rdata")

histfilelist = do.call("c",strsplit(histfilelist,",",fixed=TRUE))
projfilelist = do.call("c",strsplit(projfilelist,",",fixed=TRUE))

varin = varname

# Creating file breakdown tables
filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),9)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),3)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

#############
# 1a- historical calcs

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
years = 1981:2005

for(i in 1:length(histfilelist)){
  
  #i=1
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  #### dates to use
  if(noleap==FALSE){
    datesin = dates
    rmidx = which(substr(dates,6,10)=="02-29")
  }
  if(noleap==TRUE){
    rmidx = which(substr(dates,6,10)=="02-29")
    datesin = dates[-rmidx]
}
  
  test = nc_open(histfilelist[i])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  if(i==1){
    histdates = datesin
    ###
    # create model grid
    
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
    locstart = c(pointarea$R-floor(blocksize/2),pointarea$C-floor(blocksize/2))
    locend = c(pointarea$R+floor(blocksize/2),pointarea$C+floor(blocksize/2))
  }
  
  test = nc_open(histfilelist[i])
    tempdata = ncvar_get(test,varin,start=c(locstart[1],locstart[2],1),count=c(blocksize,blocksize,-1))
    nc_close(test)
    tmp = apply(tempdata,3,mean,na.rm=TRUE)*86400
    if(noleap==TRUE) tmp = tmp[-rmidx]
    if(noleap == FALSE){
      tmp = tmp[-rmidx]
      datesin = datesin[-rmidx]
    }
    rm(tempdata)

    if(i==1) yearframe = array(NA,dim=c(365,ncol=length(years),length(histfilelist)))
    
    for(y in 1:length(years)){
      yearidx = which(as.numeric(substr(datesin,1,4))==years[y])
      yearframe[1:365,y,i]= tmp[yearidx]
    }

  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

histyearframe = yearframe
histdates = datesin

#############
# 1b- Future calcs

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
years = 2006:2099

for(i in 1:length(projfilelist)){
  
  #i=1
  ptm = proc.time()
  message("Starting work on file ",projfilelist[i])
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  #### dates to use
  if(noleap==FALSE){
    datesin = dates
    rmidx = which(substr(dates,6,10)=="02-29")
  }
  if(noleap==TRUE){
    rmidx = which(substr(dates,6,10)=="02-29")
    datesin = dates[-rmidx]
  }
  
  test = nc_open(projfilelist[i])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  nc_close(test)
  if(i==1){
    histdates = datesin
    ###
    # create model grid
    
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
    locend = c(pointarea$R+1,pointarea$C+1)
  }
  
  test = nc_open(projfilelist[i])
  tempdata = ncvar_get(test,varin,start=c(locstart[1],locstart[2],1),count=c(3,3,-1))
  nc_close(test)
  tmp = apply(tempdata,3,mean,na.rm=TRUE)*86400
  if(noleap==TRUE) tmp = tmp[-rmidx]
  if(noleap == FALSE){
    tmp = tmp[-rmidx]
    datesin = datesin[-rmidx]
  }
  rm(tempdata)
  
  if(i==1) yearframe = array(NA,dim=c(365,ncol=length(years),length(projfilelist)))
  
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(datesin,1,4))==years[y])
    yearframe[1:365,y,i]= tmp[yearidx]
  }
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

projyearframe = yearframe
projdates = datesin

#####

histyearmean = apply(histyearframe,c(2,3),mean,na.rm=TRUE) # in this object, rows are each year, columns are each file
histclimo = apply(histyearmean,2,mean,na.rm=TRUE) # in this object, each value is the climo for each file
histyears = 1981:2005
histdailymean = apply(histyearframe,c(1,3),mean,na.rm=TRUE) 
#histdailyanom = histdailymean-histclimo

histdailyanom = histdailymean
for(i in 1:length(histfilelist)){
  histdailyanom[,i] = histdailymean[,i]-histclimo[i]
}

histdailycumanom = histdailymean
for(f in 1:length(histfilelist)){
    for(i in 1:365){
      if(i==1){
        histdailycumanom[1,f] = histdailyanom[1,f]
      } else {
        histdailycumanom[i,f] = histdailyanom[i,f]+histdailycumanom[(i-1),f]
      }
    }
  }

for(f in 1:length(histfilelist)){
plot(histdailycumanom[,f],type="l",main = paste(histfilebreakdown$GCM[f],histfilebreakdown$DS[f],histfilebreakdown$obs[f],sep=" "))
}

#####

projyearmean = apply(projyearframe,c(2,3),mean,na.rm=TRUE) # in this object, rows are each year, columns are each file
midyears = 2041:2070
endyears = 2071:2099
projyears = 2006:2099

mididx = which(projyears %in% midyears)
endidx = which(projyears %in% endyears)

projclimo_mid = apply(projyearmean[mididx,],2,mean,na.rm=TRUE) # in this object, each value is the climo for each file
projclimo_end = apply(projyearmean[endidx,],2,mean,na.rm=TRUE) # in this object, each value is the climo for each file

projdailymean_mid = apply(projyearframe[,mididx,],c(1,3),mean,na.rm=TRUE) 
projdailymean_end = apply(projyearframe[,endidx,],c(1,3),mean,na.rm=TRUE) 

#histdailyanom = histdailymean-histclimo

projdailyanom_mid = projdailymean_mid
projdailyanom_end = projdailymean_end

for(i in 1:length(projfilelist)){
  projdailyanom_mid[,i] = projdailymean_mid[,i]-projclimo_mid[i]
  projdailyanom_end[,i] = projdailymean_end[,i]-projclimo_end[i]
}

projdailycumanom_mid = projdailymean_mid
projdailycumanom_end = projdailymean_end
for(f in 1:length(projfilelist)){
  for(i in 1:365){
    if(i==1){
      projdailycumanom_mid[1,f] = projdailyanom_mid[1,f]
      projdailycumanom_end[1,f] = projdailyanom_end[1,f]
    } else {
      projdailycumanom_mid[i,f] = projdailyanom_mid[i,f]+projdailycumanom_mid[(i-1),f]
      projdailycumanom_end[i,f] = projdailyanom_end[i,f]+projdailycumanom_end[(i-1),f]
    }
  }
}

par(mfrow=c(2,1))

for(f in 1:length(projfilelist)){
  plot(projdailycumanom_mid[,f],type="l",main = paste("2041-2070",projfilebreakdown$GCM[f],projfilebreakdown$scen[f],projfilebreakdown$DS[f],projfilebreakdown$obs[f],sep=" "))
  plot(projdailycumanom_end[,f],type="l",main = paste("2071-2099",projfilebreakdown$GCM[f],projfilebreakdown$scen[f],projfilebreakdown$DS[f],projfilebreakdown$obs[f],sep=" "))
}

save(list=c("histfilebreakdown","projfilebreakdown","histyearframe","histdailycumanom","projyearframe","projdailycumanom_mid","projdailycumanom_end","histyears","midyears","endyears","projyears"),file=paste("wetseasontesting_",loc_name,".Rdata",sep=""))







