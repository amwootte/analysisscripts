library(ncdf4)
source("/data2/3to5/I35/scripts/analysisfunctions.R")

varname = "pr"
#loc_name = "Deming"
#loc_lon = 360-107.7586
#loc_lat = 32.2687

#loc_name = "SanAntonio"
#loc_lon = 360-106.8659
#loc_lat = 33.9178

#loc_name = "Albuquerque"
#loc_lon = 360-106.6504
#loc_lat = 35.0844

blocksize = 5

load("/data2/3to5/I35/scripts/prhistfiles.Rdata")
load("/data2/3to5/I35/scripts/prprojfiles.Rdata")
varin = varname

histfilelist = paste("/data4/data/DS_proj/LOCA/",varin,"/historical/",c(
  "pr_CCSM4_r6i1p1_historical.nc","pr_MIROC5_r1i1p1_historical.nc","pr_MPI-ESM-LR_r1i1p1_historical.nc"),sep="")
  
projfilelist = paste("/data4/data/DS_proj/LOCA/",varin,"/future/",c(
  "pr_CCSM4_r6i1p1_rcp45-2071-2080.nc", "pr_CCSM4_r6i1p1_rcp45-2081-2090.nc", "pr_CCSM4_r6i1p1_rcp45-2091-2100.nc", 
  "pr_CCSM4_r6i1p1_rcp85-2071-2080.nc", "pr_CCSM4_r6i1p1_rcp85-2081-2090.nc", "pr_CCSM4_r6i1p1_rcp85-2091-2100.nc",
  "pr_MIROC5_r1i1p1_rcp45-2071-2080.nc", "pr_MIROC5_r1i1p1_rcp45-2081-2090.nc", "pr_MIROC5_r1i1p1_rcp45-2091-2100.nc",
  "pr_MIROC5_r1i1p1_rcp85-2071-2080.nc", "pr_MIROC5_r1i1p1_rcp85-2081-2090.nc", "pr_MIROC5_r1i1p1_rcp85-2091-2100.nc",
  "pr_MPI-ESM-LR_r1i1p1_rcp45-2071-2080.nc", "pr_MPI-ESM-LR_r1i1p1_rcp45-2081-2090.nc", "pr_MPI-ESM-LR_r1i1p1_rcp45-2091-2100.nc",
  "pr_MPI-ESM-LR_r1i1p1_rcp85-2071-2080.nc", "pr_MPI-ESM-LR_r1i1p1_rcp85-2081-2090.nc", "pr_MPI-ESM-LR_r1i1p1_rcp85-2091-2100.nc"),sep="")

#############
# 1a- historical calcs

dates = seq(as.Date("1950-01-01"),as.Date("2005-12-31"),by="day")
years = 1981:2005

for(i in 1:length(histfilelist)){
  
  #i=1
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])

  test = nc_open(histfilelist[i])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  nc_close(test)
  
  if(i==1){
    histdates = dates
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
    tempdata = ncvar_get(test,test$var[[1]]$name,start=c(locstart[1],locstart[2],1),count=c(blocksize,blocksize,-1))
    nc_close(test)
    tmp = apply(tempdata,3,mean,na.rm=TRUE)
    
    #### dates to use
    if(length(dates)==length(tmp)){
      rmidx = which(substr(dates,6,10)=="02-29")
      tmp = tmp[-rmidx]
      datestmp = histdates[-rmidx]
      finaldates = which(as.numeric(substr(datestmp,1,4))>=1981)
      datestmp = datestmp[finaldates]
      tmp = tmp[finaldates]
      if(i>=2) tmp = tmp*86400
    }
    rm(tempdata)

    if(i==1) yearframe = array(NA,dim=c(365,ncol=length(years),length(histfilelist)))
    
    for(y in 1:length(years)){
      yearidx = which(as.numeric(substr(datestmp,1,4))==years[y])
      yearframe[1:365,y,i]= tmp[yearidx]
    }

  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

histyearframe = yearframe
histdates = datestmp

#############
# 1b- Future calcs

years = 2071:2099

GCMS = c("CCSM4","MIROC5","MPI-ESM-LR")

for(i in 1:length(GCMS)){
  
  filesin = grep(GCMS[i],projfilelist)
  
  ptm = proc.time()
  message("Starting work on file ",GCMS[i])
  
  test = nc_open(projfilelist[filesin[1]])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  times = ncvar_get(test,"time")
  nc_close(test)
  
  if(i==1){
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
  
  dates = seq(as.Date("2071-01-01"),as.Date("2100-12-31"),by="day")
  tmp = c()
  for(j in 1:length(filesin)){
    test = nc_open(projfilelist[filesin[j]])
    tempdata = ncvar_get(test,test$var[[1]]$name,start=c(locstart[1],locstart[2],1),count=c(blocksize,blocksize,-1))
    nc_close(test)
    tmp = c(tmp,apply(tempdata,3,mean,na.rm=TRUE))
    rm(tempdata)
  }
  
  tmp1 = tmp[1:length(dates)]
  tmp2 = tmp[(length(dates)+1):length(tmp)]
  
  if(length(dates)==length(tmp1)){
    rmidx = which(substr(dates,6,10)=="02-29")
    tmp1 = tmp1[-rmidx]
    tmp2 = tmp2[-rmidx]
    datestmp = dates[-rmidx]
    finaldates = which(as.numeric(substr(datestmp,1,4))>=2071 & as.numeric(substr(datestmp,1,4))<=2099)
    datestmp = datestmp[finaldates]
    tmp1 = tmp1[finaldates]
    tmp2 = tmp2[finaldates]
  }
  
  if(i==1) yearframe = array(NA,dim=c(365,ncol=length(years),length(GCMS)*2))
  
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(datestmp,1,4))==years[y])
    yearframe[1:365,y,i+(i-1)]= tmp1[yearidx]
    yearframe[1:365,y,i+(i)]= tmp2[yearidx]
  }
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(GCMS))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

projyearframe = yearframe
projdates = datestmp

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
plot(histdailycumanom[,f],type="l",main = paste(GCMS[f],"LOCA",sep=" "))
}

#####

projyearmean = apply(projyearframe,c(2,3),mean,na.rm=TRUE) # in this object, rows are each year, columns are each file
endyears = 2071:2099
projyears = 2071:2099

endidx = which(projyears %in% endyears)
projclimo_end = apply(projyearmean[endidx,],2,mean,na.rm=TRUE) # in this object, each value is the climo for each file
projdailymean_end = apply(projyearframe[,endidx,],c(1,3),mean,na.rm=TRUE) 

#histdailyanom = histdailymean-histclimo

projdailyanom_end = projdailymean_end

GCMscens = paste(rep(c("rcp45","rcp85"),each=3),rep(GCMS,2),sep="_")

for(i in 1:length(GCMscens)){
  projdailyanom_end[,i] = projdailymean_end[,i]-projclimo_end[i]
}

projdailycumanom_end = projdailymean_end
for(f in 1:length(GCMscens)){
  for(i in 1:365){
    if(i==1){
      projdailycumanom_end[1,f] = projdailyanom_end[1,f]
    } else {
      projdailycumanom_end[i,f] = projdailyanom_end[i,f]+projdailycumanom_end[(i-1),f]
    }
  }
}

par(mfrow=c(2,1))

for(f in 1:length(GCMscens)){
  plot(projdailycumanom_end[,f],type="l",main = paste("2071-2099",GCMscens[f],"LOCA",sep=" "))
}

save(list=c("GCMS","GCMscens","histdailycumanom","projdailycumanom_end","histyears","endyears","projyears"),file=paste("wetseasontesting_LOCA_",loc_name,".Rdata",sep=""))







