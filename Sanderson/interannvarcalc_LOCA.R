#######################
#
# Monthly climo calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = "pr"

filelist = read.table(paste("/home/woot0002/RCMES/GCMfilelist_",var,"_v2.csv",sep=""),sep=",",header=TRUE,colClasses = "character")
load("/home/woot0002/DS_ind/manuscript1/GCMlist.Rdata")

filelistin = filelist[which(filelist$model %in% GCMlist),]

GCMtable = filelistin
GCMtable$GCMs = paste(GCMtable$model,GCMtable$ensemble,sep="_")

######
# set up for regrid

filein = "/home/woot0002/monthlyclimo/tasmax_day_I35txdetrp1-DeltaSD-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231_mon.nc"
nctest = nc_open(filein)
latobs = ncvar_get(nctest,"lat")
lonobs = ncvar_get(nctest,"lon")
varobs = ncvar_get(nctest,nctest$var[[1]]$name,start=c(1,1,1),count=c(-1,-1,1)) # serves to help with a mask
nc_close(nctest)

histdatrg = array(NA,dim=c(length(lonobs),length(latobs),length(1981:2005),26))
projdatrg = array(NA,dim=c(length(lonobs),length(latobs),length(2070:2099),26))

######
# Historical

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",var,"/historical/*.nc",sep=""),intern=TRUE)

LOCAtable = do.call("rbind",strsplit(filesin,split="_",fixed=TRUE))

GCMtablehist = GCMtable[which(GCMtable$experiment=="historical"),]
GCMs = unique(paste(GCMtablehist$model,GCMtablehist$ensemble,sep="_"))

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
monyear = dates[which(substr(dates,9,10)=="01")]
years = 1981:2005
months = 1:12

for(g in 1:length(GCMs)){
  
  idx=grep(GCMs[g],filesin)
  
  test = nc_open(filesin[idx])
  lon = ncvar_get(test,"lon")
  lat = ncvar_get(test,"lat")
  
  startdate = test$dim[[3]]$units
  startdate = as.Date(substr(startdate,12,21))
  time = ncvar_get(test,"time")
  datesin = startdate+time
  
  for(y in 1:length(years)){
    
    yearidx = which(as.numeric(substr(datesin,1,4))==years[y])
    
    vardata = ncvar_get(test,test$var[[1]]$name,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
    
    if(var=="pr"){
      tmpagg=apply(vardata,c(1,2),sum,na.rm=TRUE)
    }
    if(var=="tasmax"){
      tmpagg=apply(vardata,c(1,2),mean,na.rm=TRUE)
    }
    
    int1 = interp.surface.grid(list(x=lon,y=lat,z=tmpagg),grid.list=list(x=lonobs,y=latobs))
    newmatrix = matrix(int1[[3]],nrow=length(lonobs),ncol=length(latobs))
    if(var=="pr"){
      newmatrix = ifelse(newmatrix<0,0,newmatrix)
    }
    histdatrg[,,y,g] = ifelse(is.na(varobs)==FALSE,newmatrix,NA)
    
    message("Finished calcs for year ",y," / ",length(years))
  }
  nc_close(test)
  message("Finished calcs for GCM ",g," / ",length(GCMs))
}

######
# Future

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",var,"/future/*.nc",sep=""),intern=TRUE)

LOCAtable = do.call("rbind",strsplit(filesin,split="_",fixed=TRUE))

GCMtableproj = GCMtable[which(GCMtable$experiment=="rcp85"),]
GCMs = unique(paste(GCMtableproj$model,GCMtableproj$ensemble,sep="_"))

dates = seq(as.Date("2070-01-01"),as.Date("2099-12-31"),by="day")
monyear = dates[which(substr(dates,9,10)=="01")]
years = 2070:2099
months = 1:12

for(g in 1:length(GCMs)){
  
  idx=grep(GCMs[g],filesin)
  filesuse = filesin[idx]
  idxscen = grep("rcp85",filesuse)
  filesuse2 = filesuse[idxscen]
  
  for(y in 1:length(years)){
    
    idx2 = grep(years[y],filesuse2)
    
    test = nc_open(filesuse2[idx2])
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
    
    startdate = test$dim[[3]]$units
    startdate = as.Date(substr(startdate,12,21))
    time = ncvar_get(test,"time")
    datesin = startdate+time
    
    yearidx = which(as.numeric(substr(datesin,1,4))==years[y])
    
    vardata = ncvar_get(test,test$var[[1]]$name,start=c(1,1,1),count=c(-1,-1,-1))
    if(var=="pr"){
      vardata = vardata*86400
    }
   
    if(var=="pr"){
      tmpagg=apply(vardata,c(1,2),sum,na.rm=TRUE)
    }
    if(var=="tasmax"){
      tmpagg=apply(vardata,c(1,2),mean,na.rm=TRUE)
    }
    
    int1 = interp.surface.grid(list(x=lon,y=lat,z=tmpagg),grid.list=list(x=lonobs,y=latobs))
    newmatrix = matrix(int1[[3]],nrow=length(lonobs),ncol=length(latobs))
    if(var=="pr"){
      newmatrix = ifelse(newmatrix<0,0,newmatrix)
    }
    projdatrg[,,y,g] = ifelse(is.na(varobs)==FALSE,newmatrix,NA)
    message("Finished calcs for year ",y," / ",length(years))
    nc_close(test)
  }
  
  message("Finished calcs for GCM ",g," / ",length(GCMs))
}


#######
# yearly domain averages by model

#full domain yearly average by model

fullavghist = apply(histdatrg,c(3,4),mean,na.rm=TRUE)
# rows are years, columns are models
fullavgproj = apply(projdatrg,c(3,4),mean,na.rm=TRUE)

# masks
test = nc_open("/home/woot0002/DS_ind/louisiana_mask.nc")
LAregionmask = ncvar_get(test,"mask")
nc_close(test)

test = nc_open("/home/woot0002/DS_ind/new mexico_mask.nc")
NMregionmask = ncvar_get(test,"mask")
nc_close(test)

# mask data

histdatrgNM = histdatrgLA = histdatrg
projdatrgNM = projdatrgLA = projdatrg

for(g in 1:26){
  for(i in 1:25){
    histdatrgNM[,,i,g] = ifelse(NMregionmask==1,histdatrg[,,i,g],NA)
    histdatrgLA[,,i,g] = ifelse(LAregionmask==1,histdatrg[,,i,g],NA)
  }
}

for(g in 1:26){
  for(i in 1:30){
    projdatrgNM[,,i,g] = ifelse(NMregionmask==1,projdatrg[,,i,g],NA)
    projdatrgLA[,,i,g] = ifelse(LAregionmask==1,projdatrg[,,i,g],NA)
  }
}

fullavghistNM = apply(histdatrgNM,c(3,4),mean,na.rm=TRUE)
# rows are years, columns are models
fullavgprojNM = apply(projdatrgNM,c(3,4),mean,na.rm=TRUE)

fullavghistLA = apply(histdatrgLA,c(3,4),mean,na.rm=TRUE)
# rows are years, columns are models
fullavgprojLA = apply(projdatrgLA,c(3,4),mean,na.rm=TRUE)

###
# interannual sd calcs

GCM = rep(unique(GCMtableproj$model),3)
region = rep(c("full","NM","LA"),each=26)
variable = rep(var,26*3)
histintann = c(apply(fullavghist,2,sd),apply(fullavghistNM,2,sd),apply(fullavghistLA,2,sd))
projintann = c(apply(fullavgproj,2,sd),apply(fullavgprojNM,2,sd),apply(fullavgprojLA,2,sd))

intanntab = data.frame(GCM,region,variable,histintann,projintann)

write.table(intanntab,paste("/home/woot0002/DS_ind/intanntab_LOCA_",var,".csv",sep=""),sep=",",row.names=FALSE)

