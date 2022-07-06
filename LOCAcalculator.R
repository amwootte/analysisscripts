library(iterators)
#library(parallel)
#library(foreach)
#library(doParallel)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)
library(SPEI)
source("/data2/3to5/I35/scripts/analysisfunctions.R")

varname = "tmax95"
filesin = system(paste("ls /data4/data/DS_proj/LOCA/",varname,"/historical/*.nc",sep=""),intern=TRUE)
futureperiod = c(2036,2065)
histperiod = c(1950,2005)
GCM = exp = varlist = period = c()

colorchoiceuse = "bluetored"
Blimituse = 20

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesin[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  varin = paste(filesplit2[1],filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
  varlist[i]=filesplit2[1]
  GCM[i] = filesplit2[2]
  exp[i] = filesplit2[3]
  period[i] = substr(filesplit2[4],1,nchar(filesplit2[4])-3)
  
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
    timeoutmon= time
    
    yearmon = unique(substr(dates,1,7))
    
    if(i==1){
      DMdata = array(NA,dim=c(length(lon),length(lat),length(filesin)))
      yearidx = which(as.numeric(substr(dates,1,4))>=histperiod[1],as.numeric(substr(dates,1,4))<=histperiod[2])
    }
  }
  
  testdat = ncvar_get(nctest,varin,start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
  
  start_time <- Sys.time()
  DMout = apply(testdat,c(1,2),mean,na.rm=TRUE)
  DMdata[,,i] = ifelse(is.na(testdat[,,1])==FALSE,DMout,NA)
  end_time <- Sys.time()
  end_time - start_time
  
  nc_close(nctest)
  
  
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(filesin)," Time:",ptmend[3]," secs")
  
}

histtable=data.frame(varlist,GCM,exp,period)
DMhist = DMdata

####

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",varname,"/future/*.nc",sep=""),intern=TRUE)

GCM = exp = varlist = period = c()

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesin[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  varin = paste(filesplit2[1],filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
  varlist[i]=filesplit2[1]
  GCM[i] = filesplit2[2]
  exp[i] = filesplit2[3]
  period[i] = substr(filesplit2[4],1,nchar(filesplit2[4])-3)
  
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
    timeoutmon= time
    
    yearmon = unique(substr(dates,1,7))
    years = as.numeric(substr(dates,1,4))
    
    fileidx = which(years>=futureperiod[1] & years<=futureperiod[2])
    
    DMdata = array(NA,dim=c(length(lon),length(lat),length(filesin)))
    
  }
  
  testdat = ncvar_get(nctest,varin)
  
  start_time <- Sys.time()
  DMout = apply(testdat[,,fileidx],c(1,2),mean,na.rm=TRUE)
  DMdata[,,i] = ifelse(is.na(testdat[,,1])==FALSE,DMout,NA)
  end_time <- Sys.time()
  end_time - start_time
  
  nc_close(nctest)
  
  ptmend = proc.time()-ptm
  message("Finished calcs for file ",i," / ",length(filesin)," Time:",ptmend[3]," secs")
  
}

projtable=data.frame(varlist,GCM,exp,period)
DMproj = DMdata

####
# Calculate projected changes

rcp45idx = which(projtable$period=="rcp45")
rcp85idx = which(projtable$period=="rcp85")

rcp45change = DMproj[,,rcp45idx]-DMhist
rcp85change = DMproj[,,rcp85idx]-DMhist

rcp45meanchange = apply(rcp45change,c(1,2),mean,na.rm=TRUE)
rcp85meanchange = apply(rcp85change,c(1,2),mean,na.rm=TRUE)

####

mod.subV45 = mod.subV85 = list()
regions = list()
for(j in 1:6){
 
    shapefile=paste("/home/woot0002/shapefiles/PSS analysis unit ",j,sep="")

  test = readShapePoly(shapefile) # appropriate filename format
  projection(test) <- CRS("+init=epsg:5070")
  test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
  #projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
  
  regions[[j]]=test
  
  test.sub <- test # in here are two important details. 1) column name in the data array, 2) item in the data array to subset by
  tmpV = rcp45meanchange # need the 
  modrasV = raster(t(tmpV)[length(lat):1,])
  if(all(lon>0)){
    extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
  } else {
    extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
  }
  mod.subV45[[j]] <- crop(modrasV, extent(test.sub))
  mod.subV45[[j]] <- mask(modrasV, test.sub)
  
  tmpV = rcp85meanchange # need the 
  modrasV = raster(t(tmpV)[length(lat):1,])
  if(all(lon>0)){
    extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
  } else {
    extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
  }
  mod.subV85[[j]] <- crop(modrasV, extent(test.sub))
  mod.subV85[[j]] <- mask(modrasV, test.sub)
}

idxlon = which(lon-360>= -113 & lon-360<= -88)

rangeout = range(c(rcp45meanchange[idxlon,],rcp85meanchange[idxlon,]),na.rm=TRUE)
absoutmax = max(ceiling(abs(rangeout)))
rangeuse = c(-absoutmax,absoutmax)

diffcolorbar = colorramp(c(rcp45meanchange,rcp85meanchange),colorchoice=colorchoiceuse,Blimit=Blimituse,use_fixed_scale = TRUE,fixed_scale=c(0,130),type="difference") #,fixed_scale=fixed_scale

pdf(paste("/home/woot0002/PSS_",varname,"_",futureperiod[1],"-",futureperiod[2],".pdf",sep=""),width=10,height=7)

par(mfrow=c(1,2))

testsfc = list(x=lon-360,y=lat,z=rcp45meanchange)
surface(testsfc,type="I",xlim=c(-113,-88),main=paste(futureperiod[1],"-",futureperiod[2]," RCP 4.5\n Mean Projected Change ",varname,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("world",add=TRUE) # 
map("state",add=TRUE,col="darkgray")
for(j in 1:6){
  plot(regions[[j]],add=TRUE)
}

testsfc = list(x=lon-360,y=lat,z=rcp85meanchange)
surface(testsfc,type="I",xlim=c(-113,-88),main=paste(futureperiod[1],"-",futureperiod[2]," RCP 8.5\n Mean Projected Change ",varname,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("world",add=TRUE) #
map("state",add=TRUE,col="darkgray")
for(j in 1:6){
  plot(regions[[j]],add=TRUE)
}

dev.off()

####
# Analysis unit grabs for barplots, grab mean first

meanval45 = meanval85 = c()

for(j in 1:6){
testin = matrix(getValues(mod.subV45[[j]]),nrow=length(lon),ncol=length(lat))
testin = testin[,length(lat):1]
meanval45[j]=mean(testin,na.rm=TRUE)

testin = matrix(getValues(mod.subV85[[j]]),nrow=length(lon),ncol=length(lat))
testin = testin[,length(lat):1]
meanval85[j]=mean(testin,na.rm=TRUE)
}

####
# grab analysis unit values for each model and calculate standard deviation

sdval45 = sdval85 = c()

for(j in 1:6){
  
  modvals45=modvals85=c()
  
  for(m in 1:30){
    
    shapefile=paste("/home/woot0002/shapefiles/PSS analysis unit ",j,sep="")
    
    test = readShapePoly(shapefile) # appropriate filename format
    projection(test) <- CRS("+init=epsg:5070")
    test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
    #projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
    
    regions[[j]]=test
    
    test.sub <- test # in here are two important details. 1) column name in the data array, 2) item in the data array to subset by
    tmpV = rcp45change[,,m] # need the 
    modrasV = raster(t(tmpV)[length(lat):1,])
    if(all(lon>0)){
      extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
    } else {
      extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
    }
    mod.subV45 <- crop(modrasV, extent(test.sub))
    mod.subV45 <- mask(modrasV, test.sub)
    
    tmpV = rcp85change[,,m] # need the 
    modrasV = raster(t(tmpV)[length(lat):1,])
    if(all(lon>0)){
      extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
    } else {
      extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
    }
    mod.subV85 <- crop(modrasV, extent(test.sub))
    mod.subV85 <- mask(modrasV, test.sub)
    
    testin = matrix(getValues(mod.subV45),nrow=length(lon),ncol=length(lat))
    modvals45[m] = mean(testin,na.rm=TRUE)
    
    testin = matrix(getValues(mod.subV85),nrow=length(lon),ncol=length(lat))
    modvals85[m] = mean(testin,na.rm=TRUE)
  }
  sdval45[j] = sd(modvals45,na.rm=TRUE)
  sdval85[j] = sd(modvals85,na.rm=TRUE)
}

####

meanval = c(meanval45,meanval85)
sdval = c(sdval45,sdval85)
scen = rep(c("rcp45","rcp85"),each=6)
unit = rep(1:6,2)

dataindat = data.frame(unit,scen,meanval,sdval)
dataindat$ci = (dataindat$sd/sqrt(30))*1.96

ggplot(dataindat, aes(x=scen, y=meanval, fill=scen)) + 
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=meanval-ci, ymax=meanval+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  scale_fill_manual(values=c("#92c5de", "#b2182b"))+facet_wrap(unit,nrow=2,ncol=3)+
  ggtitle(paste(futureperiod[1],"-",futureperiod[2]," Projected Change in ",varname,sep=""))+xlab("scenario")+ylab("Projected Change")
ggsave(paste("/home/woot0002/PSS_barplot_",varname,"_",futureperiod[1],"-",futureperiod[2],".pdf",sep=""),width=10,height=5)


