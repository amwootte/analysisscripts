
#ptm= Sys.time()
#.libPaths('/home/amwootte/R/x86_64-redhat-linux-gnu-library/3.3/')
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

source("/data2/3to5/I35/scripts/analysisfunctions.R")

#cl <- makePSOCKcluster(20,outfile="")#makePSOCKcluster(20)
#registerDoParallel(cl)

##########
# Set input information

message("Setting inputs")

shapefile = "/home/woot0002/shapefiles/GCWA_Recovery_Units_1992"

varname = "pr" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2036,2065)

colorchoicechange = "browntogreen"
BINSchange = 20
USEFS = FALSE
colorchoicehist = "whitetogreen"
BINShist = 20

lonbnds = c(-101.5,-96.5)
latbnds = c(29,33)

varin = varname
if(varname=="tmax85" | varname=="tmax90" | varname=="tmax95" | varname=="tmax100" | varname=="gsl") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="prcptot" | varname=="prmean" | varname=="prmed" | varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd" | varname=="R99" | varname=="R95" | varname=="R90"  | varname=="R99freq" | varname=="R95freq" | varname=="R90freq") varin="pr"

########
# Find all file names

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data4/data/DS_proj/MACA_METDATA/",varin,"/historical/*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data4/data/DS_proj/MACA_METDATA/",varin,"/future/*rcp85*.nc",sep=""),intern=T)

uhistfiles = histfilelist[grep("1981",histfilelist)]
filebreakdown = do.call(rbind,strsplit(uhistfiles,"_",fixed=TRUE))
GCMs = filebreakdown[,5]
rm(filebreakdown)

uhistfiles = histfilelist[grep(GCMs[1],histfilelist)]
filebreakdown = do.call(rbind,strsplit(uhistfiles,"_",fixed=TRUE))
histyears = unique(as.numeric(substr(filebreakdown[,8],1,4)))

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = filebreakdown[,c(5)]
filebreakdown2a = filebreakdown[,c(8,9)]
filebreakdown2a = data.frame(filebreakdown2a)
filebreakdown2a[,1] = as.numeric(as.character(filebreakdown2a[,1]))
filebreakdown2a[,2] = as.numeric(as.character(filebreakdown2a[,2]))
filebreakdown2 = data.frame(filebreakdown2)
filebreakdown3 = cbind(filebreakdown2,filebreakdown2a)
names(filebreakdown3) = c("GCM","startyear","endyear")
filebreakdown3[,2] = as.numeric(filebreakdown3[,2])
filebreakdown3[,3] = as.numeric(filebreakdown3[,3])
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filecheck = c()
for(i in 1:nrow(projfilebreakdown)){
  filecheck[i] = ifelse(any(as.numeric(projfilebreakdown[i,2]):as.numeric(projfilebreakdown[i,3]) %in% futureperiod[1]:futureperiod[2])==TRUE,1,0)
}

projfilebreakdown=projfilebreakdown[which(filecheck==1),]
projfilelist=projfilelist[which(filecheck==1)]

###########

test = nc_open(histfilelist[1])
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

lonin = lon[which(lon >= lonbnds[1] & lon<=lonbnds[2])]
lonidx = which(lon >= lonbnds[1] & lon<=lonbnds[2])
latin = lat[which(lat >= latbnds[1] & lat<=latbnds[2])]
latidx = which(lat >= latbnds[1] & lat<=latbnds[2])

###
# get shapefile
message("opening URG shapefile")
test = readShapePoly(shapefile)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

###
# historical climatology calcs with MACA

climoout = array(NA,dim=c(length(lonin),length(latin),length(GCMs)))

for(i in 1:length(GCMs)){
  
  ptm=proc.time()
histfilesuse = histfilelist[grep(paste(GCMs[i],"_",sep=""),histfilelist)]
filetmp = do.call(rbind,strsplit(histfilesuse,"_",fixed=TRUE))
tmpdat = array(NA,dim=c(length(lonin),length(latin),length(histfilesuse)))  

for(j in 1:length(histfilesuse)){
  #filevar = paste(filetmp[j,4],filetmp[j,5],filetmp[j,6],filetmp[j,7],sep="_")
  #if(i==8 & j==38){
  #  filevar="precipitation"
  #}
  nctest = nc_open(histfilesuse[j])
  tmpdata = ncvar_get(nctest,nctest$var[[1]]$name,start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))
  nc_close(nctest)
  if(varname=="pr"){
    tmpdat[,,j] = apply(tmpdata,c(1,2),sum,na.rm=TRUE)
  }
  message("Finished calcs for year ",j," / ",length(histfilesuse))
}
  
  climoout[,,i] = apply(tmpdat,c(1,2),mean,na.rm=TRUE)
  ptmend = proc.time()-ptm
  message("Finished historical calcs for GCM ",i," / ",length(GCMs))
  message("Time taken is ",ptmend[3]," seconds")
}

###
# future climatology calcs with MACA

climooutp = array(NA,dim=c(length(lonin),length(latin),length(GCMs)))

for(i in 1:length(GCMs)){
  
  projfilesuse = projfilelist[grep(paste(GCMs[i],"_",sep=""),projfilelist)]
  filetmp = do.call(rbind,strsplit(projfilesuse,"_",fixed=TRUE))
  tmpdat = array(NA,dim=c(length(lonin),length(latin),length(futureperiod[1]:futureperiod[2])))  
  
  for(j in 1:length(projfilesuse)){
    #filevar = paste(filetmp[j,4],filetmp[j,5],filetmp[j,6],filetmp[j,7],sep="_")
    nctest = nc_open(projfilesuse[j])
    tmpdata = ncvar_get(nctest,nctest$var[[1]]$name,start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1))
    nc_close(nctest)
    
    years = as.numeric(filetmp[j,8]):as.numeric(filetmp[j,9])
    dates = seq(as.Date(paste(years[1],"-01-01",sep="")),as.Date(paste(years[length(years)],"-12-31",sep="")),by="day")
    
    if(length(dates)>dim(tmpdata)[[3]]){
      dates = dates[-which(substr(dates,6,10)=="02-29")]
    }
    
    for(y in 1:length(years)){
      yearidx = which(as.numeric(substr(dates,1,4))==years[y])
      if(varname=="pr"){
        tmpdat[,,(years[y]-(futureperiod[1]-1))] = apply(tmpdata[,,yearidx],c(1,2),sum,na.rm=TRUE)
      }
    }
    message("Finished calcs for ",j," / ",length(projfilesuse))
  }
  
  climooutp[,,i] = apply(tmpdat,c(1,2),mean,na.rm=TRUE)
  message("Finished projected calcs for GCM ",i," / ",length(GCMs))
}

######
# Calc anomalies

anomout = climooutp-climoout

save(list=c("anomout","climooutp","climoout","GCMs","futureperiod","varname"),file=paste("/home/woot0002/GCWregion_rcp85_pr_",futureperiod[1],"-",futureperiod[2],".pdf",sep=""))

ensmean = apply(anomout,c(1,2),mean,na.rm=TRUE)
ensmin = apply(anomout,c(1,2),min,na.rm=TRUE)
ensmax = apply(anomout,c(1,2),max,na.rm=TRUE)


######
##########
# Plot maps

pdf(paste("/home/woot0002/GCWregion_rcp85_pr_changes.pdf",sep=""),width=10,height=5,onefile=TRUE)


diffcolorbar = colorramp(c(ensmin,ensmean,ensmax),colorchoice="browntogreen",Blimit=20,type="difference",use_fixed_scale = FALSE)

testsfc = list(x=lon[lonidx],y=lat[latidx],z=ensmin)
surface(testsfc,type="I",main="GCW units overlaid with RCP8.5 2036-2065 \n Ensemble Lowest Precipitation Change",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)
map("state",add=TRUE,col="gray")

testsfc = list(x=lon[lonidx],y=lat[latidx],z=ensmean)
surface(testsfc,type="I",main="GCW units overlaid with RCP8.5 2036-2065 \n Ensemble Mean Precipitation Change",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)
map("state",add=TRUE,col="gray")

testsfc = list(x=lon[lonidx],y=lat[latidx],z=ensmax)
surface(testsfc,type="I",main="GCW units overlaid with RCP8.5 2036-2065 \n Ensemble Highest Precipitation Change",xlab="Longtitude",ylab="Latitude",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]])
plot(test,add=TRUE)
map("state",add=TRUE,col="gray")


dev.off()
