
filespr = system("ls /data2/3to5/I35/pr/EDQM/*.nc",intern=TRUE)
filestmax = system("ls /data2/3to5/I35/tasmax/EDQM/*.nc",intern=TRUE)
filestmin = system("ls /data2/3to5/I35/tasmin/EDQM/*.nc",intern=TRUE)

files = c(filespr,filestmax,filestmin)

scenin="rcp85"

#####################

library(ncdf4)
library(sp)
library(fields)
library(smacof)

filebreakdown = do.call(rbind,strsplit(files,"/",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,7],"_",fixed=TRUE))
filebreakdown3 = do.call(rbind,strsplit(filebreakdown2[,3],"-",fixed=TRUE))

filepath = paste("",filebreakdown[,2],filebreakdown[,3],filebreakdown[,4],filebreakdown[,5],filebreakdown[,6],sep="/")
filebreakdown4 = data.frame(filepath,filebreakdown2[,1:2],filebreakdown3,filebreakdown2[,4:5])

GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown4$GCM = rep(GCM,3)
filebreakdown4$obs = rep(c("Daymet","Livneh","PRISM"),27)
filebreakdown4 = filebreakdown4[,-4]
names(filebreakdown4)=c("filepath","var","tempres","DS","code","scen","experiment","GCM","obs")
filebreakdown4$filename = filebreakdown[,7]

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
dates = dates[-which(substr(dates,6,10)=="02-29")]

earlyidx = which(as.numeric(substr(dates,1,4))>=2006 & as.numeric(substr(dates,1,4))<=2015) 
lateidx = which(as.numeric(substr(dates,1,4))>=2041 & as.numeric(substr(dates,1,4))<=2070)

#####################
# Gather all data - GCM netcdfs

filebreakdownin = subset(filebreakdown4,scen==scenin)

for(f in 1:nrow(filebreakdownin)){
  ptm1 = proc.time()
  message("Data grab began at ",Sys.time()," for file: ",filebreakdownin$filename[f])
  
  var = filebreakdownin$var[f]
  
  test = nc_open(paste(filebreakdownin$filepath[f],filebreakdownin$filename[f],sep="/"))
  tempdataearly = ncvar_get(test,as.character(var),start=c(1,1,earlyidx[1]),count=c(-1,-1,length(earlyidx)))
  tempdatalate = ncvar_get(test,as.character(var),start=c(1,1,lateidx[1]),count=c(-1,-1,length(lateidx)))
  if(f==1){
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
  }
  nc_close(test)
  
  tempearly = template = array(NA,dim=c(length(lon),length(lat),12))
  
  datese = dates[earlyidx]
  datesl = dates[lateidx]
  for(m in 1:12){
    etmpidx = which(as.numeric(substr(datese,6,7))==m)
    ltmpidx = which(as.numeric(substr(datesl,6,7))==m)
    tempearly[,,m] = apply(tempdataearly[,,etmpidx],c(1,2),mean,na.rm=TRUE)
    template[,,m] = apply(tempdatalate[,,ltmpidx],c(1,2),mean,na.rm=TRUE)
  }
  if(f<=9){
    tempearly=tempearly*86400
    template=template*86400
  } 
  if(f>9 & f<=18){
    tempearly=tempearly-273.15
    template=template-273.15
  } 
  if(f>18 & f<=27){
    tempearly=tempearly-273.15
    template=template-273.15
  } 
  
  dimX <- ncdim_def("lon","degrees_east",lon)
  dimY <- ncdim_def("lat","degrees_north",lat)
  dimT <- ncdim_def("time","months",1:12)
  
  var1 <- ncvar_def(as.character(filebreakdownin$var[f]),"",dim=list(dimX,dimY,dimT),missval=999999)
  nc <- nc_create(paste("/home/woot0002/",filebreakdownin$var[f],"_monthlyclimo_",filebreakdownin$DS[f],"_",scenin,"_",filebreakdownin$GCM[f],"_",filebreakdownin$obs[f],"_2006-2015.nc",sep=""),var1)
  ncvar_put(nc,var1,tempearly)
  nc_close(nc)
  
  var1 <- ncvar_def(as.character(filebreakdownin$var[f]),"",dim=list(dimX,dimY,dimT),missval=999999)
  nc <- nc_create(paste("/home/woot0002/",filebreakdownin$var[f],"_monthlyclimo_",filebreakdownin$DS[f],"_",scenin,"_",filebreakdownin$GCM[f],"_",filebreakdownin$obs[f],"_2041-2070.nc",sep=""),var1)
  ncvar_put(nc,var1,template)
  nc_close(nc)
  
  message("Climo Calc finished at ",Sys.time())
}
