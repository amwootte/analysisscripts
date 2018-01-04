#################################
#
# Create MACA climatologies

scenin="rcp85"

library(ncdf4)

filepath = paste("/data1/MACAv2LIVNEH_sub/",scenin,sep="")

##########
# Get 3^5 lat and lon

test = nc_open("/home/woot0002/RRdomainmask.nc")
targetlon = ncvar_get(test,"lon")
targetlat = ncvar_get(test,"lat")
nc_close(test)

#########
# Early period 2006-2015

earlydate = seq(as.Date("2006-01-01"),as.Date("2025-12-31"),by="day")
pullidx = which(as.numeric(substr(earlydate,1,4))>=2006 & as.numeric(substr(earlydate,1,4))<=2015)
pulldate = earlydate[pullidx]

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_pr_CCSM4_r6i1p1_rcp85_2006_2025_CONUS_daily.nc",sep=""))
lat = ncvar_get(test1,"lat")
lon = ncvar_get(test1,"lon")-360
lonidx = which(lon> -110 & lon< -89)
latidx = which(lat> 25 & lat< 41)
prCCSM4 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],pullidx[1]),count=c(length(lonidx),length(latidx),length(pullidx))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_pr_MIROC5_r1i1p1_rcp85_2006_2025_CONUS_daily.nc",sep=""))
prMIROC5 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],pullidx[1]),count=c(length(lonidx),length(latidx),length(pullidx))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmax_CCSM4_r6i1p1_rcp85_2006_2025_CONUS_daily.nc",sep=""))
tmaxCCSM4 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx[1]),count=c(length(lonidx),length(latidx),length(pullidx))) # in K
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmax_MIROC5_r1i1p1_rcp85_2006_2025_CONUS_daily.nc",sep=""))
tmaxMIROC5 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx[1]),count=c(length(lonidx),length(latidx),length(pullidx))) # in K
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmin_CCSM4_r6i1p1_rcp85_2006_2025_CONUS_daily.nc",sep=""))
tminCCSM4 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx[1]),count=c(length(lonidx),length(latidx),length(pullidx))) # in K
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmin_MIROC5_r1i1p1_rcp85_2006_2025_CONUS_daily.nc",sep=""))
tminMIROC5 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx[1]),count=c(length(lonidx),length(latidx),length(pullidx))) # in K
nc_close(test1)

#########
# Calculate climos for everything

pravgCCSM4 = pravgMIROC5 = tmaxavgCCSM4 = tmaxavgMIROC5 = tminavgCCSM4 = tminavgMIROC5 = array(NA,dim=c(length(lonidx),length(latidx),12))
for(i in 1:12){
  dateidx = which(as.numeric(substr(pulldate,6,7))==i)
  pravgCCSM4[,,i] = apply(prCCSM4[,,dateidx],c(1,2),mean,na.rm=TRUE)
  pravgMIROC5[,,i] = apply(prMIROC5[,,dateidx],c(1,2),mean,na.rm=TRUE)
  tmaxavgCCSM4[,,i] = apply(tmaxCCSM4[,,dateidx],c(1,2),mean,na.rm=TRUE)
  tmaxavgMIROC5[,,i] = apply(tmaxMIROC5[,,dateidx],c(1,2),mean,na.rm=TRUE)
  tminavgCCSM4[,,i] = apply(tminCCSM4[,,dateidx],c(1,2),mean,na.rm=TRUE)
  tminavgMIROC5[,,i] = apply(tminMIROC5[,,dateidx],c(1,2),mean,na.rm=TRUE)
}
rm(prCCSM4) ; rm(prMIROC5) ; rm(tmaxCCSM4) ; rm(tmaxMIROC5) ; rm(tminCCSM4) ; rm(tminMIROC5)

###########
# Regrid to 3^5 grid

library(fields)
library(sp)

newrownum = length(targetlon) # number of cells in x direction
newcolnum = length(targetlat)  # number of cells in y direction

prCCSM4out = prMIROC5out = tmaxCCSM4out = tmaxMIROC5out = tminCCSM4out = tminMIROC5out = array(NA,dim=c(newrownum,newcolnum,12))

for(j in 1:12){
  ptm = proc.time()
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=pravgCCSM4[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  prCCSM4out[,,j] = newvar
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=pravgMIROC5[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  prMIROC5out[,,j] = newvar
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tmaxavgCCSM4[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tmaxCCSM4out[,,j] = newvar
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tmaxavgMIROC5[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tmaxMIROC5out[,,j] = newvar
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tminavgCCSM4[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tminCCSM4out[,,j] = newvar
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tminavgMIROC5[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tminMIROC5out[,,j] = newvar
  
  message("Regridded Time slice ",j," / ",12)
  ptmend = proc.time()-ptm
  message("Time = ",ptmend[3]," secs")
  
}

library(maps)
testsfc = list(x=targetlon,y=targetlat,z=tminCCSM4out[,,7])
surface(testsfc,type="I")
map("state",add=TRUE)

#############
# Create netcdfs

dimX <- ncdim_def("lon","degrees_east",targetlon)
dimY <- ncdim_def("lat","degrees_north",targetlat)
dimT <- ncdim_def("time","months",1:12)

####
# CCSM4

var1 <- ncvar_def("pr","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_pr_",scenin,"CCSM4_2006-2015.nc",sep=""),var1)
ncvar_put(nc,var1,prCCSM4out)
nc_close(nc)

var1 <- ncvar_def("tasmax","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmax_",scenin,"CCSM4_2006-2015.nc",sep=""),var1)
ncvar_put(nc,var1,tmaxCCSM4out)
nc_close(nc)

var1 <- ncvar_def("tasmin","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmin_",scenin,"CCSM4_2006-2015.nc",sep=""),var1)
ncvar_put(nc,var1,tminCCSM4out)
nc_close(nc)

####
# MIROC5

var1 <- ncvar_def("pr","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_pr_",scenin,"MIROC5_2006-2015.nc",sep=""),var1)
ncvar_put(nc,var1,prMIROC5out)
nc_close(nc)

var1 <- ncvar_def("tasmax","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmax_",scenin,"MIROC5_2006-2015.nc",sep=""),var1)
ncvar_put(nc,var1,tmaxMIROC5out)
nc_close(nc)

var1 <- ncvar_def("tasmin","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmin_",scenin,"MIROC5_2006-2015.nc",sep=""),var1)
ncvar_put(nc,var1,tminMIROC5out)
nc_close(nc)













