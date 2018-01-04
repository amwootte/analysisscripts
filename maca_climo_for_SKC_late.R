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
# late period 2041-2070

dates1 = seq(as.Date("2026-01-01"),as.Date("2045-12-31"),by="day")
pullidx1 = which(as.numeric(substr(dates1,1,4))>=2041)
pulldate1 = dates1[pullidx1]

dates2 = seq(as.Date("2066-01-01"),as.Date("2085-12-31"),by="day")
pullidx2 = which(as.numeric(substr(dates2,1,4))<=2070)
pulldate2 = dates2[pullidx2]

dates = seq(as.Date("2041-01-01"),as.Date("2070-12-31"),by="day")

########################################################
##############
# CCSM4 precip

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_pr_CCSM4_r6i1p1_rcp85_2026_2045_CONUS_daily.nc",sep=""))
lat = ncvar_get(test1,"lat")
lon = ncvar_get(test1,"lon")-360
lonidx = which(lon> -110 & lon< -89)
latidx = which(lat> 25 & lat< 41)
prCCSM4_p1 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],pullidx1[1]),count=c(length(lonidx),length(latidx),length(pullidx1))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_pr_CCSM4_r6i1p1_rcp85_2046_2065_CONUS_daily.nc",sep=""))
prCCSM4_p2 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1)) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_pr_CCSM4_r6i1p1_rcp85_2066_2085_CONUS_daily.nc",sep=""))
prCCSM4_p3 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],pullidx2[1]),count=c(length(lonidx),length(latidx),length(pullidx2))) # in mm
nc_close(test1)

prCCSM4 = abind(prCCSM4_p1,prCCSM4_p2,along=3)
rm(prCCSM4_p1); rm(prCCSM4_p2); gc()
prCCSM4 = abind(prCCSM4,prCCSM4_p3,along=3)
rm(prCCSM4_p3)
gc()

###
# Calculate climo
pravgCCSM4 = array(NA,dim=c(length(lonidx),length(latidx),12))
for(i in 1:12){
  dateidx = which(as.numeric(substr(dates,6,7))==i)
  pravgCCSM4[,,i] = apply(prCCSM4[,,dateidx],c(1,2),mean,na.rm=TRUE)
}
rm(prCCSM4) ; gc()

####
# regrid to 3^5

newrownum = length(targetlon) # number of cells in x direction
newcolnum = length(targetlat)  # number of cells in y direction

prCCSM4out = array(NA,dim=c(newrownum,newcolnum,12))

for(j in 1:12){
  ptm = proc.time()
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=pravgCCSM4[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  prCCSM4out[,,j] = newvar

  message("Regridded Time slice ",j," / ",12)
  ptmend = proc.time()-ptm
  message("Time = ",ptmend[3]," secs")
  
}

rm(pravgCCSM4); gc()

##############
# CCSM4 tasmax

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmax_CCSM4_r6i1p1_rcp85_2026_2045_CONUS_daily.nc",sep=""))
tmaxCCSM4_p1 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx1[1]),count=c(length(lonidx),length(latidx),length(pullidx1))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmax_CCSM4_r6i1p1_rcp85_2046_2065_CONUS_daily.nc",sep=""))
tmaxCCSM4_p2 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1)) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmax_CCSM4_r6i1p1_rcp85_2066_2085_CONUS_daily.nc",sep=""))
tmaxCCSM4_p3 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx2[1]),count=c(length(lonidx),length(latidx),length(pullidx2))) # in mm
nc_close(test1)

tmaxCCSM4 = abind(tmaxCCSM4_p1,tmaxCCSM4_p2,along=3)
rm(tmaxCCSM4_p1); rm(tmaxCCSM4_p2); gc()
tmaxCCSM4 = abind(tmaxCCSM4,tmaxCCSM4_p3,along=3)
rm(tmaxCCSM4_p3)
gc()

###
# Calculate climo
tmaxavgCCSM4 = array(NA,dim=c(length(lonidx),length(latidx),12))
for(i in 1:12){
  dateidx = which(as.numeric(substr(dates,6,7))==i)
  tmaxavgCCSM4[,,i] = apply(tmaxCCSM4[,,dateidx],c(1,2),mean,na.rm=TRUE)
}
rm(tmaxCCSM4) ; gc()

####
# regrid to 3^5

tmaxCCSM4out = array(NA,dim=c(newrownum,newcolnum,12))

for(j in 1:12){
  ptm = proc.time()
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tmaxavgCCSM4[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tmaxCCSM4out[,,j] = newvar
  
  message("Regridded Time slice ",j," / ",12)
  ptmend = proc.time()-ptm
  message("Time = ",ptmend[3]," secs")
  
}

rm(tmaxavgCCSM4); gc()

##############
# CCSM4 tasmin

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmin_CCSM4_r6i1p1_rcp85_2026_2045_CONUS_daily.nc",sep=""))
tminCCSM4_p1 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx1[1]),count=c(length(lonidx),length(latidx),length(pullidx1))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmin_CCSM4_r6i1p1_rcp85_2046_2065_CONUS_daily.nc",sep=""))
tminCCSM4_p2 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1)) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/CCSM4/macav2livneh_tasmin_CCSM4_r6i1p1_rcp85_2066_2085_CONUS_daily.nc",sep=""))
tminCCSM4_p3 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx2[1]),count=c(length(lonidx),length(latidx),length(pullidx2))) # in mm
nc_close(test1)

tminCCSM4 = abind(tminCCSM4_p1,tminCCSM4_p2,along=3)
rm(tminCCSM4_p1); rm(tminCCSM4_p2); gc()
tminCCSM4 = abind(tminCCSM4,tminCCSM4_p3,along=3)
rm(tminCCSM4_p3)
gc()

###
# Calculate climo
tminavgCCSM4 = array(NA,dim=c(length(lonidx),length(latidx),12))
for(i in 1:12){
  dateidx = which(as.numeric(substr(dates,6,7))==i)
  tminavgCCSM4[,,i] = apply(tminCCSM4[,,dateidx],c(1,2),mean,na.rm=TRUE)
}
rm(tminCCSM4) ; gc()

####
# regrid to 3^5

tminCCSM4out = array(NA,dim=c(newrownum,newcolnum,12))

for(j in 1:12){
  ptm = proc.time()
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tminavgCCSM4[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tminCCSM4out[,,j] = newvar
  
  message("Regridded Time slice ",j," / ",12)
  ptmend = proc.time()-ptm
  message("Time = ",ptmend[3]," secs")
  
}

rm(tminavgCCSM4); gc()


########################################################
##############
# MIROC5 precip

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_pr_MIROC5_r1i1p1_rcp85_2026_2045_CONUS_daily.nc",sep=""))
prMIROC5_p1 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],pullidx1[1]),count=c(length(lonidx),length(latidx),length(pullidx1))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_pr_MIROC5_r1i1p1_rcp85_2046_2065_CONUS_daily.nc",sep=""))
prMIROC5_p2 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1)) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_pr_MIROC5_r1i1p1_rcp85_2066_2085_CONUS_daily.nc",sep=""))
prMIROC5_p3 = ncvar_get(test1,"precipitation",start=c(lonidx[1],latidx[1],pullidx2[1]),count=c(length(lonidx),length(latidx),length(pullidx2))) # in mm
nc_close(test1)

prMIROC5 = abind(prMIROC5_p1,prMIROC5_p2,along=3)
rm(prMIROC5_p1); rm(prMIROC5_p2); gc()
prMIROC5 = abind(prMIROC5,prMIROC5_p3,along=3)
rm(prMIROC5_p3)
gc()

###
# Calculate climo
pravgMIROC5 = array(NA,dim=c(length(lonidx),length(latidx),12))
for(i in 1:12){
  dateidx = which(as.numeric(substr(dates,6,7))==i)
  pravgMIROC5[,,i] = apply(prMIROC5[,,dateidx],c(1,2),mean,na.rm=TRUE)
}
rm(prMIROC5) ; gc()

####
# regrid to 3^5

newrownum = length(targetlon) # number of cells in x direction
newcolnum = length(targetlat)  # number of cells in y direction

prMIROC5out = array(NA,dim=c(newrownum,newcolnum,12))

for(j in 1:12){
  ptm = proc.time()
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=pravgMIROC5[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  prMIROC5out[,,j] = newvar
  
  message("Regridded Time slice ",j," / ",12)
  ptmend = proc.time()-ptm
  message("Time = ",ptmend[3]," secs")
  
}

rm(pravgMIROC5); gc()

##############
# MIROC5 tasmax

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmax_MIROC5_r1i1p1_rcp85_2026_2045_CONUS_daily.nc",sep=""))
tmaxMIROC5_p1 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx1[1]),count=c(length(lonidx),length(latidx),length(pullidx1))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmax_MIROC5_r1i1p1_rcp85_2046_2065_CONUS_daily.nc",sep=""))
tmaxMIROC5_p2 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1)) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmax_MIROC5_r1i1p1_rcp85_2066_2085_CONUS_daily.nc",sep=""))
tmaxMIROC5_p3 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx2[1]),count=c(length(lonidx),length(latidx),length(pullidx2))) # in mm
nc_close(test1)

tmaxMIROC5 = abind(tmaxMIROC5_p1,tmaxMIROC5_p2,along=3)
rm(tmaxMIROC5_p1); rm(tmaxMIROC5_p2); gc()
tmaxMIROC5 = abind(tmaxMIROC5,tmaxMIROC5_p3,along=3)
rm(tmaxMIROC5_p3)
gc()

###
# Calculate climo
tmaxavgMIROC5 = array(NA,dim=c(length(lonidx),length(latidx),12))
for(i in 1:12){
  dateidx = which(as.numeric(substr(dates,6,7))==i)
  tmaxavgMIROC5[,,i] = apply(tmaxMIROC5[,,dateidx],c(1,2),mean,na.rm=TRUE)
}
rm(tmaxMIROC5) ; gc()

####
# regrid to 3^5

tmaxMIROC5out = array(NA,dim=c(newrownum,newcolnum,12))

for(j in 1:12){
  ptm = proc.time()
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tmaxavgMIROC5[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tmaxMIROC5out[,,j] = newvar
  
  message("Regridded Time slice ",j," / ",12)
  ptmend = proc.time()-ptm
  message("Time = ",ptmend[3]," secs")
  
}

rm(tmaxavgMIROC5); gc()

##############
# MIROC5 tasmin

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmin_MIROC5_r1i1p1_rcp85_2026_2045_CONUS_daily.nc",sep=""))
tminMIROC5_p1 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx1[1]),count=c(length(lonidx),length(latidx),length(pullidx1))) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmin_MIROC5_r1i1p1_rcp85_2046_2065_CONUS_daily.nc",sep=""))
tminMIROC5_p2 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],1),count=c(length(lonidx),length(latidx),-1)) # in mm
nc_close(test1)

test1 = nc_open(paste(filepath,"/MIROC5/macav2livneh_tasmin_MIROC5_r1i1p1_rcp85_2066_2085_CONUS_daily.nc",sep=""))
tminMIROC5_p3 = ncvar_get(test1,"air_temperature",start=c(lonidx[1],latidx[1],pullidx2[1]),count=c(length(lonidx),length(latidx),length(pullidx2))) # in mm
nc_close(test1)

tminMIROC5 = abind(tminMIROC5_p1,tminMIROC5_p2,along=3)
rm(tminMIROC5_p1); rm(tminMIROC5_p2); gc()
tminMIROC5 = abind(tminMIROC5,tminMIROC5_p3,along=3)
rm(tminMIROC5_p3)
gc()

###
# Calculate climo
tminavgMIROC5 = array(NA,dim=c(length(lonidx),length(latidx),12))
for(i in 1:12){
  dateidx = which(as.numeric(substr(dates,6,7))==i)
  tminavgMIROC5[,,i] = apply(tminMIROC5[,,dateidx],c(1,2),mean,na.rm=TRUE)
}
rm(tminMIROC5) ; gc()

####
# regrid to 3^5

tminMIROC5out = array(NA,dim=c(newrownum,newcolnum,12))

for(j in 1:12){
  ptm = proc.time()
  
  test2 = interp.surface.grid(list(x=lon[lonidx],y=lat[latidx],z=tminavgMIROC5[,,j]),grid.list=list(x=targetlon,y=targetlat))
  newvar = matrix(test2[[3]],nrow=newrownum,ncol=newcolnum)
  tminMIROC5out[,,j] = newvar
  
  message("Regridded Time slice ",j," / ",12)
  ptmend = proc.time()-ptm
  message("Time = ",ptmend[3]," secs")
  
}

rm(tminavgMIROC5); gc()

#############
# Create netcdfs

dimX <- ncdim_def("lon","degrees_east",targetlon)
dimY <- ncdim_def("lat","degrees_north",targetlat)
dimT <- ncdim_def("time","months",1:12)

####
# CCSM4

var1 <- ncvar_def("pr","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_pr_",scenin,"CCSM4_2041-2070.nc",sep=""),var1)
ncvar_put(nc,var1,prCCSM4out)
nc_close(nc)

var1 <- ncvar_def("tasmax","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmax_",scenin,"CCSM4_2041-2070.nc",sep=""),var1)
ncvar_put(nc,var1,tmaxCCSM4out)
nc_close(nc)

var1 <- ncvar_def("tasmin","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmin_",scenin,"CCSM4_2041-2070.nc",sep=""),var1)
ncvar_put(nc,var1,tminCCSM4out)
nc_close(nc)

####
# MIROC5

var1 <- ncvar_def("pr","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_pr_",scenin,"MIROC5_2041-2070.nc",sep=""),var1)
ncvar_put(nc,var1,prMIROC5out)
nc_close(nc)

var1 <- ncvar_def("tasmax","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmax_",scenin,"MIROC5_2041-2070.nc",sep=""),var1)
ncvar_put(nc,var1,tmaxMIROC5out)
nc_close(nc)

var1 <- ncvar_def("tasmin","",dim=list(dimX,dimY,dimT),missval=999999)
nc <- nc_create(paste("/home/woot0002/maca_monthlyclimo_tasmin_",scenin,"MIROC5_2041-2070.nc",sep=""),var1)
ncvar_put(nc,var1,tminMIROC5out)
nc_close(nc)
