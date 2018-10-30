library(ncdf4)
library(maps)
library(fields)
library(sp)
source("/home/woot0002/scripts/analysisfunctions.R")

Dfile = "/home/woot0002/3to5/pr_day_daymet_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
Lfile = "/home/woot0002/3to5/pr_day_livneh_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
Pfile = "/home/woot0002/3to5/pr_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"

DQDMfile = "/data2/3to5/I35/pr/EDQM/pr_day_I35prp1-QDM-A10D01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
LQDMfile = "/data2/3to5/I35/pr/EDQM/pr_day_I35prp1-QDM-A10L01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
PQDMfile = "/data2/3to5/I35/pr/EDQM/pr_day_I35prp1-QDM-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"

nctest = nc_open(DQDMfile)
prDQDM = ncvar_get(nctest,"pr")
lon = ncvar_get(nctest,"lon")
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
datesidx = which(substr(dates,6,10)=="02-29")
if(dim(prDQDM)[3]<length(dates)){
  dates=dates[-datesidx]
}

years = 1981:2005
prDQDM_yearly = array(NA,dim=c(length(lon),length(lat),length(years)))

for(y in 1:length(years)){
  yearidx = which(as.numeric(substr(dates,1,4))==years[y])
  tmp = ifelse(prDQDM[,,yearidx]<(1/86400),NA,1)
  if(y==1) mask = ifelse(is.na(prDQDM[,,1])==TRUE,0,1)
  tmp2 = apply(tmp,c(1,2),sum,na.rm=TRUE)
  prDQDM_yearly[,,y] = ifelse(mask==1,tmp2,NA)
}

#######

nctest = nc_open(LQDMfile)
prLQDM = ncvar_get(nctest,"pr")
nc_close(nctest)

prLQDM_yearly = array(NA,dim=c(length(lon),length(lat),length(years)))

for(y in 1:length(years)){
  yearidx = which(as.numeric(substr(dates,1,4))==years[y])
  tmp = ifelse(prLQDM[,,yearidx]<(1/86400),NA,1)
  if(y==1) mask = ifelse(is.na(prLQDM[,,1])==TRUE,0,1)
  tmp2 = apply(tmp,c(1,2),sum,na.rm=TRUE)
  prLQDM_yearly[,,y] = ifelse(mask==1,tmp2,NA)
}

#######

nctest = nc_open(PQDMfile)
prPQDM = ncvar_get(nctest,"pr")
nc_close(nctest)

prPQDM_yearly = array(NA,dim=c(length(lon),length(lat),length(years)))

for(y in 1:length(years)){
  yearidx = which(as.numeric(substr(dates,1,4))==years[y])
  tmp = ifelse(prPQDM[,,yearidx]<(1/86400),NA,1)
  if(y==1) mask = ifelse(is.na(prPQDM[,,1])==TRUE,0,1)
  tmp2 = apply(tmp,c(1,2),sum,na.rm=TRUE)
  prPQDM_yearly[,,y] = ifelse(mask==1,tmp2,NA)
}

rm(prDQDM)
rm(prLQDM)
rm(prPQDM)
gc()

########

nctest = nc_open(Dfile)
prD = ncvar_get(nctest,"pr")
nc_close(nctest)

if(dim(prD)[3]>length(dates)){
  prD = prD[,,-datesidx]
}

prD_yearly = array(NA,dim=c(length(lon),length(lat),length(years)))
for(y in 1:length(years)){
  yearidx = which(as.numeric(substr(dates,1,4))==years[y])
  tmp = ifelse(prD[,,yearidx]<(1/86400),NA,1)
  if(y==1) mask = ifelse(is.na(prD[,,1])==TRUE,0,1)
  tmp2 = apply(tmp,c(1,2),sum,na.rm=TRUE)
  prD_yearly[,,y] = ifelse(mask==1,tmp2,NA)
}

########

nctest = nc_open(Lfile)
prL = ncvar_get(nctest,"pr")
nc_close(nctest)

if(dim(prL)[3]>length(dates)){
  prL = prL[,,-datesidx]
}

prL_yearly = array(NA,dim=c(length(lon),length(lat),length(years)))
for(y in 1:length(years)){
  yearidx = which(as.numeric(substr(dates,1,4))==years[y])
  tmp = ifelse(prL[,,yearidx]<(1/86400),NA,1)
  if(y==1) mask = ifelse(is.na(prL[,,1])==TRUE,0,1)
  tmp2 = apply(tmp,c(1,2),sum,na.rm=TRUE)
  prL_yearly[,,y] = ifelse(mask==1,tmp2,NA)
}

########

nctest = nc_open(Pfile)
prP = ncvar_get(nctest,"pr")
nc_close(nctest)

if(dim(prP)[3]>length(dates)){
  prP = prP[,,-datesidx]
}

prP_yearly = array(NA,dim=c(length(lon),length(lat),length(years)))
for(y in 1:length(years)){
  yearidx = which(as.numeric(substr(dates,1,4))==years[y])
  tmp = ifelse(prP[,,yearidx]<(1/86400),NA,1)
  if(y==1) mask = ifelse(is.na(prP[,,1])==TRUE,0,1)
  tmp2 = apply(tmp,c(1,2),sum,na.rm=TRUE)
  prP_yearly[,,y] = ifelse(mask==1,tmp2,NA)
}

rm(prD)
rm(prL)
rm(prP)
gc()

##########

prDQDM_climo = apply(prDQDM_yearly,c(1,2),mean,na.rm=TRUE)
prLQDM_climo = apply(prLQDM_yearly,c(1,2),mean,na.rm=TRUE)
prPQDM_climo = apply(prPQDM_yearly,c(1,2),mean,na.rm=TRUE)
prD_climo = apply(prD_yearly,c(1,2),mean,na.rm=TRUE)
prL_climo = apply(prL_yearly,c(1,2),mean,na.rm=TRUE)
prP_climo = apply(prP_yearly,c(1,2),mean,na.rm=TRUE)

rm(prDQDM_yearly)
rm(prLQDM_yearly)
rm(prPQDM_yearly)
rm(prD_yearly)
rm(prL_yearly)
rm(prP_yearly)
gc()

##########

prD_climo = ifelse(is.na(prL_climo)==FALSE,prD_climo,NA)
prP_climo = ifelse(is.na(prL_climo)==FALSE,prP_climo,NA)

prDQDM_climo = ifelse(is.na(prL_climo)==FALSE,prDQDM_climo,NA)
prLQDM_climo = ifelse(is.na(prL_climo)==FALSE,prLQDM_climo,NA)
prPQDM_climo = ifelse(is.na(prL_climo)==FALSE,prPQDM_climo,NA)

diffcolorbar = colorramp(c(prD_climo,prL_climo,prP_climo,prDQDM_climo,prLQDM_climo,prPQDM_climo),colorchoice="whitetogreen",Blimit=20,type="raw")

testsfc1 = list(x=lon-360,y=lat,z=prD_climo)
surface(testsfc1,type="I",main="Daymet",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=prL_climo)
surface(testsfc1,type="I",main="Livneh",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=prP_climo)
surface(testsfc1,type="I",main="PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

##############

rnnmm = c(prD_climo,prL_climo,prP_climo,prDQDM_climo,prLQDM_climo,prPQDM_climo)
modname = rep(c("Daymet","Livneh","PRISM","CCSM4_QDM_Daymet","CCSM4_QDM_Livneh","CCSM4_QDM_PRISM"),each=length(prD_climo))
moddat = data.frame(modname,rnnmm)

moddat$plotorder = NA
moddat$plotorder[which(moddat$modname=="Daymet")]=1
moddat$plotorder[which(moddat$modname=="Livneh")]=2
moddat$plotorder[which(moddat$modname=="PRISM")]=3

moddat$plotorder[which(moddat$modname=="CCSM4_QDM_Daymet")]=4
moddat$plotorder[which(moddat$modname=="CCSM4_QDM_Livneh")]=5
moddat$plotorder[which(moddat$modname=="CCSM4_QDM_PRISM")]=6


boxplot(rnnmm~plotorder,data=moddat,main="rnnmm Obs, GCM, and DS Values across I35 region grid points",ylab="days",at=c(1,2,3,5,6,7),xaxt="n",col=c("red","purple","magenta","red","purple","magenta"))
text(x =  c(1,2,3,5,6,7), y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = unique(moddat$modname), xpd = TRUE,cex=0.75)





