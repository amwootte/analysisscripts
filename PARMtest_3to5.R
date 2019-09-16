
source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

setwd("/home/woot0002/")

PARMpastfile = "tasmin_day_I35tndetrp1-PARM-Bd10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
PARMfuturefile = "tasmin_day_I35tndetrp1-PARM-Bd18P01K00_rcp85_r6i1p1_I35Land_20060101-20991231.nc"

EDQMpastfile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
EDQMfuturefile = "/data2/3to5/I35/tasmin/EDQM/tasmin_day_I35tnp1-EDQM-A18P01K00_rcp85_r6i1p1_I35Land_20060101-20991231.nc"

PRISMpastfile = "/data2/3to5/I35/tasmin/PRISM/tasmin_day_prism_historical_r0i0p0_SCCSC0p1_19810101-20051231.nc"
Deltafuturefile = "/data2/3to5/I35/tasmin/DeltaSD/tasmin_day_I35tndetrp1-DeltaSD-A18P01K00_rcp85_r6i1p1_I35Land_20060101-20991231.nc"


###
# PARM data

test = nc_open(PARMpastfile)
PARMtmaxhist = ncvar_get(test,"tasmin")
lat = ncvar_get(test,"lat")
lon = ncvar_get(test,"lon")
nc_close(test)

test = nc_open(PARMfuturefile)
PARMtmaxfut = ncvar_get(test,"tasmin")
nc_close(test)

histdates= seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
futdates= seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

if(length(histdates)>dim(PARMtmaxhist)[3]){
  histdates2 = histdates[-which(substr(histdates,6,10)=="02-29")]
}

if(length(futdates)>dim(PARMtmaxfut)[3]){
  futdates2 = futdates[-which(substr(futdates,6,10)=="02-29")]
}

###
# EDQM data

test = nc_open(EDQMpastfile)
EDQMtmaxhist = ncvar_get(test,"tasmin")
nc_close(test)

test = nc_open(EDQMfuturefile)
EDQMtmaxfut = ncvar_get(test,"tasmin")
nc_close(test)

###
# PRISM and Delta data

test = nc_open(PRISMpastfile)
PRISMtmaxhist = ncvar_get(test,"tasmin")
nc_close(test)

test = nc_open(Deltafuturefile)
Deltatmaxfut = ncvar_get(test,"tasmin")
nc_close(test)

####

futstartyears = c(2006,2040,2070)
futendyears = c(2039,2069,2099)

EDQMhistlist = PARMhistlist = PRISMhistlist = list()
EDQMlist_early = PARMlist_early = Deltalist_early = EDQMlist_mid = PARMlist_mid = Deltalist_mid = EDQMlist_late = PARMlist_late = Deltalist_late = list()
vars = c("tmax_min","tmax_05","tmax_50","tmax_mean","tmax_95","tmax_max")

mask = ifelse(is.na(EDQMtmaxhist[,,1])==FALSE,1,0)
earlyidx = which(as.numeric(substr(futdates2,1,4))>=2006 & as.numeric(substr(futdates2,1,4))<=2039)
mididx = which(as.numeric(substr(futdates2,1,4))>=2040 & as.numeric(substr(futdates2,1,4))<=2069)
lateidx = which(as.numeric(substr(futdates2,1,4))>=2070 & as.numeric(substr(futdates2,1,4))<=2099)
              
for(i in 1:6){
  
  if(vars[i]=="tmax_min"){
    EDQMhistlist[[i]]=apply(EDQMtmaxhist,c(1,2),min,na.rm=TRUE)
    PARMhistlist[[i]]=apply(PARMtmaxhist,c(1,2),min,na.rm=TRUE)
    PRISMhistlist[[i]]=apply(PRISMtmaxhist,c(1,2),min,na.rm=TRUE)
    
    EDQMlist_early[[i]]=apply(EDQMtmaxfut[,,earlyidx],c(1,2),min,na.rm=TRUE)
    PARMlist_early[[i]]=apply(PARMtmaxfut[,,earlyidx],c(1,2),min,na.rm=TRUE)
    Deltalist_early[[i]]=apply(Deltatmaxfut[,,earlyidx],c(1,2),min,na.rm=TRUE)
    
    EDQMlist_mid[[i]]=apply(EDQMtmaxfut[,,mididx],c(1,2),min,na.rm=TRUE)
    PARMlist_mid[[i]]=apply(PARMtmaxfut[,,mididx],c(1,2),min,na.rm=TRUE)
    Deltalist_mid[[i]]=apply(Deltatmaxfut[,,mididx],c(1,2),min,na.rm=TRUE)
    
    EDQMlist_late[[i]]=apply(EDQMtmaxfut[,,lateidx],c(1,2),min,na.rm=TRUE)
    PARMlist_late[[i]]=apply(PARMtmaxfut[,,lateidx],c(1,2),min,na.rm=TRUE)
    Deltalist_late[[i]]=apply(Deltatmaxfut[,,lateidx],c(1,2),min,na.rm=TRUE)
  }
  
  if(vars[i]=="tmax_mean"){
    EDQMhistlist[[i]]=apply(EDQMtmaxhist,c(1,2),mean,na.rm=TRUE)
    PARMhistlist[[i]]=apply(PARMtmaxhist,c(1,2),mean,na.rm=TRUE)
    PRISMhistlist[[i]]=apply(PRISMtmaxhist,c(1,2),mean,na.rm=TRUE)
    
    EDQMlist_early[[i]]=apply(EDQMtmaxfut[,,earlyidx],c(1,2),mean,na.rm=TRUE)
    PARMlist_early[[i]]=apply(PARMtmaxfut[,,earlyidx],c(1,2),mean,na.rm=TRUE)
    Deltalist_early[[i]]=apply(Deltatmaxfut[,,earlyidx],c(1,2),mean,na.rm=TRUE)
    
    EDQMlist_mid[[i]]=apply(EDQMtmaxfut[,,mididx],c(1,2),mean,na.rm=TRUE)
    PARMlist_mid[[i]]=apply(PARMtmaxfut[,,mididx],c(1,2),mean,na.rm=TRUE)
    Deltalist_mid[[i]]=apply(Deltatmaxfut[,,mididx],c(1,2),mean,na.rm=TRUE)
    
    EDQMlist_late[[i]]=apply(EDQMtmaxfut[,,lateidx],c(1,2),mean,na.rm=TRUE)
    PARMlist_late[[i]]=apply(PARMtmaxfut[,,lateidx],c(1,2),mean,na.rm=TRUE)
    Deltalist_late[[i]]=apply(Deltatmaxfut[,,lateidx],c(1,2),mean,na.rm=TRUE)
  }
  
  if(vars[i]=="tmax_max"){
    EDQMhistlist[[i]]=apply(EDQMtmaxhist,c(1,2),max,na.rm=TRUE)
    PARMhistlist[[i]]=apply(PARMtmaxhist,c(1,2),max,na.rm=TRUE)
    PRISMhistlist[[i]]=apply(PRISMtmaxhist,c(1,2),max,na.rm=TRUE)
    
    EDQMlist_early[[i]]=apply(EDQMtmaxfut[,,earlyidx],c(1,2),max,na.rm=TRUE)
    PARMlist_early[[i]]=apply(PARMtmaxfut[,,earlyidx],c(1,2),max,na.rm=TRUE)
    Deltalist_early[[i]]=apply(Deltatmaxfut[,,earlyidx],c(1,2),max,na.rm=TRUE)
    
    EDQMlist_mid[[i]]=apply(EDQMtmaxfut[,,mididx],c(1,2),max,na.rm=TRUE)
    PARMlist_mid[[i]]=apply(PARMtmaxfut[,,mididx],c(1,2),max,na.rm=TRUE)
    Deltalist_mid[[i]]=apply(Deltatmaxfut[,,mididx],c(1,2),max,na.rm=TRUE)
  
    EDQMlist_late[[i]]=apply(EDQMtmaxfut[,,lateidx],c(1,2),max,na.rm=TRUE)
    PARMlist_late[[i]]=apply(PARMtmaxfut[,,lateidx],c(1,2),max,na.rm=TRUE)
    Deltalist_late[[i]]=apply(Deltatmaxfut[,,lateidx],c(1,2),max,na.rm=TRUE)
  }
  
  if(vars[i]=="tmax_05"){
    EDQMhistlist[[i]]=apply(EDQMtmaxhist,c(1,2),quantile,probs=0.05,na.rm=TRUE)
    PARMhistlist[[i]]=apply(PARMtmaxhist,c(1,2),quantile,probs=0.05,na.rm=TRUE)
    PRISMhistlist[[i]]=apply(PRISMtmaxhist,c(1,2),quantile,probs=0.05,na.rm=TRUE)
    
    EDQMlist_early[[i]]=apply(EDQMtmaxfut[,,earlyidx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    PARMlist_early[[i]]=apply(PARMtmaxfut[,,earlyidx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    Deltalist_early[[i]]=apply(Deltatmaxfut[,,earlyidx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    
    EDQMlist_mid[[i]]=apply(EDQMtmaxfut[,,mididx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    PARMlist_mid[[i]]=apply(PARMtmaxfut[,,mididx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    Deltalist_mid[[i]]=apply(Deltatmaxfut[,,mididx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    
    EDQMlist_late[[i]]=apply(EDQMtmaxfut[,,lateidx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    PARMlist_late[[i]]=apply(PARMtmaxfut[,,lateidx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
    Deltalist_late[[i]]=apply(Deltatmaxfut[,,lateidx],c(1,2),quantile,probs=0.05,na.rm=TRUE)
  }
  
  if(vars[i]=="tmax_50"){
    EDQMhistlist[[i]]=apply(EDQMtmaxhist,c(1,2),quantile,probs=0.5,na.rm=TRUE)
    PARMhistlist[[i]]=apply(PARMtmaxhist,c(1,2),quantile,probs=0.5,na.rm=TRUE)
    PRISMhistlist[[i]]=apply(PRISMtmaxhist,c(1,2),quantile,probs=0.5,na.rm=TRUE)
    
    EDQMlist_early[[i]]=apply(EDQMtmaxfut[,,earlyidx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    PARMlist_early[[i]]=apply(PARMtmaxfut[,,earlyidx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    Deltalist_early[[i]]=apply(Deltatmaxfut[,,earlyidx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    
    EDQMlist_mid[[i]]=apply(EDQMtmaxfut[,,mididx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    PARMlist_mid[[i]]=apply(PARMtmaxfut[,,mididx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    Deltalist_mid[[i]]=apply(Deltatmaxfut[,,mididx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    
    EDQMlist_late[[i]]=apply(EDQMtmaxfut[,,lateidx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    PARMlist_late[[i]]=apply(PARMtmaxfut[,,lateidx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
    Deltalist_late[[i]]=apply(Deltatmaxfut[,,lateidx],c(1,2),quantile,probs=0.5,na.rm=TRUE)
  }
  
  if(vars[i]=="tmax_95"){
    EDQMhistlist[[i]]=apply(EDQMtmaxhist,c(1,2),quantile,probs=0.95,na.rm=TRUE)
    PARMhistlist[[i]]=apply(PARMtmaxhist,c(1,2),quantile,probs=0.95,na.rm=TRUE)
    PRISMhistlist[[i]]=apply(PRISMtmaxhist,c(1,2),quantile,probs=0.95,na.rm=TRUE)
    
    EDQMlist_early[[i]]=apply(EDQMtmaxfut[,,earlyidx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    PARMlist_early[[i]]=apply(PARMtmaxfut[,,earlyidx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    Deltalist_early[[i]]=apply(Deltatmaxfut[,,earlyidx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    
    EDQMlist_mid[[i]]=apply(EDQMtmaxfut[,,mididx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    PARMlist_mid[[i]]=apply(PARMtmaxfut[,,mididx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    Deltalist_mid[[i]]=apply(Deltatmaxfut[,,mididx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    
    EDQMlist_late[[i]]=apply(EDQMtmaxfut[,,lateidx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    PARMlist_late[[i]]=apply(PARMtmaxfut[,,lateidx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
    Deltalist_late[[i]]=apply(Deltatmaxfut[,,lateidx],c(1,2),quantile,probs=0.95,na.rm=TRUE)
  }
  
  EDQMhistlist[[i]] = ifelse(mask==1,EDQMhistlist[[i]],NA)
  PARMhistlist[[i]] = ifelse(mask==1,PARMhistlist[[i]],NA)
  PRISMhistlist[[i]] = ifelse(mask==1,PRISMhistlist[[i]],NA)
  
  EDQMlist_early[[i]] = ifelse(mask==1,EDQMlist_early[[i]],NA)
  PARMlist_early[[i]] = ifelse(mask==1,PARMlist_early[[i]],NA)
  Deltalist_early[[i]] = ifelse(mask==1,Deltalist_early[[i]],NA)
  
  EDQMlist_mid[[i]] = ifelse(mask==1,EDQMlist_mid[[i]],NA)
  PARMlist_mid[[i]] = ifelse(mask==1,PARMlist_mid[[i]],NA)
  Deltalist_mid[[i]] = ifelse(mask==1,Deltalist_mid[[i]],NA)
  
  EDQMlist_late[[i]] = ifelse(mask==1,EDQMlist_late[[i]],NA)
  PARMlist_late[[i]] = ifelse(mask==1,PARMlist_late[[i]],NA)
  Deltalist_late[[i]] = ifelse(mask==1,Deltalist_late[[i]],NA)
  
message("Finished calcs for variable ",vars[i])
}

##########
# Coldest Tmax - PARM vs. PRISM, EDQM

diffhist_PP = PARMhistlist[[1]]-PRISMhistlist[[1]]
diffhist_PE = PARMhistlist[[1]]-EDQMhistlist[[1]]

diffcolorbar = colorramp(c(diffhist_PP,diffhist_PE),colorchoice="bluetored",Blimit=30)

par(mfrow=c(2,1))
testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PP)
surface(testsfc1,type="I",main="Coldest Tmin PARM-PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PE)
surface(testsfc1,type="I",main="Coldest Tmax PARM-EDQM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

##########
# Median Tmax - PARM vs. PRISM, EDQM

diffhist_PP = PARMhistlist[[3]]-PRISMhistlist[[3]]
diffhist_PE = PARMhistlist[[3]]-EDQMhistlist[[3]]

diffcolorbar = colorramp(c(diffhist_PP,diffhist_PE),colorchoice="bluetored",Blimit=30)
par(mfrow=c(2,1))
testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PP)
surface(testsfc1,type="I",main="Median Tmax PARM-PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PE)
surface(testsfc1,type="I",main="Median Tmax PARM-EDQM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

##########
# Mean Tmax - PARM vs. PRISM, EDQM

diffhist_PP = PARMhistlist[[4]]-PRISMhistlist[[4]]
diffhist_PE = PARMhistlist[[4]]-EDQMhistlist[[4]]

diffcolorbar = colorramp(c(diffhist_PP,diffhist_PE),colorchoice="bluetored",Blimit=30)
par(mfrow=c(2,1))
testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PP)
surface(testsfc1,type="I",main="Mean Tmax PARM-PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PE)
surface(testsfc1,type="I",main="Mean Tmax PARM-EDQM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

##########
# Max Tmax - PARM vs. PRISM, EDQM

diffhist_PP = PARMhistlist[[6]]-PRISMhistlist[[6]]
diffhist_PE = PARMhistlist[[6]]-EDQMhistlist[[6]]

diffcolorbar = colorramp(c(diffhist_PP,diffhist_PE),colorchoice="bluetored",Blimit=30)
par(mfrow=c(2,1))
testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PP)
surface(testsfc1,type="I",main="Max Tmax PARM-PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diffhist_PE)
surface(testsfc1,type="I",main="Max Tmax PARM-EDQM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)


##########
# Coldest Tmax - PARM vs. Delta, EDQM

diff_PDe = PARMlist_early[[1]]-Deltalist_early[[1]]
diff_PEe = PARMlist_early[[1]]-EDQMlist_early[[1]]

diff_PDm = PARMlist_mid[[1]]-Deltalist_mid[[1]]
diff_PEm = PARMlist_mid[[1]]-EDQMlist_mid[[1]]

diff_PDl = PARMlist_late[[1]]-Deltalist_late[[1]]
diff_PEl = PARMlist_late[[1]]-EDQMlist_late[[1]]

diffcolorbar = colorramp(c(diff_PDe,diff_PEe,diff_PDm,diff_PEm,diff_PDl,diff_PEl),colorchoice="bluetored",Blimit=30)

par(mfrow=c(3,2))
testsfc1 = list(x=(lon-360),y=lat,z=diff_PDe)
surface(testsfc1,type="I",main="Coldest Tmax PARM-Delta 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEe)
surface(testsfc1,type="I",main="Coldest Tmax PARM-EDQM 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDm)
surface(testsfc1,type="I",main="Coldest Tmax PARM-Delta 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEm)
surface(testsfc1,type="I",main="Coldest Tmax PARM-EDQM 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDl)
surface(testsfc1,type="I",main="Coldest Tmax PARM-Delta 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEl)
surface(testsfc1,type="I",main="Coldest Tmax PARM-EDQM 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)


##########
# Median Tmax - PARM vs. Delta, EDQM

diff_PDe = PARMlist_early[[3]]-Deltalist_early[[3]]
diff_PEe = PARMlist_early[[3]]-EDQMlist_early[[3]]

diff_PDm = PARMlist_mid[[3]]-Deltalist_mid[[3]]
diff_PEm = PARMlist_mid[[3]]-EDQMlist_mid[[3]]

diff_PDl = PARMlist_late[[3]]-Deltalist_late[[3]]
diff_PEl = PARMlist_late[[3]]-EDQMlist_late[[3]]

diffcolorbar = colorramp(c(diff_PDe,diff_PEe,diff_PDm,diff_PEm,diff_PDl,diff_PEl),colorchoice="bluetored",Blimit=30)

par(mfrow=c(3,2))
testsfc1 = list(x=(lon-360),y=lat,z=diff_PDe)
surface(testsfc1,type="I",main="Median Tmax PARM-Delta 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEe)
surface(testsfc1,type="I",main="Median Tmax PARM-EDQM 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDm)
surface(testsfc1,type="I",main="Median Tmax PARM-Delta 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEm)
surface(testsfc1,type="I",main="Median Tmax PARM-EDQM 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDl)
surface(testsfc1,type="I",main="Median Tmax PARM-Delta 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEl)
surface(testsfc1,type="I",main="Median Tmax PARM-EDQM 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

##########
# Mean Tmax - PARM vs. Delta, EDQM

diff_PDe = PARMlist_early[[4]]-Deltalist_early[[4]]
diff_PEe = PARMlist_early[[4]]-EDQMlist_early[[4]]

diff_PDm = PARMlist_mid[[4]]-Deltalist_mid[[4]]
diff_PEm = PARMlist_mid[[4]]-EDQMlist_mid[[4]]

diff_PDl = PARMlist_late[[4]]-Deltalist_late[[4]]
diff_PEl = PARMlist_late[[4]]-EDQMlist_late[[4]]

diffcolorbar = colorramp(c(diff_PDe,diff_PEe,diff_PDm,diff_PEm,diff_PDl,diff_PEl),colorchoice="bluetored",Blimit=30)

par(mfrow=c(3,2))
testsfc1 = list(x=(lon-360),y=lat,z=diff_PDe)
surface(testsfc1,type="I",main="Mean Tmax PARM-Delta 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEe)
surface(testsfc1,type="I",main="Mean Tmax PARM-EDQM 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDm)
surface(testsfc1,type="I",main="Mean Tmax PARM-Delta 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEm)
surface(testsfc1,type="I",main="Mean Tmax PARM-EDQM 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDl)
surface(testsfc1,type="I",main="Mean Tmax PARM-Delta 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEl)
surface(testsfc1,type="I",main="Mean Tmax PARM-EDQM 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

##########
# Max Tmax - PARM vs. Delta, EDQM

diff_PDe = PARMlist_early[[6]]-Deltalist_early[[6]]
diff_PEe = PARMlist_early[[6]]-EDQMlist_early[[6]]

diff_PDm = PARMlist_mid[[6]]-Deltalist_mid[[6]]
diff_PEm = PARMlist_mid[[6]]-EDQMlist_mid[[6]]

diff_PDl = PARMlist_late[[6]]-Deltalist_late[[6]]
diff_PEl = PARMlist_late[[6]]-EDQMlist_late[[6]]

diffcolorbar = colorramp(c(diff_PDe,diff_PEe,diff_PDm,diff_PEm,diff_PDl,diff_PEl),colorchoice="bluetored",Blimit=30)

par(mfrow=c(3,2))
testsfc1 = list(x=(lon-360),y=lat,z=diff_PDe)
surface(testsfc1,type="I",main="Max Tmax PARM-Delta 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEe)
surface(testsfc1,type="I",main="Max Tmax PARM-EDQM 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDm)
surface(testsfc1,type="I",main="Max Tmax PARM-Delta 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEm)
surface(testsfc1,type="I",main="Max Tmax PARM-EDQM 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PDl)
surface(testsfc1,type="I",main="Max Tmax PARM-Delta 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=diff_PEl)
surface(testsfc1,type="I",main="Max Tmax PARM-EDQM 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

####
#####################

# historical boxplot
i=5
if(i==1) comp="coldest"
if(i==2) comp="5th percentile"
if(i==3) comp="median"
if(i==4) comp="mean"
if(i==5) comp="95th percentile"
if(i==6) comp="max"

DS = rep("hist_PARM",length(as.vector(PARMhistlist[[i]])))
PARMhistdf = data.frame(DS,as.vector(PARMhistlist[[i]]))
names(PARMhistdf) = c("DS","tmax")

DS = rep("hist_EDQM",length(as.vector(EDQMhistlist[[i]])))
EDQMhistdf = data.frame(DS,as.vector(EDQMhistlist[[i]]))
names(EDQMhistdf) = c("DS","tmax")
boxdf = rbind(PARMhistdf,EDQMhistdf)

DS = rep("hist_PRISM",length(as.vector(PRISMhistlist[[i]])))
PRISMhistdf = data.frame(DS,as.vector(PRISMhistlist[[i]]))
names(PRISMhistdf) = c("DS","tmax")
boxdf = rbind(boxdf,PRISMhistdf)

DS = rep("early_PARM",length(as.vector(PARMlist_early[[i]])))
PARMdf = data.frame(DS,as.vector(PARMlist_early[[i]]))
names(PARMdf) = c("DS","tmax")
boxdf = rbind(boxdf,PARMdf)

DS = rep("early_EDQM",length(as.vector(EDQMlist_early[[i]])))
EDQMdf = data.frame(DS,as.vector(EDQMlist_early[[i]]))
names(EDQMdf) = c("DS","tmax")
boxdf = rbind(boxdf,EDQMdf)

DS = rep("early_Delta",length(as.vector(Deltalist_early[[i]])))
Deltadf = data.frame(DS,as.vector(Deltalist_early[[i]]))
names(Deltadf) = c("DS","tmax")
boxdf = rbind(boxdf,Deltadf)

DS = rep("mid_PARM",length(as.vector(PARMlist_mid[[i]])))
PARMdf = data.frame(DS,as.vector(PARMlist_mid[[i]]))
names(PARMdf) = c("DS","tmax")
boxdf = rbind(boxdf,PARMdf)

DS = rep("mid_EDQM",length(as.vector(EDQMlist_mid[[i]])))
EDQMdf = data.frame(DS,as.vector(EDQMlist_mid[[i]]))
names(EDQMdf) = c("DS","tmax")
boxdf = rbind(boxdf,EDQMdf)

DS = rep("mid_Delta",length(as.vector(Deltalist_mid[[i]])))
Deltadf = data.frame(DS,as.vector(Deltalist_mid[[i]]))
names(Deltadf) = c("DS","tmax")
boxdf = rbind(boxdf,Deltadf)

DS = rep("late_PARM",length(as.vector(PARMlist_late[[i]])))
PARMdf = data.frame(DS,as.vector(PARMlist_late[[i]]))
names(PARMdf) = c("DS","tmax")
boxdf = rbind(boxdf,PARMdf)

DS = rep("late_EDQM",length(as.vector(EDQMlist_late[[i]])))
EDQMdf = data.frame(DS,as.vector(EDQMlist_late[[i]]))
names(EDQMdf) = c("DS","tmax")
boxdf = rbind(boxdf,EDQMdf)

DS = rep("late_Delta",length(as.vector(Deltalist_late[[i]])))
Deltadf = data.frame(DS,as.vector(Deltalist_late[[i]]))
names(Deltadf) = c("DS","tmax")
boxdf = rbind(boxdf,Deltadf)

boxdf$DS<-factor(boxdf$DS, levels=c("hist_PRISM", "hist_EDQM", "hist_PARM","early_Delta","early_EDQM","early_PARM","mid_Delta","mid_EDQM","mid_PARM","late_Delta","late_EDQM","late_PARM"))

boxplot(tmax~DS,data=boxdf,xaxt="n",main=paste(comp," Tmax",sep=""))
lablist = c("hist_PRISM", "hist_EDQM", "hist_PARM","early_Delta","early_EDQM","early_PARM","mid_Delta","mid_EDQM","mid_PARM","late_Delta","late_EDQM","late_PARM")
text(seq(1, 12, by=1), par("usr")[3]-2, labels = lablist, srt = 90, pos = 1, xpd = TRUE)

####
# projected change plots
i=5
if(i==1) comp="coldest"
if(i==2) comp="5th percentile"
if(i==3) comp="median"
if(i==4) comp="mean"
if(i==5) comp="95th percentile"
if(i==6) comp="max"

DS = rep("early_PARM",length(as.vector(PARMlist_early[[i]]-PARMhistlist[[i]])))
PARMdf = data.frame(DS,as.vector(PARMlist_early[[i]]-PARMhistlist[[i]]))
names(PARMdf) = c("DS","tmax")

DS = rep("early_EDQM",length(as.vector(EDQMlist_early[[i]]-EDQMhistlist[[i]])))
EDQMdf = data.frame(DS,as.vector(EDQMlist_early[[i]]-EDQMhistlist[[i]]))
names(EDQMdf) = c("DS","tmax")
boxchangedf = rbind(PARMdf,EDQMdf)

DS = rep("early_Delta",length(as.vector(Deltalist_early[[i]]-PRISMhistlist[[i]])))
Deltadf = data.frame(DS,as.vector(Deltalist_early[[i]]-PRISMhistlist[[i]]))
names(Deltadf) = c("DS","tmax")
boxchangedf = rbind(boxchangedf,Deltadf)

DS = rep("mid_PARM",length(as.vector(PARMlist_mid[[i]]-PARMhistlist[[i]])))
PARMdf = data.frame(DS,as.vector(PARMlist_mid[[i]]-PARMhistlist[[i]]))
names(PARMdf) = c("DS","tmax")
boxchangedf = rbind(boxchangedf,PARMdf)

DS = rep("mid_EDQM",length(as.vector(EDQMlist_mid[[i]]-EDQMhistlist[[i]])))
EDQMdf = data.frame(DS,as.vector(EDQMlist_mid[[i]]-EDQMhistlist[[i]]))
names(EDQMdf) = c("DS","tmax")
boxchangedf = rbind(boxchangedf,EDQMdf)

DS = rep("mid_Delta",length(as.vector(Deltalist_mid[[i]]-PRISMhistlist[[i]])))
Deltadf = data.frame(DS,as.vector(Deltalist_mid[[i]]-PRISMhistlist[[i]]))
names(Deltadf) = c("DS","tmax")
boxchangedf = rbind(boxchangedf,Deltadf)

DS = rep("late_PARM",length(as.vector(PARMlist_late[[i]]-PARMhistlist[[i]])))
PARMdf = data.frame(DS,as.vector(PARMlist_late[[i]]-PARMhistlist[[i]]))
names(PARMdf) = c("DS","tmax")
boxchangedf = rbind(boxchangedf,PARMdf)

DS = rep("late_EDQM",length(as.vector(EDQMlist_late[[i]]-EDQMhistlist[[i]])))
EDQMdf = data.frame(DS,as.vector(EDQMlist_late[[i]]-EDQMhistlist[[i]]))
names(EDQMdf) = c("DS","tmax")
boxchangedf = rbind(boxchangedf,EDQMdf)

DS = rep("late_Delta",length(as.vector(Deltalist_late[[i]]-PRISMhistlist[[i]])))
Deltadf = data.frame(DS,as.vector(Deltalist_late[[i]]-PRISMhistlist[[i]]))
names(Deltadf) = c("DS","tmax")
boxchangedf = rbind(boxchangedf,Deltadf)

boxchangedf$DS<-factor(boxchangedf$DS, levels=c("early_Delta","early_EDQM","early_PARM","mid_Delta","mid_EDQM","mid_PARM","late_Delta","late_EDQM","late_PARM"))

boxplot(tmax~DS,data=boxchangedf,xaxt="n",main=paste("projected change ",comp," Tmax",sep=""))
lablist = c("early_Delta","early_EDQM","early_PARM","mid_Delta","mid_EDQM","mid_PARM","late_Delta","late_EDQM","late_PARM")
abline(h=0,lty=3)
text(seq(1, 9, by=1), par("usr")[3]-0.3, labels = lablist, srt = 90, pos = 1, xpd = TRUE)

##################
# MAE calculations

histerror = abs(PARMtmaxhist-PRISMtmaxhist[,,-which(substr(histdates,6,10)=="02-29")])
futerror = abs(PARMtmaxfut-Deltatmaxfut)

MAEhist = apply(histerror,c(1,2),mean,na.rm=TRUE)
MAEearly = apply(futerror[,,earlyidx],c(1,2),mean,na.rm=TRUE)
MAEmid = apply(futerror[,,mididx],c(1,2),mean,na.rm=TRUE)
MAElate = apply(futerror[,,lateidx],c(1,2),mean,na.rm=TRUE)

diffcolorbar = colorramp(c(MAEhist,MAEearly,MAEmid,MAElate),colorchoice="yellowtored",type="raw",Blimit=30)

testsfc1 = list(x=(lon-360),y=lat,z=MAEhist)
surface(testsfc1,type="I",main="MAE PARM vs. PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAEearly)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAEmid)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAElate)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)



MAEhist = abs(PARMhistlist[[6]]-PRISMhistlist[[6]])
MAEearly = abs(PARMlist_early[[6]]-Deltalist_early[[6]])
MAEmid = abs(PARMlist_mid[[6]]-Deltalist_mid[[6]])
MAElate = abs(PARMlist_late[[6]]-Deltalist_late[[6]])

diffcolorbar = colorramp(c(MAEhist,MAEearly,MAEmid,MAElate),colorchoice="yellowtored",type="raw",Blimit=20)

testsfc1 = list(x=(lon-360),y=lat,z=MAEhist)
surface(testsfc1,type="I",main="MAE PARM vs. PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAEearly)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAEmid)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAElate)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

#########

histerror = PARMtmaxhist-PRISMtmaxhist[,,-which(substr(histdates,6,10)=="02-29")]
futerror = PARMtmaxfut-Deltatmaxfut

MAEhist = apply(histerror,c(1,2),mean,na.rm=TRUE)
MAEearly = apply(futerror[,,earlyidx],c(1,2),mean,na.rm=TRUE)
MAEmid = apply(futerror[,,mididx],c(1,2),mean,na.rm=TRUE)
MAElate = apply(futerror[,,lateidx],c(1,2),mean,na.rm=TRUE)

diffcolorbar = colorramp(c(MAEhist,MAEearly,MAEmid,MAElate),colorchoice="bluetored",type="difference",Blimit=30)

testsfc1 = list(x=(lon-360),y=lat,z=MAEhist)
surface(testsfc1,type="I",main="MAE PARM vs. PRISM",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAEearly)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAEmid)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=MAElate)
surface(testsfc1,type="I",main="MAE PARM vs. Delta 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

######

PRISMtmaxhist = PRISMtmaxhist[,,-which(substr(histdates,6,10)=="02-29")]

years = as.numeric(unique(substr(histdates2,1,4)))
yearmon = unique(substr(histdates2,1,7))

PRISMyearly = PARMyearly = EDQMyearly = array(NA,dim=c(length(lon),length(lat),length(years)))
PRISMyearmon = PARMyearmon = EDQMyearmon = array(NA,dim=c(length(lon),length(lat),length(yearmon)))

for(m in 1:length(yearmon)){
  idx = which(substr(histdates2,1,7)==yearmon[m])
  PARMyearmon[,,m] = apply(PARMtmaxhist[,,idx],c(1,2),mean,na.rm=TRUE)
  PRISMyearmon[,,m] = apply(PRISMtmaxhist[,,idx],c(1,2),mean,na.rm=TRUE)
  EDQMyearmon[,,m] = apply(EDQMtmaxhist[,,idx],c(1,2),mean,na.rm=TRUE)
}

for(y in 1:length(years)){
  idx = which(as.numeric(substr(histdates2,1,4))==years[y])
  PARMyearly[,,y] = apply(PARMtmaxhist[,,idx],c(1,2),mean,na.rm=TRUE)
  PRISMyearly[,,y] = apply(PRISMtmaxhist[,,idx],c(1,2),mean,na.rm=TRUE)
  EDQMyearly[,,y] = apply(EDQMtmaxhist[,,idx],c(1,2),mean,na.rm=TRUE)
}

PARMlmfit = EDQMlmfit = matrix(NA,nrow=length(lon),ncol=length(lat))
for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    
    if(is.na(PARMyearly[r,c,1])==FALSE){
    diff1 = PARMyearly[r,c,]-PRISMyearly[r,c,]
    diff2 = EDQMyearly[r,c,]-PRISMyearly[r,c,]
    lmfit1 = lm(diff1~years)
    lmfit2 = lm(diff2~years)
    PARMlmfit[r,c] = lmfit1$coefficients[2]
    EDQMlmfit[r,c] = lmfit2$coefficients[2]
    }
  }
  message("Finished calcs for row ",r," / ",length(lon))
}

diffcolorbar = colorramp(c(PARMlmfit,EDQMlmfit),colorchoice="bluetored",type="difference",Blimit=30)

testsfc1 = list(x=(lon-360),y=lat,z=PARMlmfit)
surface(testsfc1,type="I",main="PARM-PRISM trend",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMlmfit)
surface(testsfc1,type="I",main="EDQM-PRISM trend",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

################

PARMvar = apply(PARMyearly,c(1,2),var,na.rm=TRUE)
PRISMvar = apply(PRISMyearly,c(1,2),var,na.rm=TRUE)
EDQMvar = apply(EDQMyearly,c(1,2),var,na.rm=TRUE)

vardiff1 = PARMvar-PRISMvar
vardiff2 = EDQMvar-PRISMvar

diffcolorbar = colorramp(c(vardiff1,vardiff2),colorchoice="bluetored",type="difference",Blimit=30)

testsfc1 = list(x=(lon-360),y=lat,z=vardiff1)
surface(testsfc1,type="I",main="PARM-PRISM variance",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=vardiff2)
surface(testsfc1,type="I",main="EDQM-PRISM variance",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

############

PARMvar = apply(PARMyearly,3,sd,na.rm=TRUE)
PRISMvar = apply(PRISMyearly,3,sd,na.rm=TRUE)
EDQMvar = apply(EDQMyearly,3,sd,na.rm=TRUE)

vardiff1 = PARMvar-PRISMvar
vardiff2 = EDQMvar-PRISMvar

plot(vardiff1~years,type="b",pch=19,col="blue")
plot(vardiff2~years,type="b",pch=19,col="blue")

#############

library(e1071)

PARMskew = EDQMskew = PRISMskew = matrix(NA,nrow=length(lon),ncol=length(lat))
PARMskewearly = EDQMskewearly = Deltaskewearly = matrix(NA,nrow=length(lon),ncol=length(lat))
PARMskewmid = EDQMskewmid = Deltaskewmid = matrix(NA,nrow=length(lon),ncol=length(lat))
PARMskewlate = EDQMskewlate = Deltaskewlate = matrix(NA,nrow=length(lon),ncol=length(lat))

for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    
    if(is.na(PARMtmaxhist[r,c,1])==FALSE){
      PARMskew[r,c] = skewness(PARMtmaxhist[r,c,])
      EDQMskew[r,c] = skewness(EDQMtmaxhist[r,c,])
      PRISMskew[r,c] = skewness(PRISMtmaxhist[r,c,])
      
      PARMskewearly[r,c] = skewness(PARMtmaxfut[r,c,earlyidx])
      EDQMskewearly[r,c] = skewness(EDQMtmaxfut[r,c,earlyidx])
      Deltaskewearly[r,c] = skewness(Deltatmaxfut[r,c,earlyidx])
      
      PARMskewmid[r,c] = skewness(PARMtmaxfut[r,c,mididx])
      EDQMskewmid[r,c] = skewness(EDQMtmaxfut[r,c,mididx])
      Deltaskewmid[r,c] = skewness(Deltatmaxfut[r,c,mididx])
      
      PARMskewlate[r,c] = skewness(PARMtmaxfut[r,c,lateidx])
      EDQMskewlate[r,c] = skewness(EDQMtmaxfut[r,c,lateidx])
      Deltaskewlate[r,c] = skewness(Deltatmaxfut[r,c,lateidx])
    }
  }
  message("Finished calcs for row ",r," / ",length(lon))
}


diffcolorbar = colorramp(c(PARMskew,EDQMskew,PRISMskew,PARMskewearly,EDQMskewearly,Deltaskewearly,PARMskewmid,EDQMskewmid,Deltaskewmid,PARMskewlate,EDQMskewlate,Deltaskewlate),colorchoice="bluetored",type="difference",Blimit=30)

testsfc1 = list(x=(lon-360),y=lat,z=PARMskew)
surface(testsfc1,type="I",main="PARM skew",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMskew)
surface(testsfc1,type="I",main="EDQM skew",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PRISMskew)
surface(testsfc1,type="I",main="PRISM skew",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMskewearly)
surface(testsfc1,type="I",main="PARM skew 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMskewearly)
surface(testsfc1,type="I",main="EDQM skew 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=Deltaskewearly)
surface(testsfc1,type="I",main="Delta skew 2006-2039",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMskewmid)
surface(testsfc1,type="I",main="PARM skew 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMskewmid)
surface(testsfc1,type="I",main="EDQM skew 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=Deltaskewmid)
surface(testsfc1,type="I",main="Delta skew 2040-2069",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=PARMskewlate)
surface(testsfc1,type="I",main="PARM skew 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=EDQMskewlate)
surface(testsfc1,type="I",main="EDQM skew 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)

testsfc1 = list(x=(lon-360),y=lat,z=Deltaskewlate)
surface(testsfc1,type="I",main="Delta skew 2070-2099",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)