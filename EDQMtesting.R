##################
# Perfect model testing

library(ncdf4)
library(maps)
library(maptools)
library(mapdata)
library(fields)
library(sp)

setwd("/Users/amwootte/Documents")

source("PMPP/colorrampfunctions.R")
source("PMPP/spelllengthfreq.R")

test=nc_open("PMPP/BCQM/amip/pr_day_PMgCRprp1-BCQM-A00cve2X01K02_amip_DSCVh_US48_19790101-20081231_ensmbl.nc")

lat = ncvar_get(test,"lat")
lon = ncvar_get(test,"lon")

pr = ncvar_get(test,"pr")

times = ncvar_get(test,"time")

nc_close(test)

dates = seq(as.Date("1979-01-01"),as.Date("2008-12-31"),by="day")

prens1 = array(pr[,,1,],dim=c(length(lon),length(lat),length(dates)))
prens2 = array(pr[,,2,],dim=c(length(lon),length(lat),length(dates)))

prens1 = prens1*86400
prens2 = prens2*86400

par(mfrow=c(1,1))

testsfc = list(x=lon-360,y=lat,z=prens1[,,1])
surface(testsfc,type="I")
map("state",add=TRUE)

years = as.numeric(substr(dates,1,4))
yrs = unique(years)

prens1yearly = array(NA,dim=c(length(lon),length(lat),length(yrs)))
prens2yearly = array(NA,dim=c(length(lon),length(lat),length(yrs)))

for(y in 1:length(yrs)){
  yearidx = which(years==yrs[y])
  prens1yearly[,,y]=apply(prens1[,,yearidx],c(1,2),sum) #sum for annual total
  prens2yearly[,,y]=apply(prens2[,,yearidx],c(1,2),sum) #mean for daily mean, spell_length_calc,type="gte",spell_length=7
  message("Complete calcs for year: ",yrs[y])
}

plot(prens1yearly[100,50,]~yrs,type="b")
lines(prens2yearly[100,50,]~yrs,type="b",col="red")

prclimo1 = apply(prens1yearly,c(1,2),mean)
prclimo2 = apply(prens2yearly,c(1,2),mean)

testsfc = list(x=lon-360,y=lat,z=prclimo1)
surface(testsfc,type="I")
map("world",add=TRUE)

testsfc = list(x=lon-360,y=lat,z=prclimo2)
surface(testsfc,type="I")
map("world",add=TRUE)

testsfc = list(x=lon-360,y=lat,z=prclimo2-prclimo1)
surface(testsfc,type="I")
map("world",add=TRUE)

#####
# Future periods

test=nc_open("PMPP/BCQM/sst2090/pr_day_PMgCRprp1-BCQM-A01r1e3X01K00_sst2090_r1to3i1p1_US48_20860101-20951231_ensmbl.nc")
prfut1 = ncvar_get(test,"pr") #ESM
nc_close(test)

test=nc_open("PMPP/BCQM/sst2090/pr_day_PMgCRprp1-BCQM-A02r1e3X01K00_sst2090_r1to3i1p2_US48_20860101-20951231_ensmbl.nc")
prfut2 = ncvar_get(test,"pr") #CM3
nc_close(test)

datesfut = seq(as.Date("2086-01-01"),as.Date("2095-12-31"),by="day")

prfut1 = prfut1*86400
prfut2 = prfut2*86400

prfut1ens1 = array(prfut1[,,1,],dim=c(length(lon),length(lat),length(datesfut)))
prfut1ens2 = array(prfut1[,,2,],dim=c(length(lon),length(lat),length(datesfut)))
prfut1ens3 = array(prfut1[,,3,],dim=c(length(lon),length(lat),length(datesfut)))

prfut2ens1 = array(prfut2[,,1,],dim=c(length(lon),length(lat),length(datesfut)))
prfut2ens2 = array(prfut2[,,2,],dim=c(length(lon),length(lat),length(datesfut)))
prfut2ens3 = array(prfut2[,,3,],dim=c(length(lon),length(lat),length(datesfut)))

yearsfut=as.numeric(substr(datesfut,1,4))
yrsf = unique(yearsfut)

prfut1ens1yearly = array(NA,dim=c(length(lon),length(lat),length(yrsf)))
prfut1ens2yearly = array(NA,dim=c(length(lon),length(lat),length(yrsf)))
prfut1ens3yearly = array(NA,dim=c(length(lon),length(lat),length(yrsf)))
prfut2ens1yearly = array(NA,dim=c(length(lon),length(lat),length(yrsf)))
prfut2ens2yearly = array(NA,dim=c(length(lon),length(lat),length(yrsf)))
prfut2ens3yearly = array(NA,dim=c(length(lon),length(lat),length(yrsf)))

#for(r in 1:length(lon)){
 # for(c in 1:length(lat)){
  #  spell_length_calc(prfut1ens3[r,c,yearidx],type="dry",spell_length=7)
  #}
#}


for(y in 1:length(yrsf)){
  yearidx = which(yearsfut==yrsf[y])
  prfut1ens1yearly[,,y]=apply(prfut1ens1[,,yearidx],c(1,2),sum) # sum for annual total
  prfut1ens2yearly[,,y]=apply(prfut1ens2[,,yearidx],c(1,2),sum) # mean for daily mean
  prfut1ens3yearly[,,y]=apply(prfut1ens3[,,yearidx],c(1,2),sum) #spell_length_calc,type="wet",spell_length=7
  prfut2ens1yearly[,,y]=apply(prfut2ens1[,,yearidx],c(1,2),sum)
  prfut2ens2yearly[,,y]=apply(prfut2ens2[,,yearidx],c(1,2),sum)
  prfut2ens3yearly[,,y]=apply(prfut2ens3[,,yearidx],c(1,2),sum)
  message("Calculations complete for year ",yrsf[y])
}

prfut1climo1 = apply(prfut1ens1yearly,c(1,2),mean)
prfut1climo2 = apply(prfut1ens2yearly,c(1,2),mean)
prfut1climo3 = apply(prfut1ens3yearly,c(1,2),mean)
prfut2climo1 = apply(prfut2ens1yearly,c(1,2),mean)
prfut2climo2 = apply(prfut2ens2yearly,c(1,2),mean)
prfut2climo3 = apply(prfut2ens3yearly,c(1,2),mean)

prhistclimos = array(NA,dim=c(length(lon),length(lat),2))
prfut1climos = array(NA,dim=c(length(lon),length(lat),3))
prfut2climos = array(NA,dim=c(length(lon),length(lat),3))

prhistclimos[,,1]=prclimo1
prhistclimos[,,2]=prclimo2

prfut1climos[,,1]=prfut1climo1
prfut1climos[,,2]=prfut1climo2
prfut1climos[,,3]=prfut1climo3

prfut2climos[,,1]=prfut2climo1
prfut2climos[,,2]=prfut2climo2
prfut2climos[,,3]=prfut2climo3

prhclimo = apply(prhistclimos,c(1,2),mean)
prf1climo = apply(prfut1climos,c(1,2),mean)
prf2climo = apply(prfut2climos,c(1,2),mean)

data(wrld_simpl)
points <- expand.grid(lon-360,lat)
pts <- SpatialPoints(points, proj4string=CRS(proj4string(wrld_simpl)))
ii <- !is.na(over(pts, wrld_simpl)$FIPS)
plot(wrld_simpl)
points(pts, col=1+ii, pch=16)
points$mask = ii
points$R = rep(1:length(lon),length(lat))
points$C = rep(1:length(lat),each=length(lon))
landseamask = matrix(NA,nrow=length(lon),ncol=length(lat))

for(i in 1:nrow(points)){
  landseamask[points$R[i],points$C[i]]=points$mask[i]
}

fut1diff = (prf1climo-prhclimo) # E-hist
fut2diff = (prf2climo-prhclimo) # C-hist
fut12diff = (prf2climo-prf1climo) # C-E

testsfc = list(x=lon-360,y=lat,z=fut12diff)
surface(testsfc,type="I")
map("world",add=TRUE)
map("state",add=TRUE)

par(mfrow=c(2,2))

obscolorbar = obs_colorramp(prhclimo,colorchoice="rainbow",Blimit=40)

testsfc = list(x=lon-360,y=lat,z=prhclimo)
surface(testsfc,type="I",zlim=obscolorbar[[1]],breaks=obscolorbar[[2]],col=obscolorbar[[3]])
map("world",add=TRUE)
map("state",add=TRUE)

#fut1 = ESM
#fut2 = CM3

C_sub_H_land = mean(ifelse(landseamask==TRUE,fut2diff,NA),na.rm=TRUE)
E_sub_H_land = mean(ifelse(landseamask==TRUE,fut1diff,NA),na.rm=TRUE)
C_sub_E_land = mean(ifelse(landseamask==TRUE,fut12diff,NA),na.rm=TRUE)
  
C_sub_H_ocean = mean(ifelse(landseamask==FALSE,fut2diff,NA),na.rm=TRUE)
E_sub_H_ocean = mean(ifelse(landseamask==FALSE,fut1diff,NA),na.rm=TRUE)
C_sub_E_ocean = mean(ifelse(landseamask==FALSE,fut12diff,NA),na.rm=TRUE)
diffvals = c(C_sub_H_land,E_sub_H_land,C_sub_E_land,C_sub_H_ocean,E_sub_H_ocean,C_sub_E_ocean)
land = rep(c("land","ocean"),each=3)
group = rep(c("C-Hist","E-Hist","C-E"),2)
Dvals = data.frame(land,group,diffvals)
Dvals = Dvals[c(1,4,2,5,3,6),]
# round(range(Dvals[,3]),2)
barplot(Dvals[,3],ylim=c(-100,100),names.arg=Dvals$group,col=c("darkgreen","blue"),legend.text=c("Land","Ocean"),args.legend=list(x="topleft",col=c("darkgreen","blue")))
abline(h=0)

diffcolorbar = diff_colorramp(c(fut1diff,fut2diff),colorchoice="browntogreen",Blimit=40)

testsfc = list(x=lon-360,y=lat,z=fut1diff)
surface(testsfc,type="I",zlim=diffcolorbar[[1]],breaks=diffcolorbar[[2]],col=diffcolorbar[[3]])
map("world",add=TRUE)
map("state",add=TRUE)

testsfc = list(x=lon-360,y=lat,z=fut2diff)
surface(testsfc,type="I",zlim=diffcolorbar[[1]],breaks=diffcolorbar[[2]],col=diffcolorbar[[3]])
map("world",add=TRUE)
map("state",add=TRUE)



#################################################################################