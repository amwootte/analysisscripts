#################
#
# Southern Appalachians shapefile testing

library(sp)
library(fields)
library(raster)
library(rasterVis)
library(maptools)
library(maps)
library(ncdf)


varname="pr25"
anfilepath = "analysisid"

test = readShapePoly("physio_shp/physio")
projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

#plot(test)
#map("state",add=TRUE,col="red")
#lines(test.sub,col="blue")

#################
# Get Component Data

test2 = open.ncdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_components_fix2.nc",sep=""))
lon = get.var.ncdf(test2,"lon")
lat = get.var.ncdf(test2,"lat")
time = 1950:2099
V = get.var.ncdf(test2,"V")
M = get.var.ncdf(test2,"M_smooth")
D = get.var.ncdf(test2,"D_smooth")
S = get.var.ncdf(test2,"S_smooth")
G = get.var.ncdf(test2,"G_smooth")
T = get.var.ncdf(test2,"T_smooth")
close.ncdf(test2)

################
# Subset to future period only

timeidx = which(time>=2006)

M=M[,,timeidx]
D=D[,,timeidx]
S=S[,,timeidx]
G=G[,,timeidx]
T=T[,,timeidx]

#################
# cropping calculation

earlyidx = which(time[timeidx]==2020)
mididx = which(time[timeidx]==2055)
lateidx = which(time[timeidx]==2090)

CM = CD = CS = CV = M

for(i in 1:dim(M)[3]){
  CM[,,i] = (M[,,i]/T[,,i])*100
  CD[,,i] = (D[,,i]/T[,,i])*100
  CS[,,i] = (S[,,i]/T[,,i])*100
  CV[,,i] = (V/T[,,i])*100
}


###
# corresponds to the following regions in the physio shapefile
# 5a - Northern Blue Ridge
# 5b - Southern Blue Ridge
# 6b - Middle Appalachian Highland Valley and Ridge
# 8d - Alleghany Mountain, Appalchian Plateaus.

test.sub <- test[as.character(test@data$FENCODE) %in% c("5a","5b","6b","8d"), ]

tmp = CD[,,mididx]
modras = raster(t(tmp)[length(lat):1,])
extent(modras) = c(min(lon),max(lon),min(lat),max(lat))

#plot(NA,xlim=range(lon),ylim=range(lat),xlab="Longitude",ylab="Latitude",main="Southeast U.S. domain")
#lines(test,add=TRUE)
#lines(test.sub,col="red",lwd=2)
#map("state",add=TRUE)

#x1 = c(-83.5,-83.5,-80,-80,-83.5)
#y1 = c(27,29,29,27,27)
#x2 = c(-95.5,-95.5,-93,-93,-95.5)
#y2 = c(35.5,39,39,35.5,35.5)

#lines(y1~x1,col="red",lwd=2)
#lines(y2~x2,col="red",lwd=2)

#map("state",add=TRUE,col="gray")


#plot(mod.sub)
#lines(test,add=TRUE)
#lines(test.sub,add=TRUE,col="red",lwd=2)

Vts = c()
Mts = c()
Dts = c()
Sts = c()
Tts = c()
Gts = c()

for(i in 1:dim(M)[3]){

  test.sub <- test[as.character(test@data$FENCODE) %in% c("5a","5b","6b","8d"), ]
  #test.sub <- test[as.character(test@data$FENCODE) %in% c("5b"), ]
tmpV = V
tmpM = M[,,i]
tmpD = D[,,i]
tmpS = S[,,i]
tmpT = T[,,i]
tmpG = G[,,i]

modrasV = raster(t(tmpV)[length(lat):1,])
extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))

modrasM = raster(t(tmpM)[length(lat):1,])
extent(modrasM) = c(min(lon),max(lon),min(lat),max(lat))

modrasD = raster(t(tmpD)[length(lat):1,])
extent(modrasD) = c(min(lon),max(lon),min(lat),max(lat))

modrasS = raster(t(tmpS)[length(lat):1,])
extent(modrasS) = c(min(lon),max(lon),min(lat),max(lat))

modrasG = raster(t(tmpG)[length(lat):1,])
extent(modrasG) = c(min(lon),max(lon),min(lat),max(lat))

modrasT = raster(t(tmpT)[length(lat):1,])
extent(modrasT) = c(min(lon),max(lon),min(lat),max(lat))

mod.subV <- crop(modrasV, extent(test.sub))
mod.subV <- mask(modrasV, test.sub)

mod.subM <- crop(modrasM, extent(test.sub))
mod.subM <- mask(modrasM, test.sub)

mod.subD <- crop(modrasD, extent(test.sub))
mod.subD <- mask(modrasD, test.sub)

mod.subS <- crop(modrasS, extent(test.sub))
mod.subS <- mask(modrasS, test.sub)

mod.subT <- crop(modrasT, extent(test.sub))
mod.subT <- mask(modrasT, test.sub)

mod.subG <- crop(modrasG, extent(test.sub))
mod.subG <- mask(modrasG, test.sub)

#plot(mod.sub)
#lines(test,add=TRUE)
#lines(test.sub,add=TRUE,col="red",lwd=2)
#vals = getValues(mod.sub)
#modtemp = matrix(vals,nrow=length(lon),ncol=length(lat))
#modtemp = modtemp[,length(lat):1]

Vts[i] = mean(getValues(mod.subV),na.rm=TRUE)
Mts[i] = mean(getValues(mod.subM),na.rm=TRUE)
Dts[i] = mean(getValues(mod.subD),na.rm=TRUE)
Sts[i] = mean(getValues(mod.subS),na.rm=TRUE)
Tts[i] = mean(getValues(mod.subT),na.rm=TRUE)
Gts[i] = mean(getValues(mod.subG),na.rm=TRUE)

}

#################
# Fractional Uncertainty Calc

FM = abs((1.65*sqrt(Mts))/Gts)
FD = abs((1.65*sqrt(Dts))/Gts)
FS = abs((1.65*sqrt(Sts))/Gts)
FV = abs((1.65*sqrt(Vts))/Gts)
FT = abs((1.65*sqrt(Tts))/Gts)

###################
# Plotting results

if(varname=="tmax"){
  titlevar = "Average High Temperature"
  titlevarunits = "(C)"
}

if(varname=="tmax95"){
  titlevar = "Average Number of Days over 95F"
  titlevarunits = "days"
}

if(varname=="tmin"){
  titlevar = "Average Low Temperature"
  titlevarunits = "(C)"
}

if(varname=="tmin32"){
  titlevar = "Average Number of Days under 32F"
  titlevarunits = "days"
}


if(varname=="pr"){
  titlevar = "Average Total Precipitation"
  titlevarunits = "(%)"
}

if(varname=="pr25"){
  titlevar = "Average Number of Days with Rainfall > 1 in"
  titlevarunits = "days"
}


pdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_timeplotschunk_fix2_APP.pdf",sep=""),height=6,width=6,onefile=TRUE)

###########
# Fractional Uncertainty

year = 2006:2099

#Frac = data.frame(year,Fvd,Fvm,Fm,Fd,Fs,Ft)
Frac = data.frame(year,FV,FM,FD,FS,FT)

Frac2 = subset(Frac,year>2010)
Frac3 = Frac2

# normalized
plot(FT/max(FT,na.rm=TRUE)~year,data=Frac3,type="l",lwd=2,ylim=c(0,1),main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
lines(FV/max(FT,na.rm=TRUE)~year,data=Frac3,col="orange",lwd=2)
lines(FM/max(FT,na.rm=TRUE)~year,data=Frac3,col="darkblue",lwd=2)
lines(FD/max(FT,na.rm=TRUE)~year,data=Frac3,col="red",lwd=2)
lines(FS/max(FT,na.rm=TRUE)~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","brown","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))

# non normalized
plot(FT~year,data=Frac3,type="l",ylim=c(0,max(Frac3$FT,na.rm=TRUE)),lwd=2,main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
#lines(Fvm~year,data=Frac3,col="orange",lwd=2)
#lines(Fvd~year,data=Frac3,col="brown",lwd=2)
lines(FV~year,data=Frac3,col="orange",lwd=2)
lines(FM~year,data=Frac3,col="darkblue",lwd=2)
lines(FD~year,data=Frac3,col="red",lwd=2)
lines(FS~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","chocolate4","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))


###########################
# Fraction of total plot

combodat2 = Frac = data.frame(year,Vts,Mts,Dts,Sts,Tts)
datain2 = subset(combodat2,year>=2000)
datain2$leadtime = datain2$year-2000

names(datain2) = c("year","V","M","D","S","T","leadtime")

datain2$top = 100
datain2$mV = datain2$top - (datain2$V/datain2$T)*100
datain2$mVS = datain2$mV - (datain2$S/datain2$T)*100
datain2$mM = datain2$mVS - (datain2$M/datain2$T)*100

leadtime = datain2$leadtime
leadidx = which(leadtime>=5 & leadtime<95)
leadtime = leadtime[leadidx]

y.top = datain2$top[leadidx]
y.mV = datain2$mV[leadidx]
y.mVS = datain2$mVS[leadidx]
y.mM = datain2$mM[leadidx]

plot(y.top~leadtime,col="orange",type="l",ylim=c(0,100),xlim=c(0,100),xlab="Leadtime in years from 2000",ylab="Fraction of Total Uncertainty (%)",main=paste("Appalachians Fraction of Total Uncertainty\nfor ",titlevar,sep=""))
lines(y.mV~leadtime,col="darkgreen")
lines(y.mVS~leadtime,col="darkblue")
lines(y.mM~leadtime,col="red")

polygon(c(rev(leadtime), leadtime), c(rev(y.top), y.mV),col = "orange", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mV),y.mVS),col = "darkgreen", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mVS),y.mM),col = "darkblue", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mM),rep(0,length(y.top))),col = "red", border = NA)

dev.off()
