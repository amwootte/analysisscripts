################
#
# Uncertainty analysis domain time series

############
# library and function load

############
# load libraries

library(ncdf)
library(fields)
library(maps)

############
# load calculator functions

source("componentfunctions.R")

###
# set variable name

varname = "tmax95"
anfilepath = "analysis"
#latloc = 34.70734 #MCBCL
#lonloc = -77.44516
latloc = 35.17083
lonloc = -79.0145
locname = "Bragg"

coords = c(latloc,lonloc)

distfunc = function(coords = c(lat,lon),dataset){

lat = coords[1]
lon = coords[2]

dist = 3963 * acos(sin( lat / 57.2958 ) * sin( dataset$lat / 57.2958 ) + cos( lat / 57.2958 ) * cos( dataset$lat /57.2958 ) * cos( ( dataset$lon / 57.2958) - ( lon / 57.2958 ) ))

dataset[which(dist==min(dist)),]

}

###
# set extra options

printPDF = TRUE # if you want some plots printed to PDF set this to TRUE


##############
# Read in components from netcdf file

compfile = paste("uncertainty/SE/",anfilepath,"/",varname,"_components_fix2.nc",sep="")

# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file

test = open.ncdf(compfile)

lon = get.var.ncdf(test,"lon")
lat = get.var.ncdf(test,"lat")
time = 1950:2099

LON = rep(lon,each=length(lat))
R = rep(1:length(lon),each=length(lat))
LAT = rep(lat,length(lon))
C = rep(1:length(lat),length(lon))

gridmeta = data.frame(LON,LAT,R,C)
names(gridmeta)[1:2] = c("lon","lat")

locgrab = distfunc(coords,gridmeta)

if(locname=="MCBCL"){
xloc = locgrab$R[1]-1 # MCBCL
yloc = locgrab$C[1]
} else {
xloc = locgrab$R[1] # Bragg
yloc = locgrab$C[1]
}

V = get.var.ncdf(test,"V")
M = get.var.ncdf(test,"M_smooth")
D = get.var.ncdf(test,"D_smooth")
S = get.var.ncdf(test,"S_smooth")
G = get.var.ncdf(test,"G_smooth")
T = get.var.ncdf(test,"T_smooth")
close.ncdf(test)

################
# Subset to future period only, grab only for desired points

timeidx = which(time>=2006)

Vts = V[xloc,yloc]
Mts=M[xloc,yloc,timeidx]
Dts=D[xloc,yloc,timeidx]
Sts=S[xloc,yloc,timeidx]
Gts=G[xloc,yloc,timeidx]
Tts=T[xloc,yloc,timeidx]

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

pdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_timeplots_fix2_",locname,".pdf",sep=""),height=6,width=6,onefile=TRUE)

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

plot(y.top~leadtime,col="orange",type="l",ylim=c(0,100),xlim=c(0,100),xlab="Leadtime in years from 2000",ylab="Fraction of Total Uncertainty (%)",main=paste("Fraction of Total Uncertainty\nfor ",titlevar,sep=""))
lines(y.mV~leadtime,col="darkgreen")
lines(y.mVS~leadtime,col="darkblue")
lines(y.mM~leadtime,col="red")

polygon(c(rev(leadtime), leadtime), c(rev(y.top), y.mV),col = "orange", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mV),y.mVS),col = "darkgreen", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mVS),y.mM),col = "darkblue", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mM),rep(0,length(y.top))),col = "red", border = NA)

dev.off()




