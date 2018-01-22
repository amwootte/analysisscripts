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

varname = "pr25"
anfilepath = "analysis"

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

V = get.var.ncdf(test,"V")
M = get.var.ncdf(test,"M_smooth")
D = get.var.ncdf(test,"D_smooth")
S = get.var.ncdf(test,"S_smooth")
G = get.var.ncdf(test,"G_smooth")
T = get.var.ncdf(test,"T_smooth")
close.ncdf(test)

################
# Subset to future period only

timeidx = which(time>=2006)

M=M[,,timeidx]
D=D[,,timeidx]
S=S[,,timeidx]
G=G[,,timeidx]
T=T[,,timeidx]

################
# Get multiple domain tseries

# Central Florida= -84 -> -80, 27 -> 29
# Southern Appalachians= -85.5 -> -81, 34.5 -> 37
# Southwest Missouri= -95.5 -> -93, 35.5 -> 39

lonidx1 = which(lon>= -84 & lon<= -80)
latidx1 = which(lat>= 27 & lat>= 29)

lonidx2 = which(lon>= -85.5 & lon<= -81)
latidx2 = which(lat>= 34.5 & lat>= 34)

lonidx3 = which(lon>= -95.5 & lon<= -93)
latidx3 = which(lat>= 35.5 & lat>= 39)

Vts1 = mean(V[lonidx1,latidx1],na.rm=TRUE) # central Florida time series
Mts1 = apply(M[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Dts1 = apply(D[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Sts1 = apply(S[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Gts1 = apply(G[lonidx1,latidx1,],3,mean,na.rm=TRUE)
Tts1 = apply(T[lonidx1,latidx1,],3,mean,na.rm=TRUE)

Vts2 = mean(V[lonidx2,latidx2],na.rm=TRUE) # Southern Appalachians time series
Mts2 = apply(M[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Dts2 = apply(D[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Sts2 = apply(S[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Gts2 = apply(G[lonidx2,latidx2,],3,mean,na.rm=TRUE)
Tts2 = apply(T[lonidx2,latidx2,],3,mean,na.rm=TRUE)

Vts3 = mean(V[lonidx3,latidx3],na.rm=TRUE) # Southwest Missouri time series
Mts3 = apply(M[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Dts3 = apply(D[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Sts3 = apply(S[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Gts3 = apply(G[lonidx3,latidx3,],3,mean,na.rm=TRUE)
Tts3 = apply(T[lonidx3,latidx3,],3,mean,na.rm=TRUE)

#################
# Fractional Uncertainty Calc

FM1 = abs((1.65*sqrt(Mts1))/Gts1)
FD1 = abs((1.65*sqrt(Dts1))/Gts1)
FS1 = abs((1.65*sqrt(Sts1))/Gts1)
FV1 = abs((1.65*sqrt(Vts1))/Gts1)
FT1 = abs((1.65*sqrt(Tts1))/Gts1)

FM2 = abs((1.65*sqrt(Mts2))/Gts2)
FD2 = abs((1.65*sqrt(Dts2))/Gts2)
FS2 = abs((1.65*sqrt(Sts2))/Gts2)
FV2 = abs((1.65*sqrt(Vts2))/Gts2)
FT2 = abs((1.65*sqrt(Tts2))/Gts2)

FM3 = abs((1.65*sqrt(Mts3))/Gts3)
FD3 = abs((1.65*sqrt(Dts3))/Gts3)
FS3 = abs((1.65*sqrt(Sts3))/Gts3)
FV3 = abs((1.65*sqrt(Vts3))/Gts3)
FT3 = abs((1.65*sqrt(Tts3))/Gts3)

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

pdf(paste("uncertainty/SE/",anfilepath,"/",varname,"_timeplotschunk_fix2.pdf",sep=""),height=6,width=6,onefile=TRUE)

############
# Central Florida

###########
# Fractional Uncertainty

year = 2006:2099

#Frac = data.frame(year,Fvd,Fvm,Fm,Fd,Fs,Ft)
Frac = data.frame(year,FV1,FM1,FD1,FS1,FT1)

Frac2 = subset(Frac,year>2010)
Frac3 = Frac2

# normalized
plot(FT1/max(FT1,na.rm=TRUE)~year,data=Frac3,type="l",lwd=2,ylim=c(0,1),main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
lines(FV1/max(FT1,na.rm=TRUE)~year,data=Frac3,col="orange",lwd=2)
lines(FM1/max(FT1,na.rm=TRUE)~year,data=Frac3,col="darkblue",lwd=2)
lines(FD1/max(FT1,na.rm=TRUE)~year,data=Frac3,col="red",lwd=2)
lines(FS1/max(FT1,na.rm=TRUE)~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","brown","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))

# non normalized
plot(FT1~year,data=Frac3,type="l",ylim=c(0,max(Frac3$FT,na.rm=TRUE)),lwd=2,main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
#lines(Fvm~year,data=Frac3,col="orange",lwd=2)
#lines(Fvd~year,data=Frac3,col="brown",lwd=2)
lines(FV1~year,data=Frac3,col="orange",lwd=2)
lines(FM1~year,data=Frac3,col="darkblue",lwd=2)
lines(FD1~year,data=Frac3,col="red",lwd=2)
lines(FS1~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","chocolate4","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))


###########################
# Fraction of total plot

combodat2 = Frac = data.frame(year,Vts1,Mts1,Dts1,Sts1,Tts1)
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

plot(y.top~leadtime,col="orange",type="l",ylim=c(0,100),xlim=c(0,100),xlab="Leadtime in years from 2000",ylab="Fraction of Total Uncertainty (%)",main=paste("Central Florida Fraction of Total Uncertainty\nfor ",titlevar,sep=""))
lines(y.mV~leadtime,col="darkgreen")
lines(y.mVS~leadtime,col="darkblue")
lines(y.mM~leadtime,col="red")

polygon(c(rev(leadtime), leadtime), c(rev(y.top), y.mV),col = "orange", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mV),y.mVS),col = "darkgreen", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mVS),y.mM),col = "darkblue", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mM),rep(0,length(y.top))),col = "red", border = NA)

#############
# Southern Appalachains
###########
# Fractional Uncertainty

year = 2006:2099

#Frac = data.frame(year,Fvd,Fvm,Fm,Fd,Fs,Ft)
Frac = data.frame(year,FV2,FM2,FD2,FS2,FT2)

Frac2 = subset(Frac,year>2010)
Frac3 = Frac2

# normalized
plot(FT2/max(FT2,na.rm=TRUE)~year,data=Frac3,type="l",lwd=2,ylim=c(0,1),main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
lines(FV2/max(FT2,na.rm=TRUE)~year,data=Frac3,col="orange",lwd=2)
lines(FM2/max(FT2,na.rm=TRUE)~year,data=Frac3,col="darkblue",lwd=2)
lines(FD2/max(FT2,na.rm=TRUE)~year,data=Frac3,col="red",lwd=2)
lines(FS2/max(FT2,na.rm=TRUE)~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","brown","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))

# non normalized
plot(FT2~year,data=Frac3,type="l",ylim=c(0,max(Frac3$FT,na.rm=TRUE)),lwd=2,main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
#lines(Fvm~year,data=Frac3,col="orange",lwd=2)
#lines(Fvd~year,data=Frac3,col="brown",lwd=2)
lines(FV2~year,data=Frac3,col="orange",lwd=2)
lines(FM2~year,data=Frac3,col="darkblue",lwd=2)
lines(FD2~year,data=Frac3,col="red",lwd=2)
lines(FS2~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","chocolate4","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))


###########################
# Fraction of total plot

combodat2 = Frac = data.frame(year,Vts2,Mts2,Dts2,Sts2,Tts2)
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

plot(y.top~leadtime,col="orange",type="l",ylim=c(0,100),xlim=c(0,100),xlab="Leadtime in years from 2000",ylab="Fraction of Total Uncertainty (%)",main=paste("Southern Appalachains Fraction of Total Uncertainty\nfor ",titlevar,sep=""))
lines(y.mV~leadtime,col="darkgreen")
lines(y.mVS~leadtime,col="darkblue")
lines(y.mM~leadtime,col="red")

polygon(c(rev(leadtime), leadtime), c(rev(y.top), y.mV),col = "orange", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mV),y.mVS),col = "darkgreen", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mVS),y.mM),col = "darkblue", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mM),rep(0,length(y.top))),col = "red", border = NA)

#####################
# Southwest Missouri
###########
# Fractional Uncertainty

year = 2006:2099

#Frac = data.frame(year,Fvd,Fvm,Fm,Fd,Fs,Ft)
Frac = data.frame(year,FV3,FM3,FD3,FS3,FT3)

Frac2 = subset(Frac,year>2010)
Frac3 = Frac2

# normalized
plot(FT3/max(FT3,na.rm=TRUE)~year,data=Frac3,type="l",lwd=2,ylim=c(0,1),main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
lines(FV3/max(FT3,na.rm=TRUE)~year,data=Frac3,col="orange",lwd=2)
lines(FM3/max(FT3,na.rm=TRUE)~year,data=Frac3,col="darkblue",lwd=2)
lines(FD3/max(FT3,na.rm=TRUE)~year,data=Frac3,col="red",lwd=2)
lines(FS3/max(FT3,na.rm=TRUE)~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","brown","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))

# non normalized
plot(FT3~year,data=Frac3,type="l",ylim=c(0,max(Frac3$FT,na.rm=TRUE)),lwd=2,main=paste("Fractional Uncertainty ",titlevar,sep=""),ylab="Fractional Uncertainty",xlab="year")
#lines(Fvm~year,data=Frac3,col="orange",lwd=2)
#lines(Fvd~year,data=Frac3,col="brown",lwd=2)
lines(FV3~year,data=Frac3,col="orange",lwd=2)
lines(FM3~year,data=Frac3,col="darkblue",lwd=2)
lines(FD3~year,data=Frac3,col="red",lwd=2)
lines(FS3~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","chocolate4","blue","red","green"))
legend("topright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))


###########################
# Fraction of total plot

combodat2 = Frac = data.frame(year,Vts3,Mts3,Dts3,Sts3,Tts3)
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

plot(y.top~leadtime,col="orange",type="l",ylim=c(0,100),xlim=c(0,100),xlab="Leadtime in years from 2000",ylab="Fraction of Total Uncertainty (%)",main=paste("Southwest Missouri Fraction of Total Uncertainty\nfor ",titlevar,sep=""))
lines(y.mV~leadtime,col="darkgreen")
lines(y.mVS~leadtime,col="darkblue")
lines(y.mM~leadtime,col="red")

polygon(c(rev(leadtime), leadtime), c(rev(y.top), y.mV),col = "orange", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mV),y.mVS),col = "darkgreen", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mVS),y.mM),col = "darkblue", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mM),rep(0,length(y.top))),col = "red", border = NA)

dev.off()




