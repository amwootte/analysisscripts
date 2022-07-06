################
#
# Uncertainty analysis domain time series

############
# library and function load

############
# load libraries

library(ncdf4)
library(fields)
library(maps)

############
# load calculator functions

source("/home/woot0002/scripts/uncertainty/componentfunctions_2021.R")

###
# set variable name

varname = "TAVG_12M"
datasetname = "08144500"
anfilepath = "analysis"

###
# set extra options

printPDF = TRUE # if you want some plots printed to PDF set this to TRUE


##############
# Read in components from netcdf file

compfile = paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",datasetname,"_analysis_fix2.Rdata",sep="")

# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file

load(compfile)

load(paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",datasetname,"_components_fix2.Rdata",sep=""))


pdf(paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",datasetname,"_timeplotschunk_fix2.pdf",sep=""),height=6,width=6,onefile=TRUE)

###########
# Fractional Uncertainty

year = years[which(years>=2006)]

Frac = data.frame(year,FV,FM,FD,FS,FT)
Frac2 = subset(Frac,year>2010)
Frac3 = Frac2

# normalized
plot(FT/max(FT,na.rm=TRUE)~year,data=Frac3,type="l",lwd=2,ylim=c(0,4),main=paste("Fractional Uncertainty ",datasetname,sep=""),ylab="Fractional Uncertainty",xlab="year")
lines(FV/max(FT,na.rm=TRUE)~year,data=Frac3,col="orange",lwd=2)
lines(FM/max(FT,na.rm=TRUE)~year,data=Frac3,col="darkblue",lwd=2)
lines(FD/max(FT,na.rm=TRUE)~year,data=Frac3,col="red",lwd=2)
lines(FS/max(FT,na.rm=TRUE)~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","brown","blue","red","green"))
legend("topleft",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))

# non normalized
plot(FT~year,data=Frac3,type="l",ylim=range(Frac3[,2:ncol(Frac3)]),lwd=2,main=paste("Fractional Uncertainty ",datasetname,sep=""),ylab="Fractional Uncertainty",xlab="year")
#lines(Fvm~year,data=Frac3,col="orange",lwd=2)
#lines(Fvd~year,data=Frac3,col="brown",lwd=2)
lines(FV~year,data=Frac3,col="orange",lwd=2)
lines(FM~year,data=Frac3,col="darkblue",lwd=2)
lines(FD~year,data=Frac3,col="red",lwd=2)
lines(FS~year,data=Frac3,col="darkgreen",lwd=2)

#legend("topright",legend=c("Total","Int. Var GCM","Int. Var DS","GCM","DS","Scenario"),lwd=2,col=c("black","orange","chocolate4","blue","red","green"))
legend("bottomright",legend=c("Total","Nat. Variability","GCM","Downscaling","Scenario"),lwd=2,col=c("black","orange","darkblue","red","darkgreen"))


###########################
# Fraction of total plot

combodat2 = Frac = data.frame(year,M_smooth[which(years>=2006)],D_smooth[which(years>=2006)],S_smooth[which(years>=2006)],T_smooth[which(years>=2006)])
combodat2$V = V
datain2 = subset(combodat2,year>=2000)
datain2$leadtime = datain2$year-2000

names(datain2) = c("year","M","D","S","T","V","leadtime")

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

plot(y.top~leadtime,col="orange",type="l",ylim=c(0,100),xlim=c(0,100),xlab="Leadtime in years from 2000",ylab="Fraction of Total Uncertainty (%)",main=paste("Fraction of Total Uncertainty\n Location: ",datasetname,sep=""))
lines(y.mV~leadtime,col="darkgreen")
lines(y.mVS~leadtime,col="darkblue")
lines(y.mM~leadtime,col="red")

polygon(c(rev(leadtime), leadtime), c(rev(y.top), y.mV),col = "orange", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mV),y.mVS),col = "darkgreen", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mVS),y.mM),col = "darkblue", border = NA)
polygon(c(rev(leadtime), leadtime), c(rev(y.mM),rep(0,length(y.top))),col = "red", border = NA)

dev.off()

plot(y.top~leadtime,col="orange",type="l",ylim=c(0,100),xlim=c(0,100),xlab="Leadtime in years from 2000",ylab="Fraction of Total Uncertainty (%)",main=paste("Fraction of Total Uncertainty\n Location: ",datasetname,sep=""))
lines(y.mV~leadtime,col="darkgreen")
lines(y.mVS~leadtime,col="darkblue")
lines(y.mM~leadtime,col="red")

polygon(c(leadtime, rev(leadtime)), c(y.top, rev(y.mV)),col = "orange", border = NA)
polygon(c(leadtime, rev(leadtime)), c(y.mV,rev(y.mVS)),col = "darkgreen", border = NA)
polygon(c(leadtime, rev(leadtime)), c(y.mVS,rev(y.mM)),col = "darkblue", border = NA)
polygon(c(leadtime, rev(leadtime)), c(y.mM,rep(0,length(y.top))),col = "red", border = NA)


