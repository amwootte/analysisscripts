source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rgdal)
#library(RPushbullet)

var = "prcp"
tempres = "monthly"
varfile = "prcp"

filesin = system(paste("ls /data4/data/OBS/DAYMETv4/na/pr/*.nc",sep=""),intern=TRUE)

filesplit = do.call(rbind,strsplit(filesin,split="/",fixed=TRUE))
filesplit2 = do.call(rbind,strsplit(filesplit[,ncol(filesplit)],split="_",fixed=TRUE))
filesplit2 = data.frame(filesplit2)
filesplit2$year = as.numeric(substr(filesplit2[,ncol(filesplit2)],1,4))

years = 1980

loclon = -98.847                    
loclat = 29.822

#loclon = -98.48395
#loclat = 29.54429

fileidx = which(filesplit2$year>=years[1] & filesplit2$year<=years[2])

nctest = nc_open(filesin[1])
lon = ncvar_get(nctest,"lon")
lat = ncvar_get(nctest,"lat")
nc_close(nctest)

LON = as.vector(lon)
LAT = as.vector(lat)
R = rep(1:dim(lon)[1],dim(lon)[2])
C = rep(1:dim(lon)[2],each=dim(lon)[1])

modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")

###
# get cells to use
loc_lon = as.numeric(loclon)
loc_lat = as.numeric(loclat)

pointarea = distfunc(loc_lon,loc_lat,modelgrid)
locstart = c(pointarea$R,pointarea$C)

####

#tmplon=lon[4704:4706,6367:6369]
#tmplat=lat[4704:4706,6367:6369]

#tmplon = tmplon+(tmplon-tmplon[2,2])/2
#tmplat = tmplat+(tmplat-tmplat[2,2])/2
  #tmplon[1,1] = tmplon[1,1]+(tmplon[1,1]-tmplon[2,2])/2
  


#####
# get data

var = "prcp"
tempres = "monthly"
varfile = "prcp"

unitsout="mm"

filesin = system(paste("ls /data4/data/OBS/DAYMETv4/na/pr/*.nc",sep=""),intern=TRUE)

outputdata = NULL

for(i in 1:length(filesin)){
  nctest = nc_open(filesin[i])
  varnames = names(nctest$var)
  times = ncvar_get(nctest,"time")
  timeunits = nctest$var[[which(varnames==varfile)]]$dim[[3]]$units
  #timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
  prdat = ncvar_get(nctest,varfile,start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  indexdate = as.Date(substr(timeunits,12,21))
  nc_close(nctest)
dates = indexdate+times  
tmpdat = data.frame(dates,prdat)
outputdata=rbind(outputdata,tmpdat)
  message("Finished grab for file ",i," / ",length(filesin))
}

###

var = "tasmax"
tempres = "monthly"
varfile = "tmax"

filesin = system(paste("ls /data4/data/OBS/DAYMETv4/na/tasmax/*.nc",sep=""),intern=TRUE)

outputdata2 = NULL

for(i in 1:length(filesin)){
  nctest = nc_open(filesin[i])
  varnames = names(nctest$var)
  times = ncvar_get(nctest,"time")
  timeunits = nctest$var[[which(varnames==varfile)]]$dim[[3]]$units
  #timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
  tmaxdat = ncvar_get(nctest,varfile,start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  indexdate = as.Date(substr(timeunits,12,21))
  nc_close(nctest)
  dates = indexdate+times  
  tmpdat = data.frame(dates,tmaxdat)
  outputdata2=rbind(outputdata2,tmpdat)
  message("Finished grab for file ",i," / ",length(filesin))
}

###

var = "tasmin"
tempres = "monthly"
varfile = "tmin"

filesin = system(paste("ls /data4/data/OBS/DAYMETv4/na/tasmin/*.nc",sep=""),intern=TRUE)

outputdata3 = NULL

for(i in 1:length(filesin)){
  nctest = nc_open(filesin[i])
  varnames = names(nctest$var)
  times = ncvar_get(nctest,"time")
  timeunits = nctest$var[[which(varnames==varfile)]]$dim[[3]]$units
  #timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
  tmindat = ncvar_get(nctest,varfile,start=c(pointarea$R[1],pointarea$C[1],1),count=c(1,1,-1))
  indexdate = as.Date(substr(timeunits,12,21))
  nc_close(nctest)
  dates = indexdate+times  
  tmpdat = data.frame(dates,tmindat)
  outputdata3=rbind(outputdata3,tmpdat)
  message("Finished grab for file ",i," / ",length(filesin))
}

####
# merge files

outputdata = merge(outputdata,outputdata2,by="dates")
outputdata = merge(outputdata,outputdata3,by="dates")

names(outputdata) = c("date","prcp","tmax","tmin")

write.table(outputdata,"EAA_north_of_rechargezone.csv",sep=",",row.names=FALSE)

#####
# read airport data
outputdata = read.table("SanAntonioAirport_DAYMET.csv",sep=",",header=TRUE,colClasses = c("character",rep("numeric",3)))

obsdat = read.table("SanAntonioAirport.csv",sep=",",header=TRUE,colClasses = c(rep("character",6)))
obsdat$TMAX = (as.numeric(obsdat$TMAX)-32)*5/9
obsdat$TMIN = (as.numeric(obsdat$TMIN)-32)*5/9
obsdat$PRCP = as.numeric(obsdat$PRCP)*25.4
names(obsdat)[3]="date"

outputdata$date = as.character(outputdata$date)

fulldata = merge(obsdat,outputdata,by="date")

fulldata$year = as.numeric(substr(fulldata$date,1,4))
fulldata$yearmon = substr(fulldata$date,1,7)

fullanntemp = aggregate(fulldata[,c(5,6,8,9)],by=list(year=fulldata$year),mean,na.rm=TRUE)
fullannprcp = aggregate(fulldata[,c(4,7)],by=list(year=fulldata$year),sum,na.rm=TRUE)

plot(density(fulldata$TMAX),col="black",lwd=2,main="PDFs for station and DAYMETv4 tmax values")
lines(density(fulldata$tmax),col="blue",lwd=2)
legend("topleft",c("station","DAYMETv4"),col=c("black","blue"),lty=1,lwd=2)

plot(density(fulldata$TMIN),col="black",lwd=2,main="PDFs for station and DAYMETv4 tmin values",ylim=c(0,0.08))
lines(density(fulldata$tmin),col="blue",lwd=2)
legend("topleft",c("station","DAYMETv4"),col=c("black","blue"),lty=1,lwd=2)

plot(density(fulldata$PRCP),col="black",lwd=2,main="PDFs for station and DAYMETv4 prcp values",ylim=c(0,0.35))
lines(density(fulldata$prcp),col="blue",lwd=2)
legend("topleft",c("station","DAYMETv4"),col=c("black","blue"),lty=1,lwd=2)

tmaxRMSE = sqrt(mean((fulldata$tmax-fulldata$TMAX)^2))
tminRMSE = sqrt(mean((fulldata$tmin-fulldata$TMIN)^2))
prcpRMSE = sqrt(mean((fulldata$prcp-fulldata$PRCP)^2))

tmaxerror = mean((fulldata$tmax-fulldata$TMAX))
tminerror = mean((fulldata$tmin-fulldata$TMIN))
prcperror = mean((fulldata$prcp-fulldata$PRCP))

tmaxcor = cor(fulldata$tmax,fulldata$TMAX)
tmincor = cor(fulldata$tmin,fulldata$TMIN)
prcpcor = cor(fulldata$prcp,fulldata$PRCP)

plot(fulldata$TMAX~as.Date(fulldata$date),type="l",lwd=2,main="Station vs. DAYMETv4 - High Temperature",ylab="degrees C",xlab="Date")
lines(fulldata$tmax~as.Date(fulldata$date),lwd=2,col="blue")
abline(h=0,lty=2)
legend("topleft",c("station","DAYMETv4"),col=c("black","blue"),lty=1,lwd=2)

plot(fulldata$TMIN~as.Date(fulldata$date),type="l",lwd=2,main="Station vs. DAYMETv4 - Low Temperature",ylab="degrees C",xlab="Date")
lines(fulldata$tmin~as.Date(fulldata$date),lwd=2,col="blue")
abline(h=0,lty=2)
legend("topleft",c("station","DAYMETv4"),col=c("black","blue"),lty=1,lwd=2)

plot(fulldata$PRCP~as.Date(fulldata$date),type="l",lwd=2,main="Station vs. DAYMETv4 - Precipitation",ylab="mm",xlab="Date")
lines(fulldata$prcp~as.Date(fulldata$date),lwd=2,col="blue")
abline(h=0,lty=2)
legend("topleft",c("station","DAYMETv4"),col=c("black","blue"),lty=1,lwd=2)

####
# rearrange to boxplot formatting

tmpdat1 = fulldata[,c(1,5)]
tmpdat2 = fulldata[,c(1,8)]

tmpdat1$group = "station"
tmpdat2$group = "DAYMETv4"
names(tmpdat1)[2]="tmax"

tmaxall = rbind(tmpdat1,tmpdat2)
tmaxall$season = NA

idx = which(as.numeric(substr(tmaxall$date,6,7))>=12 | as.numeric(substr(tmaxall$date,6,7))<=2)
tmaxall$season[idx] = "DJF"

idx = which(as.numeric(substr(tmaxall$date,6,7))>=3 & as.numeric(substr(tmaxall$date,6,7))<=5)
tmaxall$season[idx] = "MAM"

idx = which(as.numeric(substr(tmaxall$date,6,7))>=6 & as.numeric(substr(tmaxall$date,6,7))<=8)
tmaxall$season[idx] = "JJA"

idx = which(as.numeric(substr(tmaxall$date,6,7))>=9 & as.numeric(substr(tmaxall$date,6,7))<=11)
tmaxall$season[idx] = "SON"

tmaxall$group = factor(tmaxall$group,levels=c("station","DAYMETv4"),ordered=TRUE)
tmaxall$season = factor(tmaxall$season,levels=c("DJF","MAM","JJA","SON"),ordered=TRUE)

library(ggplot2)
ggplot(tmaxall, aes(x=group, y=tmax, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison tmax")+xlab("")+ylab("Deg. C")

ggplot(tmaxall, aes(x=group, y=tmax, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison tmax seasonal")+xlab("")+ylab("Deg. C")+facet_wrap(facets=vars(season),nrow=2,ncol=2)

####


tmpdat1 = fulldata[,c(1,6)]
tmpdat2 = fulldata[,c(1,9)]

tmpdat1$group = "station"
tmpdat2$group = "DAYMETv4"
names(tmpdat1)[2]="tmin"

tminall = rbind(tmpdat1,tmpdat2)

tminall$season = NA

idx = which(as.numeric(substr(tminall$date,6,7))>=12 | as.numeric(substr(tminall$date,6,7))<=2)
tminall$season[idx] = "DJF"

idx = which(as.numeric(substr(tminall$date,6,7))>=3 & as.numeric(substr(tminall$date,6,7))<=5)
tminall$season[idx] = "MAM"

idx = which(as.numeric(substr(tminall$date,6,7))>=6 & as.numeric(substr(tminall$date,6,7))<=8)
tminall$season[idx] = "JJA"

idx = which(as.numeric(substr(tminall$date,6,7))>=9 & as.numeric(substr(tminall$date,6,7))<=11)
tminall$season[idx] = "SON"

tminall$group = factor(tminall$group,levels=c("station","DAYMETv4"),ordered=TRUE)
tminall$season = factor(tminall$season,levels=c("DJF","MAM","JJA","SON"),ordered=TRUE)


ggplot(tminall, aes(x=group, y=tmin, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison tmin")+xlab("")+ylab("Deg. C")

ggplot(tminall, aes(x=group, y=tmin, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison tmin seasonal")+xlab("")+ylab("Deg. C")+facet_wrap(facets=vars(season),nrow=2,ncol=2)

###

tmpdat1 = fulldata[,c(1,4)]
tmpdat2 = fulldata[,c(1,7)]

tmpdat1$group = "station"
tmpdat2$group = "DAYMETv4"
names(tmpdat1)[2]="prcp"

prcpall = rbind(tmpdat1,tmpdat2)

prcpall$season = NA

idx = which(as.numeric(substr(prcpall$date,6,7))>=12 | as.numeric(substr(prcpall$date,6,7))<=2)
prcpall$season[idx] = "DJF"

idx = which(as.numeric(substr(prcpall$date,6,7))>=3 & as.numeric(substr(prcpall$date,6,7))<=5)
prcpall$season[idx] = "MAM"

idx = which(as.numeric(substr(prcpall$date,6,7))>=6 & as.numeric(substr(prcpall$date,6,7))<=8)
prcpall$season[idx] = "JJA"

idx = which(as.numeric(substr(prcpall$date,6,7))>=9 & as.numeric(substr(prcpall$date,6,7))<=11)
prcpall$season[idx] = "SON"

prcpall$group = factor(prcpall$group,levels=c("station","DAYMETv4"),ordered=TRUE)
prcpall$season = factor(prcpall$season,levels=c("DJF","MAM","JJA","SON"),ordered=TRUE)

ggplot(prcpall, aes(x=group, y=prcp, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison prcp")+xlab("")+ylab("mm")

ggplot(prcpall, aes(x=group, y=prcp, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison prcp seasonal")+xlab("")+ylab("mm")+facet_wrap(facets=vars(season),nrow=2,ncol=2)

######
# aggregate to yearly

tmaxall$year = as.numeric(substr(tmaxall$date,1,4))
tmaxall$year_group = paste(tmaxall$group,tmaxall$year,sep="-")
tmaxann = aggregate(tmaxall[,c(2,5)],by=list(year_group=tmaxall$year_group),mean,na.rm=TRUE)
tmaxann$group = do.call("rbind",strsplit(tmaxann$year_group,"-",fixed=TRUE))[,1]
tmaxann$group = factor(tmaxann$group,levels=c("station","DAYMETv4"),ordered=TRUE)

ggplot(tmaxann, aes(x=group, y=tmax, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison tmax yearly")+xlab("")+ylab("C")

plot(tmaxann$tmax[which(tmaxann$group=="station")]~tmaxann$year[which(tmaxann$group=="station")],type="l",lwd=2,ylim=range(tmaxann$tmax),main="Comparison - Annual Average of Daily High Temperature",ylab="degrees C",xlab="Year")
lines(tmaxann$tmax[which(tmaxann$group=="DAYMETv4")]~tmaxann$year[which(tmaxann$group=="DAYMETv4")],lwd=2,col="blue")
abline(h=0,lty=2)
legend("topleft",c(paste("station mean=",round(mean(tmaxann$tmax[which(tmaxann$group=="station")]),2)," C",sep=""),
                   paste("DAYMETv4 mean=",round(mean(tmaxann$tmax[which(tmaxann$group=="DAYMETv4")]),2)," C",sep="")),col=c("black","blue"),lty=1,lwd=2)
cor(tmaxann$tmax[which(tmaxann$group=="DAYMETv4")],tmaxann$tmax[which(tmaxann$group=="station")])

###

tminall$year = as.numeric(substr(tminall$date,1,4))
tminall$year_group = paste(tminall$group,tminall$year,sep="-")
tminann = aggregate(tminall[,c(2,5)],by=list(year_group=tminall$year_group),mean,na.rm=TRUE)
tminann$group = do.call("rbind",strsplit(tminann$year_group,"-",fixed=TRUE))[,1]
tminann$group = factor(tminann$group,levels=c("station","DAYMETv4"),ordered=TRUE)

ggplot(tminann, aes(x=group, y=tmin, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison tmin yearly")+xlab("")+ylab("C")

plot(tminann$tmin[which(tminann$group=="station")]~tminann$year[which(tminann$group=="station")],type="l",lwd=2,ylim=range(tminann$tmin),main="Comparison - Annual Average of Daily Low Temperature",ylab="degrees C",xlab="Year")
lines(tminann$tmin[which(tminann$group=="DAYMETv4")]~tminann$year[which(tminann$group=="DAYMETv4")],lwd=2,col="blue")
abline(h=0,lty=2)
legend("topleft",c(paste("station mean=",round(mean(tminann$tmin[which(tminann$group=="station")]),2)," C",sep=""),
                   paste("DAYMETv4 mean=",round(mean(tminann$tmin[which(tminann$group=="DAYMETv4")]),2)," C",sep="")),col=c("black","blue"),lty=1,lwd=2)
cor(tminann$tmin[which(tminann$group=="DAYMETv4")],tminann$tmin[which(tminann$group=="station")])

####



prcpall$year = as.numeric(substr(prcpall$date,1,4))
prcpall$year_group = paste(prcpall$group,prcpall$year,sep="-")
prcpann = aggregate(prcpall[,c(2,5)],by=list(year_group=prcpall$year_group),sum,na.rm=TRUE)
prcpann$group = do.call("rbind",strsplit(prcpann$year_group,"-",fixed=TRUE))[,1]
prcpann$group = factor(prcpann$group,levels=c("station","DAYMETv4"),ordered=TRUE)

ggplot(prcpann, aes(x=group, y=prcp, fill=group)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown","cyan")) +
  ggtitle("San Antonio Comparison prcp yearly")+xlab("")+ylab("mm")


plot(prcpann$prcp[which(prcpann$group=="station")]~tmaxann$year[which(tmaxann$group=="station")],type="l",lwd=2,ylim=range(prcpann$prcp),main="Comparison - Annual Total of Daily Total Precipitation",ylab="mm",xlab="Year")
lines(prcpann$prcp[which(prcpann$group=="DAYMETv4")]~tmaxann$year[which(tmaxann$group=="DAYMETv4")],lwd=2,col="blue")
abline(h=0,lty=2)
legend("topleft",c(paste("station mean=",round(mean(prcpann$prcp[which(prcpann$group=="station")]),2)," mm",sep=""),
                   paste("DAYMETv4 mean=",round(mean(prcpann$prcp[which(prcpann$group=="DAYMETv4")]),2)," mm",sep="")),col=c("black","blue"),lty=1,lwd=2)
cor(prcpann$prcp[which(prcpann$group=="DAYMETv4")],prcpann$prcp[which(prcpann$group=="station")])

