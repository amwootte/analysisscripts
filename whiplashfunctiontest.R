
library(iterators)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)

whiplashcalc = function(PrecipData,lowperc=0.25,highperc=0.75,spread=1){
  
  # PrecipData: array of precipitation data arranged as [longitude,latitude,months]
  # lowperc: low percentile, defaults to 0.25 to 25th percentile
  # highperc: high percentile, defaults to 0.75 to 75th percentile
  
  #PrecipData = precipmon
  #lowperc = 0.25
  #highperc = 0.75
  
  percentile_bottom = apply(PrecipData,c(1,2),quantile,probs=lowperc,na.rm=TRUE)
  percentile_top = apply(PrecipData,c(1,2),quantile,probs=highperc,na.rm=TRUE)
  
  #testsfc = list(x=lon,y=lat,z=percentile_top)
  #surface(testsfc,type="I")
  
  precip_hardcode = PrecipData
  
  for(r in 1:dim(PrecipData)[1]){
    for(c in 1:dim(PrecipData)[2]){
      if(all(is.na(PrecipData[r,c,])==FALSE)==TRUE){
        precip_hardcode[r,c,]=ifelse(PrecipData[r,c,]<=percentile_bottom[r,c] | PrecipData[r,c,]>=percentile_top[r,c],ifelse(PrecipData[r,c,]<=percentile_bottom[r,c],0,2),1)
      }
    }
  }
  
  PDcount = DPcount = matrix(NA,nrow=dim(PrecipData)[1],ncol=dim(PrecipData)[2])
  precip_PD = precip_DP = PrecipData
  
  for(r in 1:dim(PrecipData)[1]){
    for(c in 1:dim(PrecipData)[2]){
      if(all(is.na(PrecipData[r,c,])==FALSE)==TRUE){
        tmp = precip_hardcode[r,c,1:(dim(PrecipData)[3]-spread)]-precip_hardcode[r,c,(spread+1):dim(PrecipData)[3]]
        tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,spread))
        tmppd = c(ifelse(tmp==2,1,0),rep(NA,spread))
        precip_PD[r,c,] = tmppd
        precip_DP[r,c,] = tmpdp
        PDcount[r,c] = sum(tmppd,na.rm=TRUE)
        DPcount[r,c] = sum(tmpdp,na.rm=TRUE)
      }
    }
  }
  
  #testsfc = list(x=lon,y=lat,z=PDcount)
  #surface(testsfc,type="I")
  
  #testsfc = list(x=lon,y=lat,z=DPcount)
  #surface(testsfc,type="I")
  results = list(percentile_bottom=percentile_bottom,percentile_top=percentile_top,precip_hardcode=precip_hardcode,precip_PD=precip_PD,precip_DP=precip_DP,PDcount=PDcount,DPcount=DPcount)
  results
}

source("/data2/3to5/I35/scripts/analysisfunctions.R")

##########
# Set input information

message("Setting inputs")

shapefile = "/home/woot0002/shapefiles/WBD_HU8_SanSabaLlano"

varname = "pr" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2070,2099)

colorchoicechange = "browntogreen"
BINSchange = 20
USEFS = FALSE
colorchoicehist = "whitetored"
BINShist = 20

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

############

varin = varname
if(varname=="tmax85" | varname=="tmax90" | varname=="tmax95" | varname=="tmax100" | varname=="gsl") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="prcptot" | varname=="prmean" | varname=="prmed" | varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd" | varname=="R99" | varname=="R95" | varname=="R90"  | varname=="R99freq" | varname=="R95freq" | varname=="R90freq") varin="pr"

########
# Find all file names

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

projnotes = paste(projnotes,collapse=",")
histnotes = paste(histnotes,collapse=",")

histlist = paste(histfilelist,collapse=",")
projlist = paste(projfilelist,collapse=",")

##########
# Load in test data

test = nc_open(histfilelist[2])
precip = ncvar_get(test,"pr")
lat = ncvar_get(test,"lat")
lon = ncvar_get(test,"lon")-360
nc_close(test)

precip = precip*86400

##########
# precip autocorrelation map

precipac = matrix(NA,nrow=length(lon),ncol=length(lat))

for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    if(is.na(precip[r,c,1])==FALSE){
      test = acf(precip[r,c,],plot=FALSE)
      precipac[r,c] = which(test[[1]]<0.01)[1]
    }
  }
}

testsfc = list(x=lon,y=lat,z=precipac)
surface(testsfc,type="I")
mean(precipac,na.rm=TRUE)


##########
# 30 day rolling sum test - single location

S = 30

calcrollsum = function(inputdata,size=5){
  require(zoo)
  rollapply(inputdata,size,sum,na.rm=TRUE)
}

R = 95
C = 70

ptm = proc.time()
preciprollsum = calcrollsum(precip[R,C,],size=S)
ptmend=  proc.time()-ptm

acf(precip[R,C,],plot=FALSE)

pacf(precip[R,C,])

acf(preciprollsum)
pacf(preciprollsum)


plot(preciprollsum[1:100])
plot(diff(preciprollsum[1:100]))

lowperc = 0.25
highperc = 0.75

percentile_bottom = quantile(preciprollsum,probs=lowperc,na.rm=TRUE)
percentile_top = quantile(preciprollsum,probs=highperc,na.rm=TRUE)

plot(preciprollsum[1:100])
abline(h=c(percentile_bottom,percentile_top),lty=2)
abline(v=seq(0,90,by=30),lty=3)
plot(diff(preciprollsum[1:100]))


precip_hardcode=ifelse(preciprollsum<=percentile_bottom | preciprollsum>=percentile_top,ifelse(preciprollsum<=percentile_bottom,0,2),1)

spread=c(1,4,7,10,14,30)

PDcounts = DPcounts = c()
for(s in 1:length(spread)){
  tmp = precip_hardcode[1:(length(preciprollsum)-spread[s])]-precip_hardcode[(spread[s]+1):length(preciprollsum)]
  
  tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,spread[s]))
  tmppd = c(ifelse(tmp==2,1,0),rep(NA,spread[s]))
  
  PDcounts[s] = sum(tmppd,na.rm=TRUE)
  DPcounts[s] = sum(tmpdp,na.rm=TRUE)
}



idx = which(tmpdp==1)
diffidx = diff(idx)
tmpdp[idx[1]:idx[which(diffidx>1)[1]]]

pickthese = idx[1]:idx[which(diffidx>1)[1]]

preciprollsum[pickthese]

PDcount = DPcount = matrix(NA,nrow=dim(PrecipData)[1],ncol=dim(PrecipData)[2])
precip_PD = precip_DP = PrecipData

for(r in 1:dim(PrecipData)[1]){
  for(c in 1:dim(PrecipData)[2]){
    if(all(is.na(PrecipData[r,c,])==FALSE)==TRUE){
      tmp = precip_hardcode[r,c,1:(dim(PrecipData)[3]-spread)]-precip_hardcode[r,c,(spread+1):dim(PrecipData)[3]]
      tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,spread))
      tmppd = c(ifelse(tmp==2,1,0),rep(NA,spread))
      precip_PD[r,c,] = tmppd
      precip_DP[r,c,] = tmpdp
      PDcount[r,c] = sum(tmppd,na.rm=TRUE)
      DPcount[r,c] = sum(tmpdp,na.rm=TRUE)
    }
  }
}






###########
# Test with historical PRISM
# monthly

histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
histmon = unique(substr(histdates,1,7))

precipmon = array(NA,dim=c(length(lon),length(lat),length(histmon)))
for(m in 1:length(histmon)){
  idx = which(substr(histdates,1,7)==histmon[m])
  precipmon[,,m] = apply(precip[,,idx],c(1,2),sum,na.rm=TRUE)
  precipmon[,,m] = ifelse(is.na(precip[,,1])==FALSE,precipmon[,,m],NA)
}

#testsfc = list(x=lon,y=lat,z=precipmon[,,1])
#surface(testsfc,type="I")

output = whiplashcalc(precipmon)

testsfc = list(x=lon,y=lat,z=output$DPcount)
surface(testsfc,type="I")

##########
# 30 day rolling sum test

S = 30

calcrollsum = function(inputdata,size=5){
  require(zoo)
  rollapply(inputdata,size,sum,na.rm=TRUE)
}
ptm = proc.time()
preciprollsum = apply(precip,c(1,2),calcrollsum,size=S)
ptmend=  proc.time()-ptm

preciprollsum2 = array(NA,dim=c(length(lon),length(lat),dim(preciprollsum)[1]))

for(m in 1:dim(preciprollsum)[1]){
  preciprollsum2[,,m] = ifelse(is.na(precip[,,1])==FALSE,preciprollsum[m,,],NA)
}

output2 = whiplashcalc(preciprollsum2)
output3 = whiplashcalc(preciprollsum2,spread=S)

##########
# 7 day rolling sum test

S = 7

calcrollsum = function(inputdata,size=5){
  require(zoo)
  rollapply(inputdata,size,sum,na.rm=TRUE)
}
ptm = proc.time()
preciprollsum = apply(precip,c(1,2),calcrollsum,size=S)
ptmend=  proc.time()-ptm

preciprollsum2 = array(NA,dim=c(length(lon),length(lat),dim(preciprollsum)[1]))

for(m in 1:dim(preciprollsum)[1]){
  preciprollsum2[,,m] = ifelse(is.na(precip[,,1])==FALSE,preciprollsum[m,,],NA)
}

output4 = whiplashcalc(preciprollsum2)
output5 = whiplashcalc(preciprollsum2,spread=S)



#####
# plotting
shapefile = "/home/woot0002/shapefiles/WBD_HU8_SanSabaLlano"

message("opening URG shapefile")
test = readShapePoly(shapefile)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

##
# Monthly based
DPcolorbar = colorramp(c(output$DPcount,output$PDcount),colorchoice="whitetogreen",Blimit=20,type="raw",use_fixed_scale = FALSE)

testsfc = list(x=lon,y=lat,z=output$DPcount)
surface(testsfc,type="I",main="Drought-Pluvial Transitions - monthly",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar[[1]],col=DPcolorbar[[3]],breaks=DPcolorbar[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

testsfc = list(x=lon,y=lat,z=output$PDcount)
surface(testsfc,type="I",main="Pluvial-Drought Transitions - monthly",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar[[1]],col=DPcolorbar[[3]],breaks=DPcolorbar[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)


##
# 30-day (spread=1,max overlap in moving average) based
DPcolorbar2 = colorramp(c(output2$DPcount,output2$PDcount),colorchoice="whitetogreen",Blimit=20,type="raw",use_fixed_scale = FALSE)

testsfc = list(x=lon,y=lat,z=output2$DPcount)
surface(testsfc,type="I",main="Drought-Pluvial Transitions - 30 day (spread=1)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar2[[1]],col=DPcolorbar2[[3]],breaks=DPcolorbar2[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

testsfc = list(x=lon,y=lat,z=output2$PDcount)
surface(testsfc,type="I",main="Pluvial-Drought Transitions - 30 day (spread=1)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar2[[1]],col=DPcolorbar2[[3]],breaks=DPcolorbar2[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

##
# 30-day (spread=30,max overlap in moving average) based
DPcolorbar3 = colorramp(c(output3$DPcount,output3$PDcount),colorchoice="whitetogreen",Blimit=20,type="raw",use_fixed_scale = FALSE)

testsfc = list(x=lon,y=lat,z=output3$DPcount)
surface(testsfc,type="I",main="Drought-Pluvial Transitions - 30 day (spread=30)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar3[[1]],col=DPcolorbar3[[3]],breaks=DPcolorbar3[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

testsfc = list(x=lon,y=lat,z=output3$PDcount)
surface(testsfc,type="I",main="Pluvial-Drought Transitions - 30 day (spread=30)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar3[[1]],col=DPcolorbar3[[3]],breaks=DPcolorbar3[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

##
# 7-day (spread=1,max overlap in moving average) based
DPcolorbar4 = colorramp(c(output4$DPcount,output4$PDcount),colorchoice="whitetogreen",Blimit=20,type="raw",use_fixed_scale = FALSE)

testsfc = list(x=lon,y=lat,z=output4$DPcount)
surface(testsfc,type="I",main="Drought-Pluvial Transitions - 7 day (spread=1)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar4[[1]],col=DPcolorbar4[[3]],breaks=DPcolorbar4[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

testsfc = list(x=lon,y=lat,z=output4$PDcount)
surface(testsfc,type="I",main="Pluvial-Drought Transitions - 7 day (spread=1)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar4[[1]],col=DPcolorbar4[[3]],breaks=DPcolorbar4[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

##
# 7-day (spread=7,min overlap in moving average) based
DPcolorbar5 = colorramp(c(output5$DPcount,output5$PDcount),colorchoice="whitetogreen",Blimit=20,type="raw",use_fixed_scale = FALSE)

testsfc = list(x=lon,y=lat,z=output5$DPcount)
surface(testsfc,type="I",main="Drought-Pluvial Transitions - 7 day (spread=7)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar5[[1]],col=DPcolorbar5[[3]],breaks=DPcolorbar5[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)

testsfc = list(x=lon,y=lat,z=output5$PDcount)
surface(testsfc,type="I",main="Pluvial-Drought Transitions - 7 day (spread=7)",xlab="Longtitude",ylab="Latitude",zlim=DPcolorbar5[[1]],col=DPcolorbar5[[3]],breaks=DPcolorbar5[[2]])
plot(test,add=TRUE)
map("state",add=TRUE)
