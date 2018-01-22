###################
#
# METDATA data grabber for Task 4 work


# libraries and set up
library(ncdf4)
library(sp)
library(fields)

# set arguments
res = "50km"
location = "Bragg"

###
# distance function

distfunc = function(gridmeta,lat,lon){
dist = 3963*acos(sin( gridmeta$lat / 57.2958 ) * sin( lat / 57.2958 ) + cos( gridmeta$lat / 57.2958 ) * cos( lat /57.2958 ) * cos( ( lon / 57.2958) - ( gridmeta$lon / 57.2958 ) ))
distidx = which(dist==min(dist))
results = gridmeta[distidx,]
results
}

#############
# set file path

if(res=="4km") path="idaho_v2/"
if(res=="10km") path="idaho_v2/res10km/"
if(res=="15km") path="idaho_v2/res15km/"
if(res=="20km") path="idaho_v2/res20km/"
if(res=="30km") path="idaho_v2/res30km/"
if(res=="40km") path="idaho_v2/res40km/"
if(res=="50km") path="idaho_v2/res50km/"

#############
# set location coordinates

if(location=="MCBCL"){
latval=34.70734
lonval=-77.44516
}

if(location=="Bragg"){
latval=35.17083
lonval=-79.0145
}

#############
# get gridmeta for distance calcs from the first file, and get high temperature data

test = nc_open(paste(path,"tmmx.nc",sep=""))
lonname = test$var[[1]]$dim[[1]]$name
latname = test$var[[1]]$dim[[2]]$name
lon = ncvar_get(test,lonname)
lat = ncvar_get(test,latname)
time = ncvar_get(test,test$var[[1]]$dim[[3]]$name)
timeunits = test$var[[1]]$dim[[3]]$units
spunits = test$var[[1]]$dim[[1]]$units
missingval = test$var[[1]]$missval

R = 1:length(lon)
C = 1:length(lat)
r = rep(R,each=length(C))
c = rep(C,length(R))
lons = rep(lon,each=length(lat))
lats = rep(lat,length(lon))

gridmeta = data.frame(lons,lats,r,c)
names(gridmeta) = c("lon","lat","r","c")

inputpoint = distfunc(gridmeta,latval,lonval)

tmmxdat = ncvar_get(test,"air_temperature",start=c(inputpoint$r[1],inputpoint$c[1],1),count=c(1,1,-1))

nc_close(test)

#############
# get tmin and precipitation data

test = nc_open(paste(path,"tmmn.nc",sep=""))
tmmndat = ncvar_get(test,"air_temperature",start=c(inputpoint$r[1],inputpoint$c[1],1),count=c(1,1,-1))
nc_close(test)

test = nc_open(paste(path,"pr.nc",sep=""))
prdat = ncvar_get(test,"precipitation_amount",start=c(inputpoint$r[1],inputpoint$c[1],1),count=c(1,1,-1))
nc_close(test)

#############
# make dates vector

starttime = as.Date("1900-01-01")
dates = starttime+time
dates = as.character(dates)

############
# create time series data frame

tmmxdat = tmmxdat-273.15
tmmndat = tmmndat-273.15

tsdat = data.frame(dates,tmmxdat,tmmndat,prdat)
names(tsdat) = c("date","tmmx","tmmn","pr")

tsdat$tmean = (tsdat$tmmx+tsdat$tmmn)/2

####################
################
# Calculate monthly aggregations

###
# mean daily temperatures

meantmpdat = aggregate(tsdat[,c(2,3,5)],by=list(yearmon=substr(tsdat[,1],1,7)),mean,na.rm=TRUE)
names(meantmpdat) = c("yearmon","meantmmx","meantmmn","meantmean")

###
# variance of daily max temperature

vartmmxdat = aggregate(tsdat[,2],by=list(yearmon=substr(tsdat[,1],1,7)),var,na.rm=TRUE)
names(vartmmxdat) = c("yearmon","vartmmx")

###
# number of days with tmin below thresholds

dtmmn_5 = ifelse(tsdat$tmmn<5,1,0)
dtmmn_n8 = ifelse(tsdat$tmmn< -8,1,0)
daystmmndat = data.frame(dates,dtmmn_5,dtmmn_n8)

tmmndaysumdat = aggregate(daystmmndat[,2:3],by=list(yearmon=substr(daystmmndat[,1],1,7)),sum,na.rm=TRUE)

###
# number of days with tmax above thresholds

dtmmx_25 = ifelse(tsdat$tmmx>25,1,0)
dtmmx_35 = ifelse(tsdat$tmmx>35,1,0)
daystmmxdat = data.frame(dates,dtmmx_25,dtmmx_35)

tmmxdaysumdat = aggregate(daystmmxdat[,2:3],by=list(yearmon=substr(daystmmxdat[,1],1,7)),sum,na.rm=TRUE)

###
# sum of daily total precipitation

prsumdat = aggregate(tsdat[,4],by=list(yearmon=substr(tsdat[,1],1,7)),sum,na.rm=TRUE)
names(prsumdat) = c("yearmon","prsum")

###
# variance of daily total precipitation

prvardat = aggregate(tsdat[,4],by=list(yearmon=substr(tsdat[,1],1,7)),var,na.rm=TRUE)
names(prvardat) = c("yearmon","varpr")

##################
# combine output and create the file

outputdata = merge(meantmpdat,vartmmxdat,by="yearmon")
outputdata = merge(outputdata,tmmxdaysumdat,by="yearmon")
outputdata = merge(outputdata,tmmndaysumdat,by="yearmon")
outputdata = merge(outputdata,prsumdat,by="yearmon")
outputdata = merge(outputdata,prvardat,by="yearmon")

outputdata$year = substr(outputdata$yearmon,1,4)
outputdata$mon = substr(outputdata$yearmon,6,7)

outputdata = outputdata[,-1]

w <- reshape(outputdata, 
  timevar = "mon",
  idvar = "year",
  direction = "wide")

filename = paste("idaho_v2/",location,"_",res,"_monthly.csv",sep="")
write.table(w,filename,sep=",",row.names=FALSE)








