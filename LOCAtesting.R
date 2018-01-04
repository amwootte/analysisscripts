####################
#
# Rscript for Schooner mount test

library(ncdf4)
library(maps)
library(fields)

inputdir = "/schooner/"
outputdir = "/home/woot0002/LOCAtest/"

futfiles = system(paste("ls ",inputdir,"LOCA/future/pr_ACCESS1-3*rcp45*.nc",sep=""),intern=TRUE)
futdat = do.call("rbind",strsplit(futfiles,"_"))
futyears = do.call("rbind",strsplit(futdat[,4],"-"))
futyears[,3] = substr(futyears[,3],1,4)
futyears[,2:3]=as.numeric(futyears[,2:3])

####
# Calculate yearly totals

for(i in 1:length(futfiles)){

  years = seq(as.numeric(futyears[i,2]),as.numeric(futyears[i,3]))
  dates = seq(as.Date(paste(years[1],"-01-01",sep="")),as.Date(paste(years[length(years)],"-12-31",sep="")),by="day")
    
  test = nc_open(futfiles[i])
  prdata = ncvar_get(test,paste("pr_",futdat[i,2],"_",futdat[i,3],"_",futyears[i,1],sep=""))*86400
  if(i==1){
    lat = ncvar_get(test,"lat")
    lon = ncvar_get(test,"lon")
    yearlydat = array(NA,dim=c(length(lon),length(lat),length(2006:2100)))
    yearlen = 0
  }
  
  nc_close(test)
  
  for(y in 1:length(years)){
    yearidx = which(as.numeric(substr(dates,1,4))==years[y])
    if(max(yearidx)>dim(prdata)[3]){
      endidx = which(yearidx==dim(prdata)[3])
      yearidx = yearidx[1:endidx]
    }
    tmp = apply(prdata[,,yearidx],c(1,2),sum,na.rm=TRUE)
    yearlydat[,,yearlen+y] = ifelse(is.na(prdata[,,1])==FALSE,tmp,NA)
  }
  yearlen = yearlen+length(years)
  message("Finished yearly sums for futfile ",i)
}

years = 2006:2100
yearidx = which(years>=2041 & years<=2070)
climo = apply(yearlydat[,,yearidx],c(1,2),mean,na.rm=TRUE)

rm(yearlydat)
gc()

############

histfile = paste(inputdir,"LOCA/historical/pr_ACCESS1-3_r1i1p1_historical.nc",sep="")
test= nc_open(histfile)

times = ncvar_get(test,"time")
startdate = as.Date("1900-01-01")
dates = startdate+times

years = 1976:2005

yearlydat = array(NA,dim=c(length(lon),length(lat),length(years)))

for(y in 1:length(years)){
  yearidx = which(as.numeric(substr(dates,1,4))==years[y])
  prdata = ncvar_get(test,"pr_ACCESS1-3_r1i1p1_historical",start=c(1,1,yearidx[1]),count=c(-1,-1,length(yearidx)))
  tmp = apply(prdata,c(1,2),sum,na.rm=TRUE)
  yearlydat[,,y] = ifelse(is.na(prdata[,,1])==FALSE,tmp,NA)
  message("Finished with calc for year ",y," / ",length(years))
}

nc_close(test)

histclimo = apply(yearlydat,c(1,2),mean,na.rm=TRUE)

###############

diffs = ((climo-histclimo)/histclimo)*100

################
# make a plot suite

pdf(paste(outputdir,"LOCAplottest.pdf",sep=""),onefile=TRUE,width=8,height=6)

testsfc1 = list(x=lon-360,y=lat,z=histclimo)
testsfc2 = list(x=lon-360,y=lat,z=climo)
testsfc3 = list(x=lon-360,y=lat,z=diffs)

surface(testsfc1,type="I",xlab="Longitude",ylab="Latitude",main="Historical (1976-2005) Average Annual Precipitation")
map("world",add=TRUE)
surface(testsfc2,type="I",xlab="Longitude",ylab="Latitude",main="Future (2041-2070) Average Annual Precipitation")
map("world",add=TRUE)
surface(testsfc3,type="I",xlab="Longitude",ylab="Latitude",main="2041-2070 Projected Change in Precipitation")
map("world",add=TRUE)

dev.off()




  







