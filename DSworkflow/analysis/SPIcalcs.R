####
# SPI Calcs - helper functions


spi.calc = function(filename,period=3,dates=seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")){
require(SPEI)
#filename = "/data2/3to5/I35/pr/DeltaSD/pr_day_I35prp1-DeltaSD-A10P01K00_historical_r6i1p1_I35Land_19810101-20051231.nc"
#period = 3
#dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

test = nc_open(filename)
prdata = ncvar_get(test,"pr")*86400
lon = ncvar_get(test,"lon")-360
lat = ncvar_get(test,"lat")
nc_close(test)

months = unique(substr(dates,1,7))
prmon = array(NA,dim=c(length(lon),length(lat),length(months)))

###
# step 1 convert precip to monthly

for(i in 1:length(months)){
  idx = which(substr(dates,1,7)==months[i])
  prmon[,,i] = apply(prdata[,,idx],c(1,2),sum,na.rm=TRUE)
  prmon[,,i] = ifelse(is.na(prdata[,,1])==FALSE,prmon[,,i],NA)
}

###
# step 2 calculate SPI

SPIout = prmon
#SPIout[r,c,]= 
for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    if(all(is.na(prdata[r,c,])==TRUE)==FALSE){
        SPIout[r,c,] = spi(prmon[r,c,],scale=period,distribution = 'Gamma',na.rm=TRUE)$fitted
        SPIout[r,c,] = ifelse(SPIout[r,c,]== -Inf | SPIout[r,c,]== Inf,NA,SPIout[r,c,])
    }
  }
}

list(lon,lat,months,SPIout)
#SPIclimo = apply(SPIout,c(1,2),mean,na.rm=TRUE)
#SPIclimo = ifelse(is.na(prdata[,,1])==FALSE,SPIclimo,NA)

}
