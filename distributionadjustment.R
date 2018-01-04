####################
#
# 30% precipitation adjustment calculation
# for NAF activity pseudo-data
#
###################################

library(ncdf4)
library(maps)
library(fields)

###########
# 1. Data Gather and conversion
setwd("/home/woot0002")
filelist = system("ls /data2/3to5/I35/pr/EDQM/ratio/*.nc",intern=T)

years = 2041:2070
prdat = list()

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
leapdays = which(substr(dates,6,10)=="02-29")
dates = dates[-leapdays]

timeidx = which(as.numeric(substr(dates,1,4))>=years[1] & as.numeric(substr(dates,1,4))<=years[2])

for(i in 1:length(filelist)){
  ptm = proc.time()
  test = nc_open(filelist[i])
  tempdata = ncvar_get(test,"pr",start = c(1,1,timeidx[1]),count=c(-1,-1,length(timeidx)))
  
  if(i==1){
    lat = ncvar_get(test,"lat")
    lon = ncvar_get(test,"lon")-360
    times = ncvar_get(test,"time")
    domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
  }
  nc_close(test)
  
  prdat[[i]]=tempdata*86400
  
  rm(tempdata)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(filelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}


distadj = function(inputdata,type="high"){
  
    if(all(is.na(inputdata)==TRUE)==FALSE){
      test=fitData(inputdata[which(inputdata>=0.254)],fit=c("exponential"),sample=1)
      meanc = 1/as.numeric(test[2,2])
      if(type=="high") test2=qexp(pexp(inputdata,rate=1/meanc),rate=1/(meanc*1.3));
      if(type=="low") test2=qexp(pexp(inputdata,rate=1/meanc),rate=1/(meanc*0.7))
      test2
    } else {
      inputdata
    }
  
  }

testarray = prdat[[1]]

for(r in 1:dim(prdat[[1]])[1]){
  for(c in 1:dim(prdat[[1]])[2]){
    testarray[r,c,]=distadj(prdat[[1]][r,c,],type="high")
    #test=apply(prdat[[1]],c(1,2),distadj,type="high")
    message("Finished row ",r," and col ",c)
  }
}
  
  
  



meanc = 1/as.numeric(test[2,2])
x=seq(0,100,by=0.0001)
plot(dexp(x,rate=1/meanc),type="l")
lines(dexp(x,rate=1/(meanc*1.3)),type="l",col="red")
lines(dexp(x,rate=1/(meanc*0.7)),type="l",col="blue")

pexp(0,rate=1/meanc)

test2=qexp(pexp(prdat[[1]][95,70,],rate=1/meanc),rate=1/(meanc*1.3))
test3=qexp(pexp(prdat[[1]][95,70,],rate=1/meanc),rate=1/(meanc*0.7))

lambda = 1/mean(prdat[[1]][95,70,])


obsecdf = ecdf(prdat[[1]][95,70,])

obsquant = quantile(prdat[[1]][95,70,],probs=seq(0.0,1,by=0.01),na.rm=TRUE)








