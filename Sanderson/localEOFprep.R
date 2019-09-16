###########################
#
# Sanderson work and EOF analysis - local
#
###########################
#setwd("/data2/3to5/I35/pr/EDQM")

filespr = system("ls /data2/3to5/I35/pr/EDQM/*.nc",intern=TRUE)
filestmax = system("ls /data2/3to5/I35/tasmax/EDQM/*.nc",intern=TRUE)
filestmin = system("ls /data2/3to5/I35/tasmin/EDQM/*.nc",intern=TRUE)

files = c(filespr,filestmax,filestmin)

scenin="rcp85"

#####################

library(ncdf4)
library(sp)
library(fields)
library(smacof)

filebreakdown = do.call(rbind,strsplit(files,"/",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,7],"_",fixed=TRUE))
filebreakdown3 = do.call(rbind,strsplit(filebreakdown2[,3],"-",fixed=TRUE))

filepath = paste("",filebreakdown[,2],filebreakdown[,3],filebreakdown[,4],filebreakdown[,5],filebreakdown[,6],sep="/")
filebreakdown4 = data.frame(filepath,filebreakdown2[,1:2],filebreakdown3,filebreakdown2[,4:5])

GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown4$GCM = rep(GCM,3)
filebreakdown4$obs = rep(c("Daymet","Livneh","PRISM"),27)
filebreakdown4 = filebreakdown4[,-4]
names(filebreakdown4)=c("filepath","var","tempres","DS","code","scen","experiment","GCM","obs")
filebreakdown4$filename = filebreakdown[,7]

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
dates = dates[-which(substr(dates,6,10)=="02-29")]

earlyidx = which(as.numeric(substr(dates,1,4))>=2006 & as.numeric(substr(dates,1,4))<=2015) 
lateidx = which(as.numeric(substr(dates,1,4))>=2041 & as.numeric(substr(dates,1,4))<=2070)

#####################
# Gather all data - GCM netcdfs

prearlydatalist = list()
tasmaxearlydatalist = list()
tasminearlydatalist = list()

prlatedatalist = list()
tasmaxlatedatalist = list()
tasminlatedatalist = list()

filebreakdownin = subset(filebreakdown4,scen==scenin)

for(f in 1:nrow(filebreakdownin)){
  ptm1 = proc.time()
  message("Data grab began at ",Sys.time()," for file: ",filebreakdownin$filename[f])
  
  var = filebreakdownin$var[f]
  
  test = nc_open(paste(filebreakdownin$filepath[f],filebreakdownin$filename[f],sep="/"))
  tempdataearly = ncvar_get(test,as.character(var),start=c(1,1,earlyidx[1]),count=c(-1,-1,length(earlyidx)))
  tempdatalate = ncvar_get(test,as.character(var),start=c(1,1,lateidx[1]),count=c(-1,-1,length(lateidx)))
  
  if(f==1){
    lon = ncvar_get(test,"lon")
    lat = ncvar_get(test,"lat")
  }
  
  nc_close(test)
  
  tempearly = template = array(NA,dim=c(length(lon),length(lat),12))
  
  datese = dates[earlyidx]
  datesl = dates[lateidx]
  for(m in 1:12){
    etmpidx = which(as.numeric(substr(datese,6,7))==m)
    ltmpidx = which(as.numeric(substr(datesl,6,7))==m)
    tempearly[,,m] = apply(tempdataearly[,,etmpidx],c(1,2),mean,na.rm=TRUE)
    template[,,m] = apply(tempdatalate[,,ltmpidx],c(1,2),mean,na.rm=TRUE)
  }
  if(f<=9){
    prearlydatalist[[f]] = tempearly
    prlatedatalist[[f]] = template
  } 
  if(f>9 & f<=18){
    tasmaxearlydatalist[[(f-9)]] = tempearly
    tasmaxlatedatalist[[(f-9)]] = template
  } 
  if(f>18 & f<=27){
    tasminearlydatalist[[(f-18)]] = tempearly
    tasminlatedatalist[[f-18]] = template
  } 
  message("Data grab ends at ",Sys.time())
}

rm(tempdatalate)
rm(tempdataearly)
rm(tempearly)
rm(template)

##########
# Add in MACA data for rcp45 and rcp85

if(scenin!="rcp26"){
  filesprmaca = system(paste("ls /home/woot0002/maca*pr*",scenin,"*.nc",sep=""),intern=TRUE)
  filestmaxmaca = system(paste("ls /home/woot0002/maca*tasmax*",scenin,"*.nc",sep=""),intern=TRUE)
  filestminmaca = system(paste("ls /home/woot0002/maca*tasmin*",scenin,"*.nc",sep=""),intern=TRUE)
  
  GCMsin = c("CCSM4","MIROC5")
  
  for(i in 1:2){
    I=i+9
    filesprin = filesprmaca[grep(GCMsin[i],filesprmaca)]
    filestmaxin = filestmaxmaca[grep(GCMsin[i],filestmaxmaca)]
    filestminin = filestminmaca[grep(GCMsin[i],filestminmaca)]
    
    test = nc_open(filesprin[1])
    prearlydatalist[[I]]=ncvar_get(test,"pr")
    nc_close(test)
    
    test = nc_open(filesprin[2])
    prlatedatalist[[I]]=ncvar_get(test,"pr")
    nc_close(test)
    
    test = nc_open(filestmaxin[1])
    tasmaxearlydatalist[[I]]=ncvar_get(test,"tasmax")
    nc_close(test)
    
    test = nc_open(filestmaxin[2])
    tasmaxlatedatalist[[I]]=ncvar_get(test,"tasmax")
    nc_close(test)
    
    test = nc_open(filestminin[1])
    tasminearlydatalist[[I]]=ncvar_get(test,"tasmin")
    nc_close(test)
    
    test = nc_open(filestminin[2])
    tasminlatedatalist[[I]]=ncvar_get(test,"tasmin")
    nc_close(test)
  }
}

#####################
# Create area weighting matrix considering latitude

latweights = cos((lat*pi)/180)
latweightmat = matrix(rep(latweights,each=length(lon)),nrow=length(lon),ncol=length(latweights))
latweightvector  = as.vector(latweightmat)

#####################
# Weight function for later

weightfunc = function(datafile,latweightmat){
  results=array(data=NA,dim=c(nrow(datafile),ncol(datafile),12))
  for(i in 1:12){
    results[,,i]=datafile[,,i]*latweightmat
  }
  return(results)
}

#####################
# Apply area weighting to 3^5 data

prearlydatalist = lapply(prearlydatalist,"*",86400) # convert 3^5 to mm/day
prlatedatalist = lapply(prlatedatalist,"*",86400) # convert 3^5 to mm/day

prearlydatalist[[10]]=prearlydatalist[[10]]/86400 # MACA is already in mm/day so this makes sure it matches
prearlydatalist[[11]]=prearlydatalist[[11]]/86400 # MACA is already in mm/day so this makes sure it matches
prlatedatalist[[10]]=prlatedatalist[[10]]/86400 # MACA is already in mm/day so this makes sure it matches
prlatedatalist[[11]]=prlatedatalist[[11]]/86400 # MACA is already in mm/day so this makes sure it matches

tasmaxearlydatalist = lapply(tasmaxearlydatalist,"-",273.15) # convert 3^5 to C
tasmaxlatedatalist = lapply(tasmaxlatedatalist,"-",273.15) # convert 3^5 to C
tasminearlydatalist = lapply(tasminearlydatalist,"-",273.15) # convert 3^5 to C
tasminlatedatalist = lapply(tasminlatedatalist,"-",273.15) # convert 3^5 to C

prearlyweight = lapply(prearlydatalist,weightfunc,latweightmat)
tasmaxearlyweight = lapply(tasmaxearlydatalist,weightfunc,latweightmat)
tasminearlyweight = lapply(tasminearlydatalist,weightfunc,latweightmat)

prlateweight = lapply(prlatedatalist,weightfunc,latweightmat)
tasmaxlateweight = lapply(tasmaxlatedatalist,weightfunc,latweightmat)
tasminlateweight = lapply(tasminlatedatalist,weightfunc,latweightmat)

#####################
# Normalize all 3^5 data

prearlynorm = lapply(prearlyweight,"/",2.97)
tasmaxearlynorm = lapply(tasmaxearlyweight,"/",1.43)
tasminearlynorm = lapply(tasminearlyweight,"/",1.43)

prlatenorm = lapply(prlateweight,"/",2.97)
tasmaxlatenorm = lapply(tasmaxlateweight,"/",1.43)
tasminlatenorm = lapply(tasminlateweight,"/",1.43)

######################
# Combining all output

if(scenin=="rcp26"){
  namesplit = paste(filebreakdownin$GCM[1:9],filebreakdownin$DS[1:9],filebreakdownin$obs[1:9],sep="_")
} else {
  splittemp = paste(filebreakdownin$GCM[1:9],filebreakdownin$DS[1:9],filebreakdownin$obs[1:9],sep="_")
  namesplit = c(splittemp,"CCSM4_MACA_Livneh","MIROC5_MACA_Livneh")
}

prearlymat = do.call(rbind,prearlynorm)
tasmaxearlymat = do.call(rbind,tasmaxearlynorm)
tasminearlymat = do.call(rbind,tasminearlynorm)

prlatemat = do.call(rbind,prlatenorm)
tasmaxlatemat = do.call(rbind,tasmaxlatenorm)
tasminlatemat = do.call(rbind,tasminlatenorm)
# need it to check if everything is in the right order later

######################
# Calculate anomaly matrix

tasminearlyanom = tasminearlymat
tasmaxearlyanom = tasmaxearlymat
prearlyanom = prearlymat

tasminlateanom = tasminlatemat
tasmaxlateanom = tasmaxlatemat
prlateanom = prlatemat

tasminearlyens = apply(tasminearlymat,2,mean,na.rm=TRUE)
tasmaxearlyens = apply(tasmaxearlymat,2,mean,na.rm=TRUE)
prearlyens = apply(prearlymat,2,mean,na.rm=TRUE)

tasminlateens = apply(tasminlatemat,2,mean,na.rm=TRUE)
tasmaxlateens = apply(tasmaxlatemat,2,mean,na.rm=TRUE)
prlateens = apply(prlatemat,2,mean,na.rm=TRUE)

for(i in 1:nrow(tasmaxearlyanom)){
  tasmaxearlyanom[i,] = tasmaxearlymat[i,]-tasmaxearlyens
  tasminearlyanom[i,] = tasminearlymat[i,]-tasminearlyens
  prearlyanom[i,] = prearlymat[i,]-prearlyens
  tasmaxlateanom[i,] = tasmaxlatemat[i,]-tasmaxlateens
  tasminlateanom[i,] = tasminlatemat[i,]-tasminlateens
  prlateanom[i,] = prlatemat[i,]-prlateens
}

#####################
# Remove all columns where missing values are present in any row
# Result should be anomaly matrices for all variables with no missing values

###
# pr early

prearlycheck = ifelse(is.na(prearlyanom)==TRUE,1,0)
prearlycheck2 = apply(prearlycheck,2,sum)
prkeep = which(prearlycheck2==0)
prearlyanom = prearlyanom[,prkeep] # anomaly 

###
# pr late

prlatecheck = ifelse(is.na(prlateanom)==TRUE,1,0)
prlatecheck2 = apply(prlatecheck,2,sum)
prkeep = which(prlatecheck2==0)
prlateanom = prlateanom[,prkeep] # anomaly 

###
# tasmax early

tasmaxearlycheck = ifelse(is.na(tasmaxearlyanom)==TRUE,1,0)
tasmaxearlycheck2 = apply(tasmaxearlycheck,2,sum)
taskeep = which(tasmaxearlycheck2==0)
tasmaxearlyanom = tasmaxearlyanom[,taskeep] # anomaly 

###
# tasmax late

tasmaxlatecheck = ifelse(is.na(tasmaxlateanom)==TRUE,1,0)
tasmaxlatecheck2 = apply(tasmaxlatecheck,2,sum)
taskeep = which(tasmaxlatecheck2==0)
tasmaxlateanom = tasmaxlateanom[,taskeep] # anomaly 

###
# tasmin early

tasminearlycheck = ifelse(is.na(tasminearlyanom)==TRUE,1,0)
tasminearlycheck2 = apply(tasminearlycheck,2,sum)
taskeep = which(tasminearlycheck2==0)
tasminearlyanom = tasminearlyanom[,taskeep] # anomaly 

###
# tasmin late

tasminlatecheck = ifelse(is.na(tasminlateanom)==TRUE,1,0)
tasminlatecheck2 = apply(tasminlatecheck,2,sum)
taskeep = which(tasminlatecheck2==0)
tasminlateanom = tasminlateanom[,taskeep] # anomaly 


##################
# Combine into full matrix for the EOF analysis
#
# Combine in order - tas, rsut, rlut, psl

comboearlydat = cbind(prearlyanom,tasmaxearlyanom)
comboearlydat = cbind(comboearlydat,tasminearlyanom)

combolatedat = cbind(prlateanom,tasmaxlateanom)
combolatedat = cbind(combolatedat,tasminlateanom)

#######################
# Write out data so I can combine it elsewhere

models = namesplit
localearlyanom = comboearlydat
locallateanom = combolatedat
save(list = c("localearlyanom","locallateanom","models"),file=paste("/home/woot0002/",scenin,"_anomalymatrix.Rdata",sep=""))

















