
library(iterators)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)
library(zoo)

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
calcrollsum = function(inputdata,size=5){
  require(zoo)
  rollapply(inputdata,width=size,FUN=sum,na.rm=TRUE)
}

source("/data2/3to5/I35/scripts/analysisfunctions.R")

##########
# Set input information

message("Setting inputs")

varname = "pr" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2070,2099)
S = 30 # rolling sum size
spr = 10 # spread to use for whiplash calcs

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

######
# historical calcs
histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
#h=1
years = 1981:2005

for(h in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[h])
  test = nc_open(histfilelist[h])
  precip = ncvar_get(test,"pr")*86400
  if(h==1){
    lat = ncvar_get(test,"lat")
    lon = ncvar_get(test,"lon")-360
    PDyearly_hist = DPyearly_hist = array(NA,dim=c(length(lon),length(lat),length(years),length(histfilelist)))
  }
  if(length(histdates)>dim(precip)[3]){
    histdatesin = histdates[-which(substr(histdates,6,10)=="02-29")]
  } else {
    histdatesin = histdates
  }
  
  times = ncvar_get(test,"time")
  nc_close(test)

  ptm = proc.time()
  message("Calculating rolling sum")
  preciprollsum2 = array(NA,dim=c(length(lon),length(lat),(dim(precip)[3]-(S-1))))
  for(r in 1:dim(precip)[1]){
    for(c in 1:dim(precip)[2]){
      if(is.na(precip[r,c,1])==FALSE){
        preciprollsum2[r,c,]= rollapply(precip[r,c,],width=S,FUN=sum,na.rm=TRUE)
      }
    }
    message("Finished rolling sum for row ",r)
  }
  
output = whiplashcalc(preciprollsum2,spread=spr)
histdates2 = histdatesin[30:length(histdatesin)]

for(y in 1:length(years)){
  idx = which(as.numeric(substr(histdates2,1,4))==years[y])
  tmp = apply(output$precip_PD[,,idx],c(1,2),sum,na.rm=TRUE)
  tmp = ifelse(is.na(output$precip_PD[,,1])==FALSE,tmp,NA)
  PDyearly_hist[,,y,h]=tmp
  tmp = apply(output$precip_DP[,,idx],c(1,2),sum,na.rm=TRUE)
  tmp = ifelse(is.na(output$precip_DP[,,1])==FALSE,tmp,NA)
  DPyearly_hist[,,y,h]=tmp
}


ptmend = proc.time()
message("Finished with file ",h," / ",length(histfilelist))
message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

#####
# projected calcs

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
years=2006:2099
for(i in 1:length(projfilelist)){
  ptm = proc.time()
  message("Starting work on file ",projfilelist[i])
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  }else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  #idx = which(as.numeric(substr(datesin,1,4))>=futureperiod[1] & as.numeric(substr(datesin,1,4))<=futureperiod[2])
  
  test = nc_open(projfilelist[i])
  #precip = ncvar_get(test,"pr",start=c(1,1,idx[1]),count=c(-1,-1,length(idx)))*86400
  precip = ncvar_get(test,"pr",start=c(1,1,1),count=c(-1,-1,-1))*86400
  if(i==1){
    PDyearly_proj = DPyearly_proj = array(NA,dim=c(length(lon),length(lat),length(years),length(projfilelist)))
  }
  times = ncvar_get(test,"time")
  nc_close(test)
  
  #ptm = proc.time()
  #preciprollsum = apply(precip,c(1,2),calcrollsum,size=S)
  #ptmend=  proc.time()-ptm
  
  preciprollsum2 = array(NA,dim=c(length(lon),length(lat),(dim(precip)[3]-(S-1))))
  for(r in 1:dim(precip)[1]){
    for(c in 1:dim(precip)[2]){
      if(is.na(precip[r,c,1])==FALSE){
        preciprollsum2[r,c,]= rollapply(precip[r,c,],width=S,FUN=sum,na.rm=TRUE)
      }
    }
    message("Finished rolling sum for row ",r)
  }
  
  output = whiplashcalc(preciprollsum2,spread=spr)
  dates2 = datesin[30:length(datesin)]
  
  for(y in 1:length(years)){
    idx = which(as.numeric(substr(histdates2,1,4))==years[y])
    tmp = apply(output$precip_PD[,,idx],c(1,2),sum,na.rm=TRUE)
    tmp = ifelse(is.na(output$precip_PD[,,1])==FALSE,tmp,NA)
    PDyearly_proj[,,y,i]=tmp
    tmp = apply(output$precip_DP[,,idx],c(1,2),sum,na.rm=TRUE)
    tmp = ifelse(is.na(output$precip_DP[,,1])==FALSE,tmp,NA)
    DPyearly_proj[,,y,i]=tmp
  }

  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}

######
# 1-d Write out netcdf files

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

# notes on the order of the ensemble members for the output netcdf file.

step1_filePDp = paste("PDp",S,"sp",spr,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_yearly.nc",sep="")
step1_fileDPp = paste("DPp",S,"sp",spr,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_yearly.nc",sep="")

dimX <- ncdim_def( "lon", "degrees_east", lon)
dimY <- ncdim_def( "lat", "degrees_north", lat)
dimTh <- ncdim_def("histyears","year",1981:2005)
dimTf <- ncdim_def("projyears","year",2006:2099)
dimF <- ncdim_def("ens_mem_proj",projnotes,1:nrow(projfilebreakdown))
dimH <- ncdim_def("ens_mem_hist",histnotes,1:nrow(histfilebreakdown))

# Make varables of various dimensionality, for illustration purposes
mv <- 1E20 # missing value to use

#######
# Create netcdf file (one per base variable)

#PDp
var1d <- ncvar_def("histmean","",longname="Historical Mean All Members", list(dimX,dimY,dimTh,dimH), mv )
var2d <- ncvar_def("projmean","",longname="Projected Mean All Members", list(dimX,dimY,dimTf,dimF), mv )
nc <- nc_create(paste("/home/woot0002/",step1_filePDp,sep="") ,  list(var1d,var2d) )
# Write some data to the file
ncvar_put(nc, var1d, PDyearly_hist) # no start or count: write all values\
ncvar_put(nc, var2d, PDyearly_proj)
# close ncdf
nc_close(nc)

#DPp
var1d <- ncvar_def("histmean","",longname="Historical Mean All Members", list(dimX,dimY,dimTh,dimH), mv )
var2d <- ncvar_def("projmean","",longname="Projected Mean All Members", list(dimX,dimY,dimTf,dimF), mv )
nc <- nc_create(paste("/home/woot0002/",step1_fileDPp,sep="") ,  list(var1d,var2d) )
# Write some data to the file
ncvar_put(nc, var1d, DPyearly_hist) # no start or count: write all values\
ncvar_put(nc, var2d, DPyearly_proj)
# close ncdf
nc_close(nc)
