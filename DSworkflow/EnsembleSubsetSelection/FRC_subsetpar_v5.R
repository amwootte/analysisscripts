rm(list=ls(all=TRUE))
gc()


source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(mapdata)
library(maptools)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(ggplot2)
library(modi)
library(parallel)
library(foreach)
library(doParallel)
library(matrixStats)

#####
# Version of this script to run parallel but also one domain as a time.

#statelist = c("louisiana","oklahoma","north carolina","florida","kentucky","SGP-NCA","SE-NCA")
var = "pr"
enssize = 5
ensgroup = "EAA"
mask = "SGP-NCA"
climotype="annual"
histfile = paste("/home/woot0002/RCMES/",var,"_",ensgroup,"_climo_",mask,"_",climotype,".nc",sep="")
projfile = paste("/home/woot0002/RCMES/",var,"_",ensgroup,"_climo_ssp585_",mask,"_",climotype,".nc",sep="")

#####

nctest = nc_open(histfile)
lathist = ncvar_get(nctest,"lat")
lonhist = ncvar_get(nctest,"lon")

if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
  varnames = c()
  outputdata = array(NA,dim=c(length(lonhist),length(lathist),length(nctest$var)))
  for(i in 1:length(nctest$var)){
    varnames[i]=nctest$var[[i]]$name
    outputdata[,,i]=ncvar_get(nctest,varnames[i])
  }
} 

if(climotype=="monthly"){
  varnames = c()
  outputdata = array(NA,dim=c(length(lonhist),length(lathist),12,length(nctest$var)))
  for(i in 1:length(nctest$var)){
    varnames[i]=nctest$var[[i]]$name
    outputdata[,,,i]=ncvar_get(nctest,varnames[i])
  }
} 

varunits = nctest$var[[1]]$units

dimunits=c()
for(i in 1:length(nctest$dim)){
  dimunits[i]=nctest$dim[[i]]$units
}

nc_close(nctest)

varidx = which(varnames!="OBS" & varnames!="ENS")
varnames = varnames[varidx]
if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
outputdata = outputdata[,,varidx]
}
if(climotype=="monthly"){
  outputdata = outputdata[,,,varidx]
}

#####
# Future data

nctest = nc_open(projfile)
lathist = ncvar_get(nctest,"lat")
lonhist = ncvar_get(nctest,"lon")

if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
  varnames_proj = c()
  outputdata_proj = array(NA,dim=c(length(lonhist),length(lathist),length(nctest$var)))
  for(i in 1:length(nctest$var)){
    varnames_proj[i]=nctest$var[[i]]$name
    outputdata_proj[,,i]=ncvar_get(nctest,varnames_proj[i])
  }
} 

if(climotype=="monthly"){
  varnames_proj = c()
  outputdata_proj = array(NA,dim=c(length(lonhist),length(lathist),12,length(nctest$var)))
  for(i in 1:length(nctest$var)){
    varnames_proj[i]=nctest$var[[i]]$name
    outputdata_proj[,,,i]=ncvar_get(nctest,varnames_proj[i])
  }
} 


varunits = nctest$var[[1]]$units

dimunits=c()
for(i in 1:length(nctest$dim)){
  dimunits[i]=nctest$dim[[i]]$units
}

nc_close(nctest)

varidx = which(varnames_proj!="OBS" & varnames_proj!="ENS")
varnames_proj = varnames_proj[varidx]
if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
  outputdata_proj = outputdata_proj[,,varidx]
}
if(climotype=="monthly"){
  outputdata_proj = outputdata_proj[,,,varidx]
}

#####
# projected change calculation

projchange = outputdata_proj-outputdata

#####
# FRC calcs 

###
# Load things needed.

#load(paste("/home/woot0002/RCMES/",var,"_LOCA_RMSEdat_enssize",enssize,"_dom",mask,".Rdata",sep=""))
#combosused2 = combosused-1

combosused2 = combn(1:length(varnames_proj),enssize,simplify=TRUE)
idxs = 1:ncol(combosused2)

#testarray = apply(projchange,3,c)
#ensmax = rowMaxs(testarray,na.rm=TRUE)
#ensmin = rowMins(testarray,na.rm=TRUE)
if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
  ensmax = apply(projchange,c(1,2),max,na.rm=TRUE)
  ensmin = apply(projchange,c(1,2),min,na.rm=TRUE)
} 
if(climotype=="monthly"){
  ensmax = apply(projchange,c(1,2,3),max,na.rm=TRUE)
  ensmin = apply(projchange,c(1,2,3),min,na.rm=TRUE)
}
ensmax = ifelse(ensmax==-Inf,NA,ensmax)
ensmin = ifelse(ensmin==Inf,NA,ensmin)
ensrange = ensmax-ensmin

###
# FRC function

FRC = function(projchange,i,combosused2,ensrange,climotype="annual"){
  #i=idxs[n]
  #i = idxs[1]
  if(climotype=="annual" | climotype=="JAN" | climotype=="FEB" | climotype=="MAR" | climotype=="APR"  | climotype=="MAY" | climotype=="JUN" | climotype=="JUL" | climotype=="AUG" | climotype=="SEP"  | climotype=="OCT"  | climotype=="NOV"  | climotype=="DEC"){
    #testarray = apply(projchange[,,combosused2[,i]],3,c)  
    ensmaxsub = apply(projchange[,,combosused2[,i]],c(1,2),max,na.rm=TRUE)
    ensminsub = apply(projchange[,,combosused2[,i]],c(1,2),min,na.rm=TRUE)
  }
  if(climotype=="monthly"){
    #testarray = apply(projchange[,,,combosused2[,i]],3,c)
    ensmaxsub = apply(projchange[,,,combosused2[,i]],c(1,2,3),max,na.rm=TRUE)
    ensminsub = apply(projchange[,,,combosused2[,i]],c(1,2,3),min,na.rm=TRUE)
  }
  ensmaxsub = ifelse(ensmaxsub==-Inf,NA,ensmaxsub)
  ensminsub = ifelse(ensminsub==Inf,NA,ensminsub)
  #ensmaxsub = rowMaxs(testarray,na.rm=TRUE)
  #ensminsub = rowMins(testarray,na.rm=TRUE)
  #ensmaxsub = apply(projchange[,,combosused2[,i]],c(1,2),max,na.rm=TRUE)
  #ensminsub = apply(projchange[,,combosused2[,i]],c(1,2),min,na.rm=TRUE)
  #ensmaxsub = ifelse(ensmaxsub==-Inf,NA,ensmaxsub)
  #ensminsub = ifelse(ensminsub==Inf,NA,ensminsub)
  ensrangesub = ensmaxsub-ensminsub
  
  frcmat = ensrangesub/ensrange
  mean(frcmat,na.rm=TRUE) # return mean FRC value for domain
  
  #nfrcmat = (frcmat-mean(frcmat,na.rm=TRUE))/sd(frcmat,na.rm=TRUE)
  
  #frc = mean(frcmat,na.rm=TRUE)
  #nfrc = mean(nfrcmat,na.rm=TRUE)
  
  #c(frc,nfrc)
}

rm("nctest")
rm("dimunits")
rm("ensmax")
rm("ensmin")
rm("histfile")
rm("lathist") 
rm("lonhist")
rm("outputdata")
rm("outputdata_proj")
rm("projfile")

#system.time({
#FRC(projchange = projchange,i=idxs[10000],combosused2=combosused2,ensrange=ensrange)
#})

library(parallel)
#system.time({
  n_cores <- detectCores(logical=FALSE)
  cl <- makeCluster(n_cores-1, type = "FORK",methods=FALSE,useXDR=FALSE)
  registerDoParallel(cl)
  #FRClist <- foreach(n=1:length(idxs),.packages=c("matrixStats")) %dopar%  {
#filesout = c()
#for(n in 1:length(idxs)){
frc <- foreach(n=1:length(idxs),.combine="c") %dopar%  {
  FRC(projchange = projchange,i=idxs[n],combosused2=combosused2,ensrange=ensrange,climotype=climotype)
  #fout = paste("/home/woot0002/RCMES/scratch/FRC.group.",n,".",ensgroup,".",mask,".",enssize,".Rdata",sep="")
  #filesout = c(filesout,fout)
  #save("FRCres",file=fout)
  #message("Progress ",(n/length(idxs))*100,"% complete")
}    

## for error cases only!
#FRClist = list()
#for(n in 1:33649){
#  load(paste("/home/woot0002/RCMES/scratch/FRC.group.",n,".",ensgroup,".",mask,".",enssize,".Rdata",sep=""))
#  FRClist[[n]] = FRCres
#}
stopCluster(cl)

#frcout = do.call("c",frc)
#frc=frcout
nfrc = frc/sd(frc)

FRCs = data.frame(frc,nfrc)
#names(FRCs) = c("frc","nfrc")

fileout = paste("/home/woot0002/RCMES/",var,"_",ensgroup,"_FRCdat_enssize",enssize,"_dom",mask,"_",climotype,"_reduced.Rdata",sep="")
save(list=c("varnames","idxs","combosused2","FRCs","enssize"),file=fileout)

#save(list=c("idxs"),file=paste("RCMES/",var,"_LOCA_RMSEdat_enssize",enssize,"_dom",mask,"_idxs.Rdata",sep=""))
