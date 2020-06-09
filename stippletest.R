###################
#
# R stippling test


library(ncdf4)
library(maps)
library(fields)
library(sp)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

# set basic arguments
varname = "pr" # these are all required arguments for step 1
steps = c(1,2,3,4,5,6,7) # others can be set based upon what varname is.
difftype = "absolute"
futureperiod = c(2070,2099)
varunits = "mm"
changeunits = "mm"
seasonin = "ann"

# set regional calculation arguments 
lon = c(-101,-94)+360 # information needed for step 6 if regiontype = "box"
lat = c(33,35)
regiontype = "states_and_domain"
regionname = "oklahoma"

# set point calculation table information
locationtable = read.table(file="/data2/3to5/I35/scripts/location_name.csv",header=TRUE,sep=",",colClasses=c("character","numeric","numeric"))
locationtable[,3] = as.numeric(locationtable[,3])
locationtable[,2] = as.numeric(locationtable[,2])

# set map plot arguments
BINLIMIT=30
colorchoicediff = "browntogreen"
diffbartype = "difference"

# set ts pull specific arguments
useobs=  FALSE

# set format change arguments for ensemble means  
outfileformat = "GTiff" # file format for step 5

# If you aren't running step 1 you must supply the filename. If you aren't running step 2 you must also supply the filename.
step1_filename = NA # if these are NA and you are not running step1 or step2, then other options that rely on these will break
step2_filename = NA

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

############

if(6 %in% steps | 3 %in% steps | 7 %in% steps | 2 %in% steps){
  if(1 %in% steps == FALSE & is.na(step1_filename)==TRUE){
    stop("To Calculate steps 2,3,6, and 7, you must calculate step1 or provide the file name and path to the Individual members file",.call=TRUE)
  }
}

if(4 %in% steps | 5 %in% steps){
  if(2 %in% steps == FALSE & is.na(step2_filename)==TRUE){
    stop("To Calculate steps 4 and 5, you must calculate step2 or provide the file name and path to the Ensemble means file",.call=TRUE)
  }
}

varin = varname
if(varname=="tmax85" | varname=="tmax90" | varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="prcptot" | varname=="prmean" | varname=="prmed" | varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd" | varname=="R99" | varname=="R95" | varname=="R90"  | varname=="R99freq" | varname=="R95freq" | varname=="R90freq") varin="pr"

TC = FALSE 
TH = NA
cond=NA

if(varname == "tasmax"){
  TC=FALSE
  appfunc = "mean"
}

if(varname == "tasmin"){
  TC=FALSE
  appfunc = "mean"
}

if(varname == "pr"){
  TC=FALSE
  appfunc = "sum"
}

if(varname == "heatwaves"){
  TC=FALSE
  appfunc = "heatwaves"
  varin = "tasmax"
}

if(varname == "gsl"){
  TC=FALSE
  appfunc = "growing_season_length"
  varin = "tasmax"
}

if(varname == "R90"){
  TC=FALSE
  appfunc = "R90"
  varin = "pr"
}

if(varname == "R95"){
  TC=FALSE
  appfunc = "R95"
  varin = "pr"
}

if(varname == "R99"){
  TC=FALSE
  appfunc = "R99"
  varin = "pr"
}


if(varname == "R90freq"){
  TC=FALSE
  appfunc = "R90freq"
  varin = "pr"
}

if(varname == "R95freq"){
  TC=FALSE
  appfunc = "R95freq"
  varin = "pr"
}

if(varname == "R99freq"){
  TC=FALSE
  appfunc = "R99freq"
  varin = "pr"
}

if(varname == "tmax85"){
  TC=TRUE
  TH = 302.59
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmax90"){
  TC=TRUE
  TH = 305.372
  cond = "gte"
  appfunc = "sum"
}


if(varname == "tmax95"){
  TC=TRUE
  TH = 308.15
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmax100"){
  TC=TRUE
  TH = 310.928
  cond = "gte"
  appfunc = "sum"
}

if(varname == "tmin32"){
  TC=TRUE
  TH = 273.15
  cond = "lte"
  appfunc = "sum"
}

if(varname == "tmin28"){
  TC=TRUE
  TH = 270.928
  cond = "lte"
  appfunc = "sum"
}

if(varname == "frd"){
  TC=FALSE
  appfunc = "lastfreeze"
}

if(varname == "prcptot"){
  TC=FALSE
  appfunc = "prcptot"
}

if(varname == "prmean"){
  TC=FALSE
  appfunc = "prmean"
}

if(varname == "prmed"){
  TC=FALSE
  appfunc = "prmed"
}

if(varname == "mdrn"){
  TC=TRUE
  TH = 0.254
  cond = "gte"
  appfunc = "sum"
}

if(varname == "r1mm"){
  TC=TRUE
  TH = 1
  cond = "gte"
  appfunc = "sum"
}

if(varname == "pr25"){
  TC=TRUE
  TH = 25.4
  cond = "gte"
  appfunc = "sum"
}

if(varname == "pr50"){
  TC=TRUE
  TH = 50.8
  cond = "gte"
  appfunc = "sum"
}

if(varname == "rx1day"){
  TC=FALSE
  appfunc = "max"
}

if(varname == "rx5day"){
  TC=FALSE
  appfunc = "rx5day"
}

if(varname == "cdd"){
  TC=FALSE
  appfunc = "maxdryspell"
}

if(varname == "cwd"){
  TC=FALSE
  appfunc = "maxwetspell"
}

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

###########

test = nc_open("/data2/3to5/I35/all_mems/pr_allmem_percent_2070-2099_ann.nc")
memproj = ncvar_get(test,"projmeandiff")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

test = nc_open("/data2/3to5/I35/ens_means/pr_ensmean_percent_2070-2099_ann.nc")
meanproj85 = ncvar_get(test,"projmeandiff_rcp85")
meanproj45 = ncvar_get(test,"projmeandiff_rcp45")
meanproj26 = ncvar_get(test,"projmeandiff_rcp26")
nc_close(test)

rcp26idx = which(projfilebreakdown$scen=="rcp26")
rcp45idx = which(projfilebreakdown$scen=="rcp45")
rcp85idx = which(projfilebreakdown$scen=="rcp85")

####
# v1 calcs

tmppos = ifelse(memproj[,,rcp26idx]>=0,1,0)
tmpneg = ifelse(memproj[,,rcp26idx]<0,1,0)
agreepos26 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
agreeneg26 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)

tmppos = ifelse(memproj[,,rcp45idx]>=0,1,0)
tmpneg = ifelse(memproj[,,rcp45idx]<0,1,0)
agreepos45 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
agreeneg45 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)

tmppos = ifelse(memproj[,,rcp85idx]>=0,1,0)
tmpneg = ifelse(memproj[,,rcp85idx]<0,1,0)
agreepos85 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
agreeneg85 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)

agree85v1 = ifelse(agreepos85>18 | agreeneg85>18,1,0)
agree45v1 = ifelse(agreepos45>18 | agreeneg45>18,1,0)
agree26v1 = ifelse(agreepos26>18 | agreeneg26>18,1,0)

####
# v2 calcs

tmppos = ifelse(memproj[,,rcp26idx]>5,1,0)
tmpmid = ifelse(memproj[,,rcp26idx]<=5 & memproj[,,rcp26idx]>= -5,1,0)
tmpneg = ifelse(memproj[,,rcp26idx]< -5,1,0)
agreepos26 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
agreemid26 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
agreeneg26 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)

tmppos = ifelse(memproj[,,rcp45idx]>5,1,0)
tmpmid = ifelse(memproj[,,rcp45idx]<=5 & memproj[,,rcp45idx]>= -5,1,0)
tmpneg = ifelse(memproj[,,rcp45idx]< -5,1,0)
agreepos45 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
agreemid45 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
agreeneg45 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)

tmppos = ifelse(memproj[,,rcp85idx]>5,1,0)
tmpmid = ifelse(memproj[,,rcp85idx]<=5 & memproj[,,rcp85idx]>= -5,1,0)
tmpneg = ifelse(memproj[,,rcp85idx]< -5,1,0)
agreepos85 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
agreemid85 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
agreeneg85 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)

agree85v2 = ifelse(agreepos85>18 | agreeneg85>18 | agreemid85>18,1,0)
agree45v2 = ifelse(agreepos45>18 | agreeneg45>18 | agreemid45>18,1,0)
agree26v2 = ifelse(agreepos26>18 | agreeneg26>18 | agreemid26>18,1,0)

####
# v3 calcs

tmpsd = apply(memproj[,,rcp26idx],c(1,2),sd,na.rm=TRUE)
tmpmean = abs(meanproj26)
upperbd = tmpmean+2*tmpsd
lowerbd = tmpmean-2*tmpsd
agree26v3 = ifelse(lowerbd>0,1,0)

tmpsd = apply(memproj[,,rcp45idx],c(1,2),sd,na.rm=TRUE)
tmpmean = abs(meanproj45)
upperbd = tmpmean+2*tmpsd
lowerbd = tmpmean-2*tmpsd
agree45v3 = ifelse(lowerbd>0,1,0)

tmpsd = apply(memproj[,,rcp85idx],c(1,2),sd,na.rm=TRUE)
tmpmean = abs(meanproj85)
upperbd = tmpmean+2*tmpsd
lowerbd = tmpmean-2*tmpsd
agree85v3 = ifelse(lowerbd>0,1,0)

######
# plot maps

testsfc = list(x=lon-360,y=lat,z=agree85v1)
surface(testsfc,type="I")
map("state",add=TRUE)

testsfc = list(x=lon-360,y=lat,z=agree85v2)
surface(testsfc,type="I")
map("state",add=TRUE)

testsfc = list(x=lon-360,y=lat,z=agree85v3)
surface(testsfc,type="I")
map("state",add=TRUE)

agree85idx = which(agree85v1==1,arr.ind=TRUE)
agree45idx = which(agree45v1==1,arr.ind=TRUE)
agree26idx = which(agree26v1==1,arr.ind=TRUE)

diffcolorbar = colorramp(c(meanproj85,meanproj45,meanproj26),colorchoice="browntogreen",Blimit=15,use_fixed_scale = TRUE,fixed_scale=c(-20,20))

pdf("/home/woot0002/stippletestv1.pdf",height=6,width=6,onefile=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=meanproj26)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP2.6",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree26idx[,1]]-360, lat[agree26idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

testsfc1 = list(x=lon-360,y=lat,z=meanproj45)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree45idx[,1]]-360, lat[agree45idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

testsfc1 = list(x=lon-360,y=lat,z=meanproj85)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree85idx[,1]]-360, lat[agree85idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

dev.off()

###

agree85idx = which(agree85v2==1,arr.ind=TRUE)
agree45idx = which(agree45v2==1,arr.ind=TRUE)
agree26idx = which(agree26v2==1,arr.ind=TRUE)

diffcolorbar = colorramp(c(meanproj85,meanproj45,meanproj26),colorchoice="browntogreen",Blimit=15,use_fixed_scale = TRUE,fixed_scale=c(-20,20))

pdf("/home/woot0002/stippletestv2.pdf",height=6,width=6,onefile=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=meanproj26)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP2.6",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree26idx[,1]]-360, lat[agree26idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

testsfc1 = list(x=lon-360,y=lat,z=meanproj45)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree45idx[,1]]-360, lat[agree45idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

testsfc1 = list(x=lon-360,y=lat,z=meanproj85)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree85idx[,1]]-360, lat[agree85idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

dev.off()

###

agree85idx = which(agree85v3==1,arr.ind=TRUE)
agree45idx = which(agree45v3==1,arr.ind=TRUE)
agree26idx = which(agree26v3==1,arr.ind=TRUE)

diffcolorbar = colorramp(c(meanproj85,meanproj45,meanproj26),colorchoice="browntogreen",Blimit=15,use_fixed_scale = TRUE,fixed_scale=c(-20,20))

pdf("/home/woot0002/stippletestv3.pdf",height=6,width=6,onefile=TRUE)

testsfc1 = list(x=lon-360,y=lat,z=meanproj26)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP2.6",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree26idx[,1]]-360, lat[agree26idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

testsfc1 = list(x=lon-360,y=lat,z=meanproj45)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP4.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree45idx[,1]]-360, lat[agree45idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

testsfc1 = list(x=lon-360,y=lat,z=meanproj85)
surface(testsfc1,type="I",main="Projected Difference from Historical Climate\n Scen: RCP8.5",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
points(lon[agree85idx[,1]]-360, lat[agree85idx[,2]], cex = 0.05,col=rgb(0,0,0,alpha=0.4),pch=4)

dev.off()

