###################
#
# R stippling test

library(ncdf4)
library(maps)
library(fields)
library(sp)

source("/data2/3to5/I35/scripts/analysisfunctions.R")

# set basic arguments
varname = "tmin28" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2070,2099)
colorchoice = "redtoblue"
colorlimits = c(-90,90)
BINS = 20
USEFS = FALSE

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

###########

test = nc_open(paste("/data2/3to5/I35/all_mems/",varname,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_ann.nc",sep=""))
memhist = ncvar_get(test,"histmean")
memproj = ncvar_get(test,"projmeandiff")
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
nc_close(test)

test = nc_open(paste("/data2/3to5/I35/ens_means/",varname,"_ensmean_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_ann.nc",sep=""))
meanproj85 = ncvar_get(test,"projmeandiff_rcp85")
meanproj45 = ncvar_get(test,"projmeandiff_rcp45")
meanproj26 = ncvar_get(test,"projmeandiff_rcp26")
nc_close(test)

rcp26idx = which(projfilebreakdown$scen=="rcp26")
rcp45idx = which(projfilebreakdown$scen=="rcp45")
rcp85idx = which(projfilebreakdown$scen=="rcp85")

####
# v2 calcs

if(difftype=="percent"){
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
} 

if(difftype=="absolute"){
  
  chgper = (memproj[,,rcp26idx]/memhist)*100
  
  tmppos = ifelse(chgper>5,1,0)
  tmpmid = ifelse(chgper<=5 & chgper>= -5,1,0)
  tmpneg = ifelse(chgper< -5,1,0)
  agreepos26 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
  agreemid26 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
  agreeneg26 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)
  
  chgper = (memproj[,,rcp45idx]/memhist)*100
  
  tmppos = ifelse(chgper>5,1,0)
  tmpmid = ifelse(chgper<=5 & chgper>= -5,1,0)
  tmpneg = ifelse(chgper< -5,1,0)
  agreepos45 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
  agreemid45 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
  agreeneg45 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)
  
  chgper = (memproj[,,rcp85idx]/memhist)*100
  
  tmppos = ifelse(chgper>5,1,0)
  tmpmid = ifelse(chgper<=5 & chgper>= -5,1,0)
  tmpneg = ifelse(chgper< -5,1,0)
  agreepos85 = apply(tmppos,c(1,2),sum,na.rm=TRUE)
  agreemid85 = apply(tmpmid,c(1,2),sum,na.rm=TRUE)
  agreeneg85 = apply(tmpneg,c(1,2),sum,na.rm=TRUE)
} 

agree85v2 = ifelse(agreepos85>18 | agreeneg85>18 | agreemid85>18,1,0)
agree45v2 = ifelse(agreepos45>18 | agreeneg45>18 | agreemid45>18,1,0)
agree26v2 = ifelse(agreepos26>18 | agreeneg26>18 | agreemid26>18,1,0)

######
# plot maps

testsfc = list(x=lon-360,y=lat,z=agree85v2)
surface(testsfc,type="I")
map("state",add=TRUE)

###

agree85idx = which(agree85v2==1,arr.ind=TRUE)
agree45idx = which(agree45v2==1,arr.ind=TRUE)
agree26idx = which(agree26v2==1,arr.ind=TRUE)

diffcolorbar = colorramp(c(meanproj85,meanproj45,meanproj26),colorchoice=colorchoice,Blimit=BINS,use_fixed_scale = USEFS,fixed_scale=colorlimits)

pdf(paste("/home/woot0002/stippletest_",varname,".pdf",sep=""),height=6,width=6,onefile=TRUE)

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


