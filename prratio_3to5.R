####################
#
# 3^5 Analysis - Ratio of days with Heavy Rain to Days with Rain
#
#  Want to do projected change for each thing, and by a couple groups.
#  Group 1: by emissions scenario
#  Group 2: by emissions scenario and GCM
#  Group 3: by emissions scenario and training dataset
#
#  two time periods
#  1. 2006-2015
#  2. 2041-2070
#
#
#  A couple steps for each
#  1. Gather and convert daily rain into annual total rainfall
#  2. Calculate climatology for the two time periods
#  3. plot maps for the projected change
#     a. absolute change
#     b. percent change
#     c. past vs. future maps
#
###################################

###########
# 1. Data Gather and conversion
setwd("/data2/3to5/I35/pr/EDQM")
filelist = system("ls *.nc",intern=T)

pryearlylist = list()

for(i in 1:length(filelist)){
  ptm = proc.time()
  test = nc_open(filelist[i])
  tempdata = ncvar_get(test,"pr")
  
  if(i==1){
    lat = ncvar_get(test,"lat")
    lon = ncvar_get(test,"lon")-360
    times = ncvar_get(test,"time")
    startdate = as.Date(substr(test$var[[5]]$dim[[3]]$units,12,21),"%Y-%m-%d")
    timeunits = test$var[[5]]$dim[[3]]$units
    dates = startdate+times
    years = 2006:2070
    domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
  }
  nc_close(test)
  
  thres = 25.4
  pryearlydat = array(NA,dim=c(length(lon),length(lat),length(years)))
  for(y in 1:length(years)){
    ptm2 = proc.time()
    yearidx = which(substr(dates,1,4)==years[y])
    tempdat = tempdata[,,yearidx]*86400
    tempcount1 = ifelse(tempdat>=thres,1,0)
    tempcount2 = ifelse(tempdat>0,1,0)
    totalvals1 = apply(tempcount1,c(1,2),sum,na.rm=TRUE)
    totalvals2 = apply(tempcount2,c(1,2),sum,na.rm=TRUE)
    rm(tempdat)
    rm(tempcount1)
    rm(tempcount2)
    pryearlydat[,,y]=totalvals1/totalvals2
    rm(totalvals1)
    rm(totalvals2)
    ptm2end = proc.time()
    message("Finished year grab for year ",years[y])
    message("Time: ",ptm2end[3]-ptm2[3]," secs")
  }
  pryearlylist[[i]]=pryearlydat
  rm(pryearlydat)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(filelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}


##############
# 
# Past Climatology

prpast = array(NA,dim=c(length(lon),length(lat),length(filelist)))
yearidx = which(years<=2015 & years>=2006)

for(i in 1:length(filelist)){
  tempdat = apply(pryearlylist[[i]][,,yearidx],c(1,2),mean,na.rm=TRUE)
  prpast[,,i]= ifelse(domainmask==1,tempdat,NA)
}

##############
# 
# Future Climatology

prfut = array(NA,dim=c(length(lon),length(lat),length(filelist)))
yearidx = which(years>=2041 & years<=2070)

for(i in 1:length(filelist)){
  tempdat = apply(pryearlylist[[i]][,,yearidx],c(1,2),mean,na.rm=TRUE)
  prfut[,,i]=ifelse(domainmask==1,tempdat,NA)
}

####################
#
#  Plot Absolute Change by file

filebreakdown = do.call(rbind,strsplit(filelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),9)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")

############################################################
####################################################
# Full 3^5 domain

###########################
# All Individual Plots

diffpr = prfut-prpast
diffrange = range(diffpr,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfut,na.rm=TRUE),range(prpast,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.1,0.1)
breaksdiff = c(seq(-0.1,0.1,by=0.01))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.35)
breaksobs = c(seq(0,0.35,by=0.025))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/IndividualMembers_prratio.pdf",onefile=TRUE,width=8,height=8)

for(i in 1:length(filelist)){
  
  GCM = filebreakdown3$GCM[i]
  scen = filebreakdown3$scen[i]
  obs = filebreakdown3$obs[i]
  DS = filebreakdown3$DS[i]
  
  testsfc1 = list(x=lon,y=lat,z=diffpr[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  #surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpast[,,i])
  #testsfc2 = list(x=lon,y=lat,z=((prpast[,,i]-mean(prpast[,,i],na.rm=TRUE))/mean(prpast[,,i],na.rm=TRUE))*100)
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  #surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfut[,,i])
  #testsfc3 = list(x=lon,y=lat,z=((prfut[,,i]-mean(prfut[,,i],na.rm=TRUE))/mean(prfut[,,i],na.rm=TRUE))*100)
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  #surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
}

dev.off()

###########################
# Group 1 - by emissions scenario

scens = unique(filebreakdown3$scen)
prfutg1 = prpastg1 = diffprg1 = array(NA,dim=c(length(lon),length(lat),length(scens)))

for(s in 1:length(scens)){
  scenidx = which(filebreakdown3$scen==scens[s])
  prfutg1[,,s]= apply(prfut[,,scenidx],c(1,2),mean,na.rm=TRUE)
  prpastg1[,,s]= apply(prpast[,,scenidx],c(1,2),mean,na.rm=TRUE)
  diffprg1[,,s]= apply(diffpr[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

diffrange = range(diffprg1,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfutg1,na.rm=TRUE),range(prpastg1,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.04,0.04)
breaksdiff = c(seq(-0.04,0.04,by=0.005))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.2)
breaksobs = c(seq(0,0.2,by=0.05))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group1members_prratio.pdf",width=8,height=8,onefile=TRUE)

for(s in 1:length(scens)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg1[,,s])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scens[s],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg1[,,s])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scens[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg1[,,s])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scens[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

###########################
# Group 2 - by emissions scenario and GCM

filebreakdown3$scenGCM = paste(filebreakdown3$GCM,filebreakdown3$scen,sep="-")
scenGCMs = unique(filebreakdown3$scenGCM)
prfutg2 = prpastg2 = diffprg2 = array(NA,dim=c(length(lon),length(lat),length(scenGCMs)))

for(s in 1:length(scenGCMs)){
  scenidx = which(filebreakdown3$scenGCM==scenGCMs[s])
  prfutg2[,,s]= apply(prfut[,,scenidx],c(1,2),mean,na.rm=TRUE)
  prpastg2[,,s]= apply(prpast[,,scenidx],c(1,2),mean,na.rm=TRUE)
  diffprg2[,,s]= apply(diffpr[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

diffrange = range(diffprg2,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfutg2,na.rm=TRUE),range(prpastg2,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.07,0.07)
breaksdiff = c(seq(-0.07,0.07,by=0.005))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.25)
breaksobs = c(seq(0,0.25,by=0.025))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group2Members_prratio.pdf",onefile=TRUE,width=8,height=8)

for(s in 1:length(scenGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg2[,,s])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nGCM and Scen: ",scenGCMs[s],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg2[,,s])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nGCM and Scen: ",scenGCMs[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg2[,,s])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nGCM and Scen: ",scenGCMs[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()
###########################
# Group 3 - by emissions scenario and observed dataset

filebreakdown3$obsGCM = paste(filebreakdown3$obs,filebreakdown3$scen,sep="-")
obsGCMs = unique(filebreakdown3$obsGCM)
prfutg3 = prpastg3 = diffprg3 = array(NA,dim=c(length(lon),length(lat),length(obsGCMs)))

for(s in 1:length(obsGCMs)){
  scenidx = which(filebreakdown3$obsGCM==obsGCMs[s])
  prfutg3[,,s]= apply(prfut[,,scenidx],c(1,2),mean,na.rm=TRUE)
  prpastg3[,,s]= apply(prpast[,,scenidx],c(1,2),mean,na.rm=TRUE)
  diffprg3[,,s]= apply(diffpr[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

diffrange = range(diffprg3,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfutg3,na.rm=TRUE),range(prpastg3,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.06,0.06)
breaksdiff = c(seq(-0.06,0.06,by=0.005))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.35)
breaksobs = c(seq(0,0.35,by=0.05))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group3Members_prratio.pdf",onefile=TRUE,width=8,height=8)

for(s in 1:length(obsGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg3[,,s])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nObs and Scen: ",obsGCMs[s],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg2[,,s])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nObs and Scen: ",obsGCMs[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg2[,,s])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nObs and Scen: ",obsGCMs[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

############################################################
####################################################
# Mask to Oklahoma

OKoutline = data.frame(map("state",regions="oklahoma",plot=FALSE)[c("x","y")])
rows = 1:length(lon)
cols = 1:length(lat)
LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(rows,each=length(lat))
C = rep(cols,length(lon))
modgrid = data.frame(R,C,LON,LAT)
library(sp)
modinOK = point.in.polygon(modgrid$LON,modgrid$LAT,OKoutline$x,OKoutline$y)
modgrid = cbind(modgrid,modinOK)
mask = matrix(NA,nrow=length(lon),ncol=length(lat))

for(r in 1:length(lon)){
  for(c in 1:length(lat)){
    mask[r,c] = modgrid$modinOK[which(modgrid$R==r & modgrid$C==c)]
  }
}

lonOK = modgrid$LON[which(modinOK==1)]
latOK = modgrid$LAT[which(modinOK==1)]
lonOKrange = range(lonOK)
lonOKrange[1] = floor(lonOKrange[1])-0.25
lonOKrange[2] = ceiling(lonOKrange[2])

latOKrange = range(latOK)
latOKrange[1] = floor(latOKrange[1])
latOKrange[2] = ceiling(latOKrange[2])+0.25

diffprOK = diffpr
diffprg1OK = diffprg1
diffprg2OK = diffprg2
diffprg3OK = diffprg3
prfutOK = prfut
prfutg1OK = prfutg1
prfutg2OK = prfutg2
prfutg3OK = prfutg3
prpastOK = prpast
prpastg1OK = prpastg1
prpastg2OK = prpastg2
prpastg3OK = prpastg3

for(i in 1:dim(diffpr)[3]){
  diffprOK[,,i]=ifelse(mask==1,diffpr[,,i],NA)
  prfutOK[,,i]=ifelse(mask==1,prfut[,,i],NA)
  prpastOK[,,i]=ifelse(mask==1,prpast[,,i],NA)
}
for(i in 1:dim(diffprg1)[3]){
  diffprg1OK[,,i]=ifelse(mask==1,diffprg1[,,i],NA)
  prfutg1OK[,,i]=ifelse(mask==1,prfutg1[,,i],NA)
  prpastg1OK[,,i]=ifelse(mask==1,prpastg1[,,i],NA)
}
for(i in 1:dim(diffprg2)[3]){
  diffprg2OK[,,i]=ifelse(mask==1,diffprg2[,,i],NA)
  prfutg2OK[,,i]=ifelse(mask==1,prfutg2[,,i],NA)
  prpastg2OK[,,i]=ifelse(mask==1,prpastg2[,,i],NA)
}
for(i in 1:dim(diffprg3)[3]){
  diffprg3OK[,,i]=ifelse(mask==1,diffprg3[,,i],NA)
  prfutg3OK[,,i]=ifelse(mask==1,prfutg3[,,i],NA)
  prpastg3OK[,,i]=ifelse(mask==1,prpastg3[,,i],NA)
}

#############################
####################
# Individual Plots

diffrange = range(diffprOK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfutOK,na.rm=TRUE),range(prpastOK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.08,0.08)
breaksdiff = c(seq(-0.08,0.08,by=0.005))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.30)
breaksobs = c(seq(0,0.30,by=0.03))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/IndividualMembers_prratio_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(filelist)){
  
  GCM = filebreakdown3$GCM[i]
  scen = filebreakdown3$scen[i]
  obs = filebreakdown3$obs[i]
  DS = filebreakdown3$DS[i]
  
  testsfc1 = list(x=lon,y=lat,z=diffprOK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastOK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutOK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
}
dev.off()

#################################
# Group 1

diffrange = range(diffprg1OK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfutg1OK,na.rm=TRUE),range(prpastg1OK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.025,0.025)
breaksdiff = c(seq(-0.025,0.025,by=0.005))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.18)
breaksobs = c(seq(0,0.18,by=0.03))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group1Members_prratio_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(scens)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg1OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scens[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg1OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scens[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg1OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scens[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
}
dev.off()

#################################
# Group 2

diffrange = range(diffprg2OK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfutg2OK,na.rm=TRUE),range(prpastg2OK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.06,0.06)
breaksdiff = c(seq(-0.06,0.06,by=0.01))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.20)
breaksobs = c(seq(0,0.20,by=0.02))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group2Members_prratio_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(scenGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg2OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg2OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg2OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
}
dev.off()

#################################
# Group 3

diffrange = range(diffprg3OK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(prfutg3OK,na.rm=TRUE),range(prpastg3OK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-0.035,0.035)
breaksdiff = c(seq(-0.035,0.035,by=0.005))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,0.27)
breaksobs = c(seq(0,0.27,by=0.3))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group3Members_prratio_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(obsGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg3OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  #testsfc2 = list(x=lon,y=lat,z=prpastg3OK[,,i])
  #surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  #map("state",add=TRUE)
  
  #testsfc3 = list(x=lon,y=lat,z=prfutg3OK[,,i])
  #surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  #map("state",add=TRUE)
  
}
dev.off()

######################
# Extra Analyses

###
# important object names up to this point
# diffpr - projected change full domain all members
# diffprg1 - projected change full domain ensemble mean by scenario
# diffprg2 - projected change full domain ensemble mean by GCM and scenario
# diffprg3 - projected change full domain ensemble mean by Obs and scenario
# 
# diffprOK - projected change Oklahoma all members
# diffprg1OK - projected change Oklahoma ensemble mean by scenario
# diffprg2OK - projected change Oklahoma ensemble mean by GCM and scenario
# diffprg3OK - projected change Oklahoma ensemble mean by Obs and scenario
# 
# prpastXXXX, prfutXXXX - past and future climatology for number of days > 95F
#
# filebreakdown3 - the file metadata in order
# scens - unique scenarios for group1
# scenGCMs - unique GCMs and scenarios for group2
# obsGCMs - unique obs and scenarios for group3
#
###

which(is.na(diffprOK[,,1])==FALSE,arr.ind=TRUE)

for(i in 1:length(filelist)){
  
  if(i==1){
    temp1 = diffprOK[,,i]
    denstemp=density(temp1[which(is.na(temp1)==FALSE)])
    plot(denstemp$y~denstemp$x,type="l",col="red")
  } else {
    temp1 = diffprOK[,,i]
    denstemp=density(temp1[which(is.na(temp1)==FALSE)])
    lines(denstemp$y~denstemp$x,col="green")
  }
  
}




