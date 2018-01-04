####################
#
# 3^5 Analysis with High temperature change
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
#  1. Gather and convert daily tasmax into yearly number of days > 95F
#  2. Calculate climatology for the two time periods
#  3. plot maps for the projected change
#     a. absolute change
#     b. percent change
#     c. past vs. future maps
#
###################################

###########
# 1. Data Gather and conversion
setwd("/data2/3to5/I35/tasmax/EDQM")
filelist = system("ls *.nc",intern=T)

tmax95yearlylist = list()

for(i in 1:length(filelist)){
  ptm = proc.time()
  test = nc_open(filelist[i])
  tempdata = ncvar_get(test,"tasmax")
  
  if(i==1){
    lat = ncvar_get(test,"lat")
    lon = ncvar_get(test,"lon")-360
    times = ncvar_get(test,"time")
    startdate = as.Date(substr(test$dim[[4]]$units,12,21),"%Y-%m-%d")
    timeunits = test$dim[[4]]$units
    dates = startdate+times
    years = 2006:2070
    domainmask = ifelse(is.na(tempdata[,,1])==FALSE,1,0)
  }
  nc_close(test)
  
  #thres = ((95-32)*5/9)+273.15
  tempyearlydat = array(NA,dim=c(length(lon),length(lat),length(years)))
  for(y in 1:length(years)){
    ptm2 = proc.time()
    yearidx = which(substr(dates,1,4)==years[y])
    tempdat = tempdata[,,yearidx]
    #tempcount = ifelse(tempdat>=thres,1,0)
    totalvals = apply(tempdat,c(1,2),mean,na.rm=TRUE)
    tempyearlydat[,,y]=totalvals
    rm(tempdat)
    ptm2end = proc.time()
    message("Finished year grab for year ",years[y])
    message("Time: ",ptm2end[3]-ptm2[3]," secs")
  }
  tmax95yearlylist[[i]]=tempyearlydat
  rm(tempyearlydat)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(filelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}


##############
# 
# Past Climatology

tmax95past = array(NA,dim=c(length(lon),length(lat),length(filelist)))
yearidx = which(years<=2015 & years>=2006)

for(i in 1:length(filelist)){
  tempdat = apply(tmax95yearlylist[[i]][,,yearidx],c(1,2),mean,na.rm=TRUE)
  tmax95past[,,i]= ifelse(domainmask==1,tempdat,NA)
}

##############
# 
# Future Climatology

tmax95fut = array(NA,dim=c(length(lon),length(lat),length(filelist)))
yearidx = which(years>=2041 & years<=2070)

for(i in 1:length(filelist)){
  tempdat = apply(tmax95yearlylist[[i]][,,yearidx],c(1,2),mean,na.rm=TRUE)
  tmax95fut[,,i]=ifelse(domainmask==1,tempdat,NA)
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

difftmax95 = tmax95fut-tmax95past
diffrange = range(difftmax95,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95fut,na.rm=TRUE),range(tmax95past,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(275,310)
breaksobs = c(seq(275,310,by=2.5))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/IndividualMembers_tmax.pdf",onefile=TRUE,width=8,height=8)

for(i in 1:length(filelist)){
  
  GCM = filebreakdown3$GCM[i]
  scen = filebreakdown3$scen[i]
  obs = filebreakdown3$obs[i]
  DS = filebreakdown3$DS[i]
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95past[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95fut[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
}

dev.off()

###########################
# Group 1 - by emissions scenario

scens = unique(filebreakdown3$scen)
tmax95futg1 = tmax95pastg1 = difftmax95g1 = array(NA,dim=c(length(lon),length(lat),length(scens)))

for(s in 1:length(scens)){
  scenidx = which(filebreakdown3$scen==scens[s])
  tmax95futg1[,,s]= apply(tmax95fut[,,scenidx],c(1,2),mean,na.rm=TRUE)
  tmax95pastg1[,,s]= apply(tmax95past[,,scenidx],c(1,2),mean,na.rm=TRUE)
  difftmax95g1[,,s]= apply(difftmax95[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

diffrange = range(difftmax95g1,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95futg1,na.rm=TRUE),range(tmax95pastg1,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(275,310)
breaksobs = c(seq(275,310,by=2.5))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/Group1members_tmax.pdf",width=8,height=8,onefile=TRUE)

for(s in 1:length(scens)){
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95g1[,,s])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scens[s],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95pastg1[,,s])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scens[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95futg1[,,s])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scens[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

###########################
# Group 2 - by emissions scenario and GCM

filebreakdown3$scenGCM = paste(filebreakdown3$GCM,filebreakdown3$scen,sep="-")
scenGCMs = unique(filebreakdown3$scenGCM)
tmax95futg2 = tmax95pastg2 = difftmax95g2 = array(NA,dim=c(length(lon),length(lat),length(scenGCMs)))

for(s in 1:length(scenGCMs)){
  scenidx = which(filebreakdown3$scenGCM==scenGCMs[s])
  tmax95futg2[,,s]= apply(tmax95fut[,,scenidx],c(1,2),mean,na.rm=TRUE)
  tmax95pastg2[,,s]= apply(tmax95past[,,scenidx],c(1,2),mean,na.rm=TRUE)
  difftmax95g2[,,s]= apply(difftmax95[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

diffrange = range(difftmax95g2,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95futg2,na.rm=TRUE),range(tmax95pastg2,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(275,310)
breaksobs = c(seq(275,310,by=2.5))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/Group2Members_tmax.pdf",onefile=TRUE,width=8,height=8)

for(s in 1:length(scenGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95g2[,,s])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nGCM and Scen: ",scenGCMs[s],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95pastg2[,,s])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nGCM and Scen: ",scenGCMs[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95futg2[,,s])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nGCM and Scen: ",scenGCMs[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

###########################
# Group 3 - by emissions scenario and observed dataset

filebreakdown3$obsGCM = paste(filebreakdown3$obs,filebreakdown3$scen,sep="-")
obsGCMs = unique(filebreakdown3$obsGCM)
tmax95futg3 = tmax95pastg3 = difftmax95g3 = array(NA,dim=c(length(lon),length(lat),length(obsGCMs)))

for(s in 1:length(obsGCMs)){
  scenidx = which(filebreakdown3$obsGCM==obsGCMs[s])
  tmax95futg3[,,s]= apply(tmax95fut[,,scenidx],c(1,2),mean,na.rm=TRUE)
  tmax95pastg3[,,s]= apply(tmax95past[,,scenidx],c(1,2),mean,na.rm=TRUE)
  difftmax95g3[,,s]= apply(difftmax95[,,scenidx],c(1,2),mean,na.rm=TRUE)
}

diffrange = range(difftmax95g3,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95futg3,na.rm=TRUE),range(tmax95pastg3,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(275,310)
breaksobs = c(seq(275,310,by=2.5))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/Group3Members_tmax.pdf",onefile=TRUE,width=8,height=8)

for(s in 1:length(obsGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95g3[,,s])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nObs and Scen: ",obsGCMs[s],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95pastg2[,,s])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nObs and Scen: ",obsGCMs[s],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95futg2[,,s])
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

difftmax95OK = difftmax95
difftmax95g1OK = difftmax95g1
difftmax95g2OK = difftmax95g2
difftmax95g3OK = difftmax95g3
tmax95futOK = tmax95fut
tmax95futg1OK = tmax95futg1
tmax95futg2OK = tmax95futg2
tmax95futg3OK = tmax95futg3
tmax95pastOK = tmax95past
tmax95pastg1OK = tmax95pastg1
tmax95pastg2OK = tmax95pastg2
tmax95pastg3OK = tmax95pastg3

for(i in 1:dim(difftmax95)[3]){
  difftmax95OK[,,i]=ifelse(mask==1,difftmax95[,,i],NA)
  tmax95futOK[,,i]=ifelse(mask==1,tmax95fut[,,i],NA)
  tmax95pastOK[,,i]=ifelse(mask==1,tmax95past[,,i],NA)
}
for(i in 1:dim(difftmax95g1)[3]){
  difftmax95g1OK[,,i]=ifelse(mask==1,difftmax95g1[,,i],NA)
  tmax95futg1OK[,,i]=ifelse(mask==1,tmax95futg1[,,i],NA)
  tmax95pastg1OK[,,i]=ifelse(mask==1,tmax95pastg1[,,i],NA)
}
for(i in 1:dim(difftmax95g2)[3]){
  difftmax95g2OK[,,i]=ifelse(mask==1,difftmax95g2[,,i],NA)
  tmax95futg2OK[,,i]=ifelse(mask==1,tmax95futg2[,,i],NA)
  tmax95pastg2OK[,,i]=ifelse(mask==1,tmax95pastg2[,,i],NA)
}
for(i in 1:dim(difftmax95g3)[3]){
  difftmax95g3OK[,,i]=ifelse(mask==1,difftmax95g3[,,i],NA)
  tmax95futg3OK[,,i]=ifelse(mask==1,tmax95futg3[,,i],NA)
  tmax95pastg3OK[,,i]=ifelse(mask==1,tmax95pastg3[,,i],NA)
}

#############################
####################
# Individual Plots

diffrange = range(difftmax95OK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95futOK,na.rm=TRUE),range(tmax95pastOK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(290,305)
breaksobs = c(seq(290,305,by=1))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/IndividualMembers_tmax_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(filelist)){
  
  GCM = filebreakdown3$GCM[i]
  scen = filebreakdown3$scen[i]
  obs = filebreakdown3$obs[i]
  DS = filebreakdown3$DS[i]
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95pastOK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95futOK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
}
dev.off()

#################################
# Group 1

diffrange = range(difftmax95g1OK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95futg1OK,na.rm=TRUE),range(tmax95pastg1OK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(290,305)
breaksobs = c(seq(290,305,by=1))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/Group1Members_tmax_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(scens)){
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95g1OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scens[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95pastg1OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scens[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95futg1OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scens[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
}
dev.off()

#################################
# Group 2

diffrange = range(difftmax95g2OK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95futg2OK,na.rm=TRUE),range(tmax95pastg2OK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(290,305)
breaksobs = c(seq(290,305,by=1))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/Group2Members_tmax_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(scenGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95g2OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95pastg2OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95futg2OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
}
dev.off()

#################################
# Group 3

diffrange = range(difftmax95g3OK,na.rm=TRUE)
diffrange[1]=floor(diffrange[1])
diffrange[2]=ceiling(diffrange[2])

if(abs(diffrange[2])>abs(diffrange[1])) diffrange[1]= -diffrange[2]
if(abs(diffrange[1])>abs(diffrange[2])) diffrange[2]= -diffrange[1]

obsrange = range(c(range(tmax95futg3OK,na.rm=TRUE),range(tmax95pastg3OK,na.rm=TRUE)))
obsrange[1]=floor(obsrange[1])
obsrange[2]=ceiling(obsrange[2])

zlimdiff = c(-5,5)
breaksdiff = c(seq(-5,5,by=0.5))
colorbardiff = colorRampPalette(c("blue","grey96","red"))(length(breaksdiff)-1)

zlimobs = c(290,305)
breaksobs = c(seq(290,305,by=1))
colorbarobs = colorRampPalette(c("yellow2","orange","red","darkred"))(length(breaksobs)-1)

pdf("/home/woot0002/Group3Members_tmax_OK.pdf",onefile=TRUE,width=9,height=6)

for(i in 1:length(obsGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=difftmax95g3OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=tmax95pastg3OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=tmax95futg3OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude",xlim=lonOKrange,ylim=latOKrange)
  map("state",add=TRUE)
  
}
dev.off()

######################
# Extra Analyses

###
# important object names up to this point
# difftmax95 - projected change full domain all members
# difftmax95g1 - projected change full domain ensemble mean by scenario
# difftmax95g2 - projected change full domain ensemble mean by GCM and scenario
# difftmax95g3 - projected change full domain ensemble mean by Obs and scenario
# 
# difftmax95OK - projected change Oklahoma all members
# difftmax95g1OK - projected change Oklahoma ensemble mean by scenario
# difftmax95g2OK - projected change Oklahoma ensemble mean by GCM and scenario
# difftmax95g3OK - projected change Oklahoma ensemble mean by Obs and scenario
# 
# tmax95pastXXXX, tmax95futXXXX - past and future climatology for number of days > 95F
#
# filebreakdown3 - the file metadata in order
# scens - unique scenarios for group1
# scenGCMs - unique GCMs and scenarios for group2
# obsGCMs - unique obs and scenarios for group3
#
###

which(is.na(difftmax95OK[,,1])==FALSE,arr.ind=TRUE)

for(i in 1:length(filelist)){
  
  if(i==1){
    temp1 = difftmax95OK[,,i]
    denstemp=density(temp1[which(is.na(temp1)==FALSE)])
    plot(denstemp$y~denstemp$x,type="l",col="red")
  } else {
    temp1 = difftmax95OK[,,i]
    denstemp=density(temp1[which(is.na(temp1)==FALSE)])
    lines(denstemp$y~denstemp$x,col="green")
  }
  
}




