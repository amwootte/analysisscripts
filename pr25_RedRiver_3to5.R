###################
# Load RR mask

test = nc_open("/home/woot0002/RRdomainmask.nc")
mask = ncvar_get(test,"RRmask")
nc_close(test)

######################
# mask to Red River

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

zlimdiff = c(-7,7)
breaksdiff = c(seq(-7,7,by=0.5))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,30)
breaksobs = c(seq(0,30,by=3))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/IndividualMembers_pr25_RR.pdf",onefile=TRUE,width=8,height=8)

for(i in 1:length(filelist)){
  
  GCM = filebreakdown3$GCM[i]
  scen = filebreakdown3$scen[i]
  obs = filebreakdown3$obs[i]
  DS = filebreakdown3$DS[i]
  
  testsfc1 = list(x=lon,y=lat,z=diffprOK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastOK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutOK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
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

zlimdiff = c(-3,3)
breaksdiff = c(seq(-3,3,by=0.25))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,30)
breaksobs = c(seq(0,30,by=3))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group1Members_pr25_RR.pdf",onefile=TRUE,width=8,height=8)

for(i in 1:length(scens)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg1OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scens[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg1OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nScen: ",scens[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg1OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nScen: ",scens[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
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

zlimdiff = c(-6,6)
breaksdiff = c(seq(-6,6,by=0.5))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,30)
breaksobs = c(seq(0,30,by=3))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group2Members_pr25_RR.pdf",onefile=TRUE,width=8,height=8)

for(i in 1:length(scenGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg2OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg2OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg2OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nGCM and Scen: ",scenGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
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

zlimdiff = c(-3,3)
breaksdiff = c(seq(-3,3,by=0.25))
colorbardiff = colorRampPalette(c("chocolate4","grey96","darkgreen"))(length(breaksdiff)-1)

zlimobs = c(0,30)
breaksobs = c(seq(0,30,by=3))
colorbarobs = colorRampPalette(c("grey96","lightgreen","darkgreen"))(length(breaksobs)-1)

pdf("/home/woot0002/Group3Members_pr25_RR.pdf",onefile=TRUE,width=8,height=8)

for(i in 1:length(obsGCMs)){
  
  testsfc1 = list(x=lon,y=lat,z=diffprg3OK[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimdiff,col=colorbardiff,breaks=breaksdiff,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc2 = list(x=lon,y=lat,z=prpastg3OK[,,i])
  surface(testsfc2,type="I",main=paste("Historical Climate (2006-2015)\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
  testsfc3 = list(x=lon,y=lat,z=prfutg3OK[,,i])
  surface(testsfc3,type="I",main=paste("Projection (2041-2070)\nObs and Scen: ",obsGCMs[i],sep=""),zlim=zlimobs,col=colorbarobs,breaks=breaksobs,xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
  
}
dev.off()

#######################
#######################



