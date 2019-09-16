##############
#
# Distance plotter

load("/home/woot0002/RMSEmats_combo.Rdata")
numobs = 5
numtotal = nrow(normRMSEmat)
nummodels = numtotal-numobs
obsidx = (nummodels+1):numtotal
obsnames = simnames[obsidx]

metadat$modavg = NA
for(i in 1:nummodels){
  metadat$modavg[i] = mean(normRMSEmat[i,1:nummodels],na.rm=TRUE)
}

metadat$travg = metadat$outavg = NA
for(i in 1:nummodels){
  
  if(metadat$training[i]=="Daymet"){
    tridx = obsidx[1]
  }
  if(metadat$training[i]=="Livneh"){
    tridx = obsidx[2]
  }
  if(metadat$training[i]=="PRISM"){
    tridx = obsidx[3]
  }
  if(metadat$training[i]=="METDATA"){
    tridx = obsidx[4]
  }
  if(metadat$training[i]=="MAURER"){
    tridx = obsidx[5]
  }
  
  if(metadat$training[i]=="Daymet"){
    outidx = obsidx[2:5]
  }
  if(metadat$training[i]=="Livneh"){
    outidx = obsidx[c(1,3:5)]
  }
  if(metadat$training[i]=="PRISM"){
    outidx = obsidx[c(1:2,4:5)]
  }
  if(metadat$training[i]=="METDATA"){
    outidx = obsidx[c(1:3,5)]
  }
  if(metadat$training[i]=="MAURER"){
    outidx = obsidx[1:4]
  }
  
  metadat$travg[i] = normRMSEmat[i,tridx]
  metadat$outavg[i] = mean(normRMSEmat[i,outidx],na.rm=TRUE)
}

#############

plot(modavg~travg,pch=19,data=metadat)

metadat$product = paste(metadat$DS,metadat$training,sep="_")
products=unique(metadat$product)[1:7]
cols = c("red","blue","green","purple","pink","orange","black")

for(p in 1:length(products)){
  if(p==1){
    plot(modavg~travg,data=metadat[which(metadat$product==products[p]),],xlim=c(0,round(max(metadat[,5:7],na.rm=TRUE),1)),ylim=c(0,round(max(metadat[,5:7],na.rm=TRUE),1)),col=cols[p])
  } else {
    points(modavg~travg,data=metadat[which(metadat$product==products[p]),],col=cols[p]) 
  }
}
abline(coef=c(0,1),lty=2)
legend("bottomright",legend=products,col=cols,pch=1)

chulllist = list()
for(p in 1:length(products)){
  tmp = subset(metadat,product==products[p])
  chxi1 <- chull(x=tmp$travg,y=tmp$modavg)
  chx1 <- rbind(tmp[chxi1,], tmp[chxi1[1], ])
  chulllist[[p]] = chx1
  
  if(p==1){
    plot(modavg~travg,data=chx1,type="l",lwd=2,xlim=c(0.5,round(max(metadat[,5:7],na.rm=TRUE),1)),ylim=c(0.5,round(max(metadat[,5:7],na.rm=TRUE),1)),col=cols[p])
  } else {
    points(modavg~travg,data=chx1,type="l",lwd=2,col=cols[p]) 
  }
  
}
abline(coef=c(0,1),lty=2)
legend("bottomright",legend=products,col=cols,pch=1)


chulllist = list()
for(p in 1:length(products)){
  tmp = subset(metadat,product==products[p])
  chxi1 <- chull(x=tmp$outavg,y=tmp$modavg)
  chx1 <- rbind(tmp[chxi1,], tmp[chxi1[1], ])
  chulllist[[p]] = chx1
  
  if(p==1){
    plot(modavg~outavg,data=chx1,type="l",lwd=2,xlim=c(0,round(max(metadat[,5:7],na.rm=TRUE),1)),ylim=c(0,round(max(metadat[,5:7],na.rm=TRUE),1)),col=cols[p])
  } else {
    points(modavg~outavg,data=chx1,type="l",lwd=2,col=cols[p]) 
  }
  
}
abline(coef=c(0,1),lty=2)
legend("bottomright",legend=products,col=cols,pch=1)






