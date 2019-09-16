##############
#
# Sanderson weight calculator
# plus consistency

load("/home/woot0002/RMSEmats_combo.Rdata")

#load("/home/woot0002/RMSEfiles/tasmin_RMSEmats.Rdata")

numobs = 5
numtotal = nrow(normRMSEmat)
nummodels = numtotal-numobs
obsidx = (nummodels+1):numtotal
obsnames = simnames[obsidx]

##
# set radii
#Du = 0.6 # original 0.48 - ???? does this need to be 0.48 * distance between best performing model and obs.
#Dq = 0.7
#Dc = 0.7

Du = 0.48 # original 0.48 - ???? does this need to be 0.48 * distance between best performing model and obs.
Dq = 0.8
Dc = 0.8


##
# similarity weight calculation
S_top = -(normRMSEmat/Du)^2
S = exp(S_top)

Ru = c()
for(i in 1:nummodels){
  Ru[i]=1+sum(S[1:nummodels,i],na.rm=TRUE)
}

Wu = Ru^(-1)

##
# quality weight calculation

Wq = t_dist = c()
for(i in 1:nummodels){
  if(metadat$training[i]=="Daymet"){
    useidx = obsidx[1]
  }
  if(metadat$training[i]=="Livneh"){
    useidx = obsidx[2]
  }
  if(metadat$training[i]=="PRISM"){
    useidx = obsidx[3]
  }
  if(metadat$training[i]=="METDATA"){
    useidx = obsidx[4]
  }
  if(metadat$training[i]=="MAURER"){
    useidx = obsidx[5]
  }
  
  t_dist[i] = normRMSEmat[useidx,i]
  
  Q_top = -(normRMSEmat[useidx,i]/Dq)^2
  Wq[i] = exp(Q_top)
}

##
# consitency weight calculation

Wc =c_dist = c()
for(i in 1:nummodels){
  if(metadat$training[i]=="Daymet"){
    useidx = obsidx[2:5]
  }
  if(metadat$training[i]=="Livneh"){
    useidx = obsidx[c(1,3:5)]
  }
  if(metadat$training[i]=="PRISM"){
    useidx = obsidx[c(1:2,4:5)]
  }
  if(metadat$training[i]=="METDATA"){
    useidx = obsidx[c(1:3,5)]
  }
  if(metadat$training[i]=="MAURER"){
    useidx = obsidx[1:4]
  }
  
  c_dist[i]= mean(normRMSEmat[useidx,i],na.rm=TRUE)
  
  C_top = -(mean(normRMSEmat[useidx,i],na.rm=TRUE)/Dc)^2
  Wc[i] = exp(C_top)
}

########
# Calculate final weight

W = Wq*Wc*Wu
WA = W/sum(W) # normalize so weights add up to 1.

W2 = Wc*Wu
WA2 = W2/sum(W2)

#########
# 


simmetdat = metadat[1:nummodels,]
simmetdat$product = paste(simmetdat$DS,simmetdat$training,sep="_")
simmetdat$Wc = Wc
simmetdat$Wu = Wu
simmetdat$Wq = Wq

simmetdat$col = NA

products=unique(simmetdat$product)
cols = c("red","blue","green","purple","pink","orange","black")

for(p in 1:length(products)){
  if(p==1){
    plot(Wc~Wq,data=simmetdat[which(simmetdat$product==products[p]),],xlim=c(0,1),ylim=c(0,1),col=cols[p])
  } else {
    points(Wc~Wq,data=simmetdat[which(simmetdat$product==products[p]),],col=cols[p]) 
  }
}
abline(coef=c(0,1),lty=2)
legend("topright",legend=products,col=cols,pch=1)



chulllist = list()
for(p in 1:length(products)){
  tmp = subset(simmetdat,product==products[p])
  chxi1 <- chull(x=tmp$Wq,y=tmp$Wc)
  chx1 <- rbind(tmp[chxi1,], tmp[chxi1[1], ])
  chulllist[[p]] = chx1
  
  if(p==1){
    plot(Wc~Wq,data=chx1,type="l",lwd=2,xlim=c(0,0.7),ylim=c(0,0.2),col=cols[p])
  } else {
    lines(Wc~Wq,lwd=2,data=chx1,col=cols[p]) 
  }
  
}
abline(coef=c(0,1),lty=2)
legend("bottomright",legend=products,col=cols,lwd=2,cex=0.7)


chulllist = list()
for(p in 1:length(products)){
  tmp = subset(simmetdat,product==products[p])
  chxi1 <- chull(x=tmp$Wq,y=tmp$Wu)
  chx1 <- rbind(tmp[chxi1,], tmp[chxi1[1], ])
  chulllist[[p]] = chx1
  
  if(p==1){
    plot(Wu~Wq,data=chx1,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),col=cols[p])
  } else {
    lines(Wu~Wq,lwd=2,data=chx1,col=cols[p]) 
  }
  
}
abline(coef=c(0,1),lty=2)
legend("bottomright",legend=products,col=cols,lwd=2,cex=0.7)



chulllist = list()
for(p in 1:length(products)){
  tmp = subset(simmetdat,product==products[p])
  chxi1 <- chull(x=tmp$Wc,y=tmp$Wu)
  chx1 <- rbind(tmp[chxi1,], tmp[chxi1[1], ])
  chulllist[[p]] = chx1
  
  if(p==1){
    plot(Wu~Wc,data=chx1,type="l",lwd=2,xlim=c(0,0.3),ylim=c(0,1),col=cols[p])
  } else {
    lines(Wu~Wc,lwd=2,data=chx1,col=cols[p]) 
  }
  
}
abline(coef=c(0,1),lty=2)
legend("bottomright",legend=products,col=cols,lwd=2,cex=0.7)


library(ggplot2) 

dg<-qplot(Wq,Wc,colour=Wu,shape=product,data=simmetdat) 
dg + scale_shape_manual(values=1:7) + scale_colour_gradient2(low="red", high="blue",mid = "gray96", 
                            midpoint = 0.3,limits=c(0,0.6),breaks=seq(0,1,by=0.1),guide="legend") +
  xlim(0,1)+ylim(0,1)


meandat = aggregate(simmetdat[,6:8],by=list(product=simmetdat$product),mean,na.rm=TRUE)
