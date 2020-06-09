##############
#
# Sanderson weight calculator
# plus consistency

load("/home/woot0002/RMSEmats_combo_LOCA.Rdata")
load("/home/woot0002/RMSEmats_combo_GCM.Rdata")

#load("/home/woot0002/RMSEfiles/tasmin_RMSEmats.Rdata")

numobs = 1
numtotal = nrow(NRMSEmatout_GCM)
nummodels = numtotal-numobs
obsidx = (nummodels+1):numtotal
obsnames = simnames_GCM[obsidx]

##
# set radii
#Du = 0.6 # original 0.48 - ???? does this need to be 0.48 * distance between best performing model and obs.
#Dq = 0.7
#Dc = 0.7

Du = 0.48 # original 0.48 - ???? does this need to be 0.48 * distance between best performing model and obs.
Dq = 0.8

##
# similarity weight calculation
S_top_GCM = -(NRMSEmatout_GCM/Du)^2
S_GCM = exp(S_top_GCM)

S_top_LOCA = -(NRMSEmatout_LOCA/Du)^2
S_LOCA = exp(S_top_LOCA)

Ru_GCM = Ru_LOCA = c()
for(i in 1:nummodels){
  Ru_GCM[i]=1+sum(S_GCM[1:nummodels,i],na.rm=TRUE)
  Ru_LOCA[i]=1+sum(S_LOCA[1:nummodels,i],na.rm=TRUE)
}

Wu_GCM = Ru_GCM^(-1)
Wu_LOCA = Ru_LOCA^(-1)

##
# quality weight calculation

Wq_GCM = t_dist_GCM = c()
Wq_LOCA = t_dist_LOCA = c()
for(i in 1:nummodels){
  useidx = obsidx
  
  t_dist_GCM[i] = NRMSEmatout_GCM[useidx,i]
  t_dist_LOCA[i] = NRMSEmatout_LOCA[useidx,i]
  
  Q_top_GCM = -(NRMSEmatout_GCM[useidx,i]/Dq)^2
  Q_top_LOCA = -(NRMSEmatout_LOCA[useidx,i]/Dq)^2
  
  Wq_LOCA[i] = exp(Q_top_LOCA)
  Wq_GCM[i] = exp(Q_top_GCM)
}

########
# Calculate final weight

W_LOCA = Wq_LOCA*Wu_LOCA
WA_LOCA = W_LOCA/sum(W_LOCA) # normalize so weights add up to 1.

W_GCM = Wq_GCM*Wu_GCM
WA_GCM = W_GCM/sum(W_GCM) # normalize so weights add up to 1.

#########
# 

LOCAdat = LOCAdat[1:nummodels,]
LOCAdat$product = paste(LOCAdat$DS,LOCAdat$training,sep="_")
LOCAdat$Wu = Wu_LOCA
LOCAdat$Wq = Wq_LOCA

GCMdat = GCMdat[1:nummodels,]
GCMdat$product = paste(GCMdat$DS,GCMdat$training,sep="_")
GCMdat$Wu = Wu_GCM
GCMdat$Wq = Wq_GCM

chxi1 <- chull(x=LOCAdat$Wq,y=LOCAdat$Wu)
chx1 <- rbind(LOCAdat[chxi1,], LOCAdat[chxi1[1], ])
chulllist_LOCA = chx1

chxi1 <- chull(x=GCMdat$Wq,y=GCMdat$Wu)
chx1 <- rbind(GCMdat[chxi1,], GCMdat[chxi1[1], ])
chulllist_GCM = chx1


    plot(Wq~Wu,data=chulllist_GCM,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),col="black",ylab="Skill Weight",xlab="Independence Weight")
    lines(Wq~Wu,lwd=2,data=chulllist_LOCA,col="blue") 
    abline(coef=c(0,1),lty=2)
legend("topright",legend=c("GCMs","LOCA"),col=c("black","blue"),lwd=2)

product = rep(c("CMIP5","LOCA"),each=length(W_LOCA))
Wu= c(Wu_GCM,Wu_LOCA)
Wq= c(Wq_GCM,Wq_LOCA)
W = c(W_GCM,W_LOCA)
simdat = data.frame(product,Wu,Wq,W)
library(ggplot2)
dg<-qplot(Wu,Wq,colour=W,shape=product,data=simdat,size=I(3)) 
dg + scale_shape_manual(values=1:2, aesthetics="fill") + scale_colour_gradient2(low="red", high="blue",mid = "gray96", 
                            midpoint = 0.1,limits=c(0,0.2),breaks=seq(0,0.2,by=0.02),guide="legend") +
  xlim(0,1)+ylim(0,1) + ylab("Skill Weight") + xlab("Independence Weight") + theme_bw()


###


