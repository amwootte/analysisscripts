##############
#
# Sanderson weight calculator
# plus consistency

type="ann"
plotname = "mv"
statemask = "new mexico"

load(paste("/home/woot0002/DS_ind/RMSEmats_combohc_LOCA_",plotname,"_",type,"_",statemask,".Rdata",sep=""))
load(paste("/home/woot0002/DS_ind/RMSEmats_combohc_GCM_",plotname,"_",type,"_",statemask,".Rdata",sep=""))

#load("/home/woot0002/RMSEfiles/tasmin_RMSEmats.Rdata")

numobs = 1
numtotal = nrow(NRMSEhmatout_GCM)
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
S_top_GCMh = -(NRMSEhmatout_GCM/Du)^2
S_GCMh = exp(S_top_GCMh)
S_top_LOCAh = -(NRMSEhmatout_LOCA/Du)^2
S_LOCAh = exp(S_top_LOCAh)

S_top_GCMc = -(NRMSEcmatout_GCM/Du)^2
S_GCMc = exp(S_top_GCMc)
S_top_LOCAc = -(NRMSEcmatout_LOCA/Du)^2
S_LOCAc = exp(S_top_LOCAc)

Ru_GCMh = Ru_LOCAh = Ru_GCMc = Ru_LOCAc = c()
for(i in 1:nummodels){
  Ru_GCMh[i]=1+sum(S_GCMh[1:nummodels,i],na.rm=TRUE)
  Ru_LOCAh[i]=1+sum(S_LOCAh[1:nummodels,i],na.rm=TRUE)
  Ru_GCMc[i]=1+sum(S_GCMc[1:nummodels,i],na.rm=TRUE)
  Ru_LOCAc[i]=1+sum(S_LOCAc[1:nummodels,i],na.rm=TRUE)
}

Wu_GCMh = Ru_GCMh^(-1)
Wu_LOCAh = Ru_LOCAh^(-1)
Wu_GCMc = Ru_GCMc^(-1)
Wu_LOCAc = Ru_LOCAc^(-1)

##
# quality weight calculation

Wq_GCMh = t_dist_GCMh = c()
Wq_LOCAh = t_dist_LOCAh = c()
Wq_GCMc = t_dist_GCMc = c()
Wq_LOCAc = t_dist_LOCAc = c()

for(i in 1:nummodels){
  useidx = obsidx
  
  t_dist_GCMh[i] = NRMSEhmatout_GCM[useidx,i]
  t_dist_LOCAh[i] = NRMSEhmatout_LOCA[useidx,i]
  
  Q_top_GCMh = -(NRMSEhmatout_GCM[useidx,i]/Dq)^2
  Q_top_LOCAh = -(NRMSEhmatout_LOCA[useidx,i]/Dq)^2
  
  Wq_LOCAh[i] = exp(Q_top_LOCAh)
  Wq_GCMh[i] = exp(Q_top_GCMh)
  
  t_dist_GCMc[i] = NRMSEcmatout_GCM[useidx,i]
  t_dist_LOCAc[i] = NRMSEcmatout_LOCA[useidx,i]
  
  Q_top_GCMc = -(NRMSEcmatout_GCM[useidx,i]/Dq)^2
  Q_top_LOCAc = -(NRMSEcmatout_LOCA[useidx,i]/Dq)^2
  
  Wq_LOCAc[i] = exp(Q_top_LOCAc)
  Wq_GCMc[i] = exp(Q_top_GCMc)
  
}

########
# Calculate final weight

W_LOCAh = Wq_LOCAh*Wu_LOCAh
WA_LOCAh = W_LOCAh/sum(W_LOCAh) # normalize so weights add up to 1.

W_GCMh = Wq_GCMh*Wu_GCMh
WA_GCMh = W_GCMh/sum(W_GCMh) # normalize so weights add up to 1.

W_LOCAc = Wq_LOCAc*Wu_LOCAc
WA_LOCAc = W_LOCAc/sum(W_LOCAc) # normalize so weights add up to 1.

W_GCMc = Wq_GCMc*Wu_GCMc
WA_GCMc = W_GCMc/sum(W_GCMc) # normalize so weights add up to 1.

#########
# 

LOCAhdat = LOCAhdat[1:nummodels,]
LOCAhdat$product = paste(LOCAhdat$DS,LOCAhdat$training,sep="_")
LOCAhdat$Wuh = Wu_LOCAh
LOCAhdat$Wqh = Wq_LOCAh
LOCAhdat$Wuc = Wu_LOCAc
LOCAhdat$Wqc = Wq_LOCAc
W_LOCAh = Wq_LOCAh*Wu_LOCAh
LOCAhdat$Wh = W_LOCAh/sum(W_LOCAh) # normalize so weights add up to 1.
W_LOCAc = Wq_LOCAc*Wu_LOCAc
LOCAhdat$Wc = W_LOCAc/sum(W_LOCAc) # normalize so weights add up to 1.

GCMhdat = GCMhdat[1:nummodels,]
GCMhdat$product = paste(GCMhdat$DS,GCMhdat$training,sep="_")
GCMhdat$Wuh = Wu_GCMh
GCMhdat$Wqh = Wq_GCMh
GCMhdat$Wuc = Wu_GCMc
GCMhdat$Wqc = Wq_GCMc
W_GCMh = Wq_GCMh*Wu_GCMh
GCMhdat$Wh = W_GCMh/sum(W_GCMh) # normalize so weights add up to 1.
W_GCMc = Wq_GCMc*Wu_GCMc
GCMhdat$Wc = W_GCMc/sum(W_GCMc) # normalize so weights add up to 1.

GCMhdat[order(GCMhdat$Wh)[1:3],]
GCMhdat[order(GCMhdat$Wc)[1:3],]

LOCAhdat[order(LOCAhdat$Wh)[1:3],]
LOCAhdat[order(LOCAhdat$Wc)[1:3],]

save(list=c("GCMhdat","LOCAhdat"),file=paste("/home/woot0002/DS_ind/Sanderson_EnsembleWeights_",plotname,"_",type,"_",statemask,".Rdata",sep=""))


chxi1 <- chull(x=LOCAhdat$Wqh,y=LOCAhdat$Wuh)
chx1 <- rbind(LOCAhdat[chxi1,], LOCAhdat[chxi1[1], ])
chulllist_LOCAh = chx1

chxi1 <- chull(x=GCMhdat$Wqh,y=GCMhdat$Wuh)
chx1 <- rbind(GCMhdat[chxi1,], GCMhdat[chxi1[1], ])
chulllist_GCMh = chx1

chxi1 <- chull(x=LOCAhdat$Wqc,y=LOCAhdat$Wuc)
chx1 <- rbind(LOCAhdat[chxi1,], LOCAhdat[chxi1[1], ])
chulllist_LOCAc = chx1

chxi1 <- chull(x=GCMhdat$Wqc,y=GCMhdat$Wuc)
chx1 <- rbind(GCMhdat[chxi1,], GCMhdat[chxi1[1], ])
chulllist_GCMc = chx1

plot(Wqh~Wuh,data=chulllist_GCMh,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),col="black",ylab="Skill Weight",xlab="Independence Weight")
lines(Wqh~Wuh,lwd=2,data=chulllist_LOCAh,col="blue") 
abline(coef=c(0,1),lty=2)
legend("bottomright",legend=c("CMIP5","LOCA"),col=c("black","blue"),lwd=2)



    plot(Wqh~Wuh,data=chulllist_GCMh,main=statemask,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),col="black",ylab="Skill Weight",xlab="Independence Weight")
    lines(Wqh~Wuh,lwd=2,data=chulllist_LOCAh,col="blue") 
    
    lines(Wqc~Wuc,lwd=2,data=chulllist_GCMc,col="black",lty=2) 
    lines(Wqc~Wuc,lwd=2,data=chulllist_LOCAc,col="blue",lty=2) 
    
    abline(coef=c(0,1),lty=2)
legend("top",horiz=FALSE,cex=0.75,legend=c("CMIP5 historical","LOCA historical","CMIP5 change","LOCA change"),col=c("black","blue","black","blue"),lwd=2,lty=c(1,1,2,2))


product = rep(rep(c("CMIP5","LOCA"),each=length(W_LOCAh)),2)

Wu= c(Wu_GCMh,Wu_LOCAh)
Wq= c(Wq_GCMh,Wq_LOCAh)
W = c(W_GCMh,W_LOCAh)
simdat1 = data.frame(product,Wu,Wq,W)

Wu= c(Wu_GCMc,Wu_LOCAc)
Wq= c(Wq_GCMc,Wq_LOCAc)
W = c(W_GCMc,W_LOCAc)
simdat2 = data.frame(product,Wu,Wq,W)

simdat = rbind(simdat1,simdat2)
simdat$group=rep(c("historical","change"),each=nrow(simdat1))


library(ggplot2)
dg<-qplot(Wu,Wq,colour=W,shape=product,data=simdat,size=I(3)) 
dg + scale_shape_manual(values=1:2, aesthetics="fill") + scale_colour_gradient2(low="red", high="blue",mid = "gray96", 
                            midpoint = 0.4,limits=c(0,0.8),breaks=seq(0,0.8,by=0.1),guide="legend") +
  xlim(0,1)+ylim(0,1) + ylab("Skill Weight") + xlab("Independence Weight") + theme_bw() +facet_grid(rows=vars(group))


###

ggplot(simdat, aes(x=Wu, y=Wq,colour=group,shape=product)) + 
  geom_point(size=2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ggtitle("Skill vs. Independence weights")+xlab("Independence Weight")+ylab("Skill Weight")+ylim(0,1)+xlim(0,1)


aggregate(simdat[,2:4],by=list(productgroup=paste(simdat$product,simdat$group,sep="_")),mean,na.rm=TRUE)

