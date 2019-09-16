###################################
#
# SKC MDS and distance matrix analysis
#
###################################

setwd("/home/woot0002")
scenin = "rcp85"

load(paste(scenin,"anomalymatrix.Rdata",sep="_"))

library(smacof)
library(fields)
library(ggplot2)
library(reshape2)
source("SKC_functions.R")

#######################
# Check out the eigen values

svdvals1 = svd(localearlyanom)
U1 = svdvals1$u # this is an m by m matrix, at this untruncated to the t modes
V1 = svdvals1$v # this is the n by m matrix - eventually n by t also untruncated
lambda1 = svdvals1$d # diagonal eigenvalues lambda
lorder = order(lambda1,decreasing=TRUE)
lambda1 = lambda1[lorder]
U1 = U1[,lorder]
V1 = V1[,lorder]

svdvals2 = svd(locallateanom)
U2 = svdvals2$u # this is an m by m matrix, at this untruncated to the t modes
V2 = svdvals2$v # this is the n by m matrix - eventually n by t also untruncated
lambda2 = svdvals2$d # diagonal eigenvalues lambda
lorder = order(lambda2,decreasing=TRUE)
lambda2 = lambda2[lorder]
U2 = U2[,lorder]
V2 = V2[,lorder]

EOFmode=1:11

plot(lambda1~EOFmode,type="b",col="blue",pch=19,xlim=c(1,9),ylim=c(0,1200))
lines(lambda2~EOFmode,type="b",col="red",pch=19)

lambda1sum = sum(lambda1)
lambda2sum = sum(lambda2)

plot(((lambda1/lambda1sum)*100)~EOFmode,type="b",col="blue",pch=19,xlim=c(1,9),ylim=c(0,60))
lines(((lambda2/lambda2sum)*100)~EOFmode,type="b",col="red",pch=19)

lb1 = lb2 = c()

for(i in 1:11){
  if(i==1){
    lb1[i] = lambda1[i]
    lb2[i] = lambda2[i]
  } else {
    lb1[i] = lambda1[i]+lb1[(i-1)]
    lb2[i] = lambda2[i]+lb2[(i-1)]
  }
}

plot(((lb1/lambda1sum)*100)~EOFmode,type="b",col="blue",pch=19,xlim=c(1,9),ylim=c(0,100),lwd=2,main="Cumulative Variance explained by each EOF mode by Ensemble",xlab="EOF Mode")
lines(((lb2/lambda2sum)*100)~EOFmode,type="b",col="red",pch=19,lwd=2)

#######################
# SVD Analysis

U.early = SVDanalysis(localearlyanom,trunc=5)
U.late = SVDanalysis(locallateanom,trunc=5)

#######################
# Distance matrix calcs

del.early = distmatcalc(U.early)
del.late = distmatcalc(U.late)

###############
# Distance matrix visualization

rownames(del.early)=models
colnames(del.early)=models

rownames(del.late)=models
colnames(del.late)=models

get_upper_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

del.earlymat = get_upper_tri(del.early)
del.latemat = get_upper_tri(del.late)

#del.earlymat = ifelse(del.earlymat==0,NA,del.earlymat)

melted_cormat1 <- melt(del.earlymat, na.rm = TRUE)
melted_cormat2 <- melt(del.latemat,na.rm=TRUE)
### Heatmap
ggplot(data = melted_cormat1, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = 1, limit = c(0.5,1.5), breaks=seq(0.5,1.5,by=0.1), space = "Lab", 
                       name="Euclidean\nDistance") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

###
ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = 1, limit = c(0.5,1.5), breaks=seq(0.5,1.5,by=0.1), space = "Lab", 
                       name="Euclidean\nDistance") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#######################
# Model weights from model distances - using the domain average calculation only

DuE = quantile(apply(ifelse(del.early==0,NA,del.early),1,min,na.rm=TRUE),probs=0.5)
DuL = quantile(apply(ifelse(del.late==0,NA,del.late),1,min,na.rm=TRUE),probs=0.5)

S = del.early
S2 = del.late

for(i in 1:nrow(S)){
  for(j in 1:ncol(S)){
    S[i,j] = exp(-(del.early[i,j]/DuE)^2)
    S2[i,j] = exp(-(del.late[i,j]/DuL)^2)
    }
}

R = R2 = c()

for(i in 1:nrow(S)){
  R[i] = 1+sum(S[i,-i])
  R2[i] = 1+sum(S2[i,-i])
}

Wu = 1/R # uniqueness weights
Wu2 = 1/R2

#Dq = mean(ifelse(del.loc[,21]==0,NA,del.loc[,21]),na.rm=TRUE)/2
#Dq2 = mean(ifelse(del.locwrf[,21]==0,NA,del.locwrf[,21]),na.rm=TRUE)/2

#Wq = exp(-(ifelse(del.loc[,21]==0,NA,del.loc[,21])/Dq)^2) # quality weights
#Wq2 = exp(-(ifelse(del.locwrf[,21]==0,NA,del.locwrf[,21])/Dq2)^2)

W = Wu#*Wq # final model weights based on distances from EOFs between models and obs.
W2 = Wu2#*Wq2

#######################
# MDS prep

#######################
# setting random initializations

numsims = 20

inits = list()

for(i in 1:numsims){
init1 = matrix(runif(22,-1,1),nrow=11,ncol=2)
inits[[i]]=init1
}

####################################################
# MDS calcs

####################
# Early Period

mds.E = MDScalc(dist(del.early),initlist=inits)
coords.E = as.data.frame(do.call(cbind,mds.E))

####################
# Late Period

mds.L = MDScalc(dist(del.late),initlist=inits)
coords.L = as.data.frame(do.call(cbind,mds.L))

##################
# 



######################
# Weighted CDFs

month=8
load(paste("interpvalues_new_",month,".Rdata",sep=""))

library(Hmisc)

pravgsin = pravgs[c(1:20)]
pravgsin2 = pravgs[c(1:6,21,8:20)]
pravgsin3 = pravgs[c(1:6,23,8:20)]
pravgsin4 = pravgs[c(1:6,22,8:20)]
pravgsin5 = pravgs[c(1:6,24,8:20)]
pravgsobs = pravgs[25]

avgw.cdf = Ecdf(pravgsin,weights=W[1:20],pl=FALSE)
avgw2.cdf = Ecdf(pravgsin2,weights=W2[1:20],pl=FALSE)
avgw3.cdf = Ecdf(pravgsin3,weights=W3[1:20],pl=FALSE)
avgw4.cdf = Ecdf(pravgsin4,weights=W4[1:20],pl=FALSE)
avgw5.cdf = Ecdf(pravgsin5,weights=W5[1:20],pl=FALSE)

#save(list=c("avgw.cdf","avgw2.cdf","avgw3.cdf","avgw4.cdf","avgw5.cdf"),file=paste("weightedCDFS_",month,"_",locname,"_new.Rdata",sep=""))

pdf(paste("prweightedCDFs_",month,"_new.pdf",sep=""),width = 11, height=11,onefile=TRUE)
plot(avgw.cdf,xlim=c(0,7),col="blue",type="l",lwd=2,xlab="Precipitation (mm/day)",ylab="CDF",main="Weighted CDFs from Each Ensemble")
lines(avgw2.cdf,col="red",lwd=2)
lines(avgw3.cdf,col="darkred",lwd=2,lty=2)
lines(avgw4.cdf,col="green",lwd=2)
lines(avgw5.cdf,col="darkgreen",lwd=2,lty=2)
abline(h=0:1,lty=3)
abline(v=pravgsobs,col="black",lty=2,lwd=2)
legend("bottomright",c("GCMs w/ CCSM4","GCMs w/ CCSM4-WRF (2km)","GCMs w/ CCSM4-WRF (10km)","GCMs w/ CCSM4-NHM (2km)","GCMs w/ CCSM4-RSM (10km)","Obs."),col=c("blue","red","darkred","green","darkgreen","black"),lty=c(1,1,2,1,2,2),lwd=2,bg="white")
dev.off()

#####
# Wet season island average rainfall

pravgsin = pravgsin2 = pravgsin3 = pravgsin4 = pravgsin5 = matrix(NA,nrow=20,ncol=7)
pravgsobsv = c()

for(month in 4:10){
  load(paste("interpvalues_new_",month,".Rdata",sep=""))
  pravgsin[c(1:20),(month)-3] = pravgs[c(1:20)]
  pravgsin2[c(1:20),(month)-3] = pravgs[c(1:6,21,8:20)]
  pravgsin3[c(1:20),(month)-3] = pravgs[c(1:6,23,8:20)]
  pravgsin4[c(1:20),(month)-3] = pravgs[c(1:6,22,8:20)]
  pravgsin5[c(1:20),(month)-3] = pravgs[c(1:6,24,8:20)]
  pravgsobs[month-3] = c(pravgsobsv,pravgs[25])
}

#month=8

library(Hmisc)

avgw.cdf = Ecdf(apply(pravgsin,1,mean,na.rm=TRUE),weights=W[1:20],pl=FALSE)
avgw2.cdf = Ecdf(apply(pravgsin2,1,mean,na.rm=TRUE),weights=W2[1:20],pl=FALSE)
avgw3.cdf = Ecdf(apply(pravgsin3,1,mean,na.rm=TRUE),weights=W3[1:20],pl=FALSE)
avgw4.cdf = Ecdf(apply(pravgsin4,1,mean,na.rm=TRUE),weights=W4[1:20],pl=FALSE)
avgw5.cdf = Ecdf(apply(pravgsin5,1,mean,na.rm=TRUE),weights=W5[1:20],pl=FALSE)

#save(list=c("avgw.cdf","avgw2.cdf","avgw3.cdf","avgw4.cdf","avgw5.cdf"),file=paste("weightedCDFS_",month,"_",locname,"_new.Rdata",sep=""))

pdf(paste("prweightedCDFs_WET_new.pdf",sep=""),width = 11, height=11,onefile=TRUE)
plot(avgw.cdf,xlim=c(0,7),col="blue",type="l",lwd=2,xlab="Precipitation (mm/day)",ylab="CDF",main="Weighted CDFs from Each Ensemble")
lines(avgw2.cdf,col="red",lwd=2)
lines(avgw3.cdf,col="darkred",lwd=2,lty=2)
lines(avgw4.cdf,col="green",lwd=2)
lines(avgw5.cdf,col="darkgreen",lwd=2,lty=2)
abline(h=0:1,lty=3)
abline(v=mean(pravgsobs),col="black",lty=2,lwd=2)
legend("bottomright",c("GCMs w/ CCSM4","GCMs w/ CCSM4-WRF (2km)","GCMs w/ CCSM4-WRF (10km)","GCMs w/ CCSM4-NHM (2km)","GCMs w/ CCSM4-RSM (10km)","Obs."),col=c("blue","red","darkred","green","darkgreen","black"),lty=c(1,1,2,1,2,2),lwd=2,bg="white")
dev.off()




######################
# Distance changes distribution

write.table(intermoddist(coords.loc[,(ncol(coords.loc)-1):ncol(coords.loc)]),"intermodGCMloc_t3_n0obs.csv",sep=",",row.names=FALSE)
write.table(intermoddist(coords.loc2[,(ncol(coords.loc2)-1):ncol(coords.loc2)]),"intermodRCM2loc_t3_n0obs.csv",sep=",",row.names=FALSE)
write.table(intermoddist(coords.loc3[,(ncol(coords.loc3)-1):ncol(coords.loc3)]),"intermodRCM10loc_t3_n0obs.csv",sep=",",row.names=FALSE)

save(list = c("del.loc","del.loc2","del.loc3","coords.loc","coords.loc2","coords.loc3"),file="coords_and_del_loc_t3_n0obs.Rdata")
