
GCMpr = read.csv("/data/static_web/RCMES/GCMs_hist_pr.csv",sep=",",header=TRUE)
LOCApr = read.csv("/data/static_web/RCMES/LOCA_hist_pr.csv",sep=",",header=TRUE)
LIVNEHpr = read.csv("/data/static_web/RCMES/LIVNEH_hist_pr.csv",sep=",",header=TRUE)

errorLOCA = LOCApr
errorGCM = GCMpr

for(i in 1:ncol(LOCApr)){
  
  errorLOCA[,i] = ifelse(is.na(LOCApr[,i])==TRUE,0,LOCApr[,i])-ifelse(is.na(LIVNEHpr[,1])==TRUE,0,LIVNEHpr[,1])
  errorGCM[,i] = ifelse(is.na(GCMpr[,i])==TRUE,0,GCMpr[,i])-ifelse(is.na(LIVNEHpr[,1])==TRUE,0,LIVNEHpr[,1])
  
}

errsqLOCA = errorLOCA^2
meanerrLOCA = apply(errsqLOCA,2,mean,na.rm=TRUE)
RMSELOCA = sqrt(meanerrLOCA)

errsqGCM = errorGCM^2
meanerrGCM = apply(errsqGCM,2,mean,na.rm=TRUE)
RMSEGCM = sqrt(meanerrGCM)

sumerrLOCA = apply(errsqLOCA,2,sum,na.rm=TRUE)
RMSELOCA2 = sqrt(sumerrLOCA/nrow(errorLOCA))

